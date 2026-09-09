#! /usr/bin/env python

import sys
import faulthandler
from pathlib import Path
import argparse

import geopandas as gpd

import cv2

from shapely import Polygon, LineString, Point
from shapely.geometry import box
from shapely.affinity import translate
from shapely.ops import linemerge, unary_union

from beam import utils
from pyproj import CRS, Transformer

import pandas as pd
import numpy as np
from scipy.ndimage import gaussian_filter
from matplotlib import pyplot as plt

from io import StringIO

# Enable faulthandler to print Python tracebacks on fatal errors (SIGSEGV, etc.)
faulthandler.enable(file=sys.stderr, all_threads=True)

def existing_dir(path_str: str) -> str:
    path = Path(path_str)
    if not path.is_dir():
        raise argparse.ArgumentTypeError(f"Directory does not exist: {path_str}")
    return str(path)

def largest_polygon(geom):
    if geom.geom_type == 'MultiPolygon':
        return max(geom.geoms, key=lambda p: p.area)
    return geom


def resample_linestring_keep_vertices(l: LineString, spacing: float) -> LineString:
    coords = list(l.coords)
    new_coords = [coords[0]]

    for (x0, y0), (x1, y1) in zip(coords[:-1], coords[1:]):
        seg_len = np.hypot(x1 - x0, y1 - y0)
        if seg_len == 0:
            continue
        # how many extra points to insert so spacing <= requested spacing
        n_segments = max(int(np.ceil(seg_len / spacing)), 1)
        # skip t=0 since that's the previous point already added
        for i in range(1, n_segments + 1):
            t = i / n_segments
            x = x0 + t * (x1 - x0)
            y = y0 + t * (y1 - y0)
            new_coords.append((x, y))

    return LineString(new_coords) 


def iter_segments(geom):
    if geom.geom_type == "LineString":
        coords = list(geom.coords)
        for start, end in zip(coords[:-1], coords[1:]):
            yield start, end
    elif geom.geom_type == "MultiLineString":
        for part in geom.geoms:
            yield from iter_segments(part)

def local_distance(p1, p2):
    return np.hypot(p2[0] - p1[0], p2[1] - p1[1])

def segment_normal(start, end, direction="left"):
    """
    Given a segment from start to end, return the unit vector perpendicular
    to the segment, pointing either 'left' or 'right' relative to the
    direction of travel from start -> end.

    'left' = counter-clockwise 90° rotation of the segment direction
    'right' = clockwise 90° rotation of the segment direction
    """
    x0, y0 = start[:2]
    x1, y1 = end[:2]

    dx, dy = x1 - x0, y1 - y0
    length = np.hypot(dx, dy)
    if length == 0:
        raise ValueError("start and end points are identical; direction is undefined")

    # unit vector along the segment
    ux, uy = dx / length, dy / length

    if direction == "left":
        # rotate 90° counter-clockwise: (x, y) -> (-y, x)
        nx, ny = -uy, ux
    elif direction == "right":
        # rotate 90° clockwise: (x, y) -> (y, -x)
        nx, ny = uy, -ux
    else:
        raise ValueError(f"direction must be 'left' or 'right', got {direction!r}")

    return (nx, ny)

def line_intersection(p1, p2, p3, p4):
    x1, y1 = p1[:2]
    x2, y2 = p2[:2]
    x3, y3 = p3[:2]
    x4, y4 = p4[:2]

    # translate everything relative to p1 to avoid catastrophic cancellation
    # when working with large absolute coordinates (e.g. EPSG:3857 meters)
    ox, oy = x1, y1
    x1, y1 = 0.0, 0.0
    x2, y2 = x2 - ox, y2 - oy
    x3, y3 = x3 - ox, y3 - oy
    x4, y4 = x4 - ox, y4 - oy

    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if abs(denom) < 1e-9:
        return None  # parallel (or coincident) lines

    px = ((x1*y2 - y1*x2) * (x3 - x4) - (x1 - x2) * (x3*y4 - y3*x4)) / denom
    py = ((x1*y2 - y1*x2) * (y3 - y4) - (y1 - y2) * (x3*y4 - y3*x4)) / denom
    return (px + ox, py + oy)


def connect_segments(segments):
    """
    Given a list of (start, end) segment pairs, replace each junction
    (end of segment i, start of segment i+1) with the intersection point
    of the two infinite lines, and return a single connected LineString.
    """
    coords = [segments[0][0]]
    
    for (start1, end1), (start2, end2) in zip(segments[:-1], segments[1:]):
        if local_distance(end1, start2) < 1000:
            pt = ((end1[0] + start2[0]) / 2, (end1[1] + start2[1]) / 2)
        else:
            pt = line_intersection(start1, end1, start2, end2)
        if pt is None:
            # parallel segments -- no unique intersection, fall back
            # to the midpoint of the gap between them
            pt = ((end1[0] + start2[0]) / 2, (end1[1] + start2[1]) / 2)
        coords.append(pt)

    coords.append(segments[-1][1])
    return LineString(coords)

def paired_interpolation(start, end, new_start, new_end, n=10):
    ts = np.linspace(0, 1, n)
    for t in ts:
        p1 = (start[0] + t * (end[0] - start[0]),
              start[1] + t * (end[1] - start[1]))
        p2 = (new_start[0] + t * (new_end[0] - new_start[0]),
              new_start[1] + t * (new_end[1] - new_start[1]))
        yield p1, p2

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Helper for gapfiller")
    parser.add_argument("--polygon", type=str, required=True, help="Polygon file (GeoJSON)")
    parser.add_argument("--tmpdir", type=str, default="/tmp", help="Temporary directory for intermediate files. Default: /tmp")
    parser.add_argument("--keep-tmp", action="store_true", help="Keep temporary files after execution. Default: False")
    parser.add_argument(
        "--gebco-dir",
        type=existing_dir,
        default="gebco_raster/",
        help="Path to folder containing the GEBCO dataset (must exist). Default: gebco_raster/",
    )
    parser.add_argument(
        "--extinction",
        type=str,
        default="EM302nautilus.txt",
        help="Extinction curve filename or comma-separated extinction curve. Default: EM302nautilus.txt\nExample: --extinction EM302nautilus.txt or --extinction 0.0 5.6,1608.0 6.6,3000.0 3.133,4000.0 2.205,5000.0 1.644,6000.0 1.198, ..."
    )
    parser.add_argument("--swath", action="store_true", help="Emit swath in addition to centerline.", default=False)
    parser.add_argument("--unmapped_file", type=str, help="Path to output unmapped raster file")
    parser.add_argument("--output-path", type=str, help="Path to output path file (GeoJSON)")
    parser.add_argument("--output-swath", type=str, help="Path to output swath file (GeoJSON)")
    args = parser.parse_args()

    transformer_localtowgs = Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)

    polygon_file = args.polygon

    swath = args.swath

    polygon = gpd.read_file(polygon_file)

    centroid = polygon.centroid
    lon, lat = centroid.x, centroid.y

    aeqd_crs = CRS.from_proj4(
        f"+proj=aeqd +lat_0={lat} +lon_0={lon} +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    )

    m = utils.Map(polygon.envelope, args.gebco_dir, extinction_file=args.extinction, unmapped_raster_path=args.unmapped_file)


    data = m.depth_raster
    top, left = 0, 0
    bottom, right = m.depth_raster.shape
    top_coord, left_coord = m.coords_of(0, 0)
    bottom_coord, right_coord = m.coords_of(m.depth_raster.shape[1]-1, m.depth_raster.shape[0]-1)
    nrows, ncols = data.shape
    x = np.linspace(left_coord, right_coord, ncols)
    y = np.linspace(bottom_coord, top_coord, nrows)
    X, Y = np.meshgrid(x, y)
    smoothed = gaussian_filter(data, sigma=9)

    levels = np.arange(-6000, 0, 25)

    # Matplotlib’s contour expects “height” data, so we invert sign for depth if needed
    cs = plt.contour(X, Y[::-1], smoothed, levels=levels, colors='black', linewidths=0.5)

    contours = cs
    lines = []
    lvl = []

    for level, segs in zip(contours.levels, contours.allsegs):
        for seg in segs:
            v = seg  # Nx2 array of (x, y) coordinates
            if len(v) > 2:
                lines.append(LineString(v))
                lvl.append(level)

    contours_gpd = gpd.GeoDataFrame(
        {"level": lvl, "geometry": lines},
        crs="EPSG:4326"
    )

    polygon2 = gpd.GeoDataFrame(geometry = polygon.to_crs(epsg=3857).buffer(2500), crs = "EPSG:3857")
    contours_gpd = gpd.overlay(contours_gpd, polygon2.to_crs("EPSG:4326"), "intersection")

    cbl = contours_gpd.to_crs("EPSG:3857").length.sort_values()[::-1]


    lines = contours_gpd.loc[[cbl.index[0]]].explode(index_parts = True)
    longest_ind = lines.length.sort_values()[::-1].index[0]

    initial_line = lines.loc[[longest_ind]].to_crs(epsg=3857)
    initial_line['geometry'] = initial_line.geometry.simplify(tolerance=2000, preserve_topology=True)

    swath = m.simple_survey_line(initial_line)
    swath = swath.to_crs("EPSG:3857")
    swath['geometry'] = swath.simplify(2000)
    swath = swath.to_crs("EPSG:4326")

    initial_resampled = gpd.GeoDataFrame(geometry = [resample_linestring_keep_vertices(initial_line.geometry.iloc[0], 10000)], crs = initial_line.crs).to_crs("EPSG:4326")

    # joined_line = []
    # for c1, c2 in zip(new_segments[:-1], new_segments[1:]):
    #     joined_line.append(c1[0])
    #     joined_line.append(c1[1])
    #     joined_line.append(c2[0])
    # joined_line.append(c2[1])
    # joined_line = gpd.GeoDataFrame(geometry = [LineString(joined_line)], crs = "EPSG:3857")
    # joined_line_simplified = joined_line.copy()
    # joined_line_simplified['geometry'] = joined_line_simplified.simplify(7500)

    # options_ls = gpd.GeoDataFrame(geometry = options_ls, crs = "EPSG:3857")
    desired_overlap = 0.15
    ev = 1-desired_overlap

    alternatives = np.linspace(0.25*ev, 0.75*ev, 7)

    plans = [initial_resampled]
    swaths = [swath]
    options_ls = []
    for direction in ["left", "right"]:
        current_plan = initial_resampled
        current_swath = swath
        acc = 0
        prev_remaining_area = 1.0
        while acc < 25:
            new_segments = []
            for ((_, row), (_, row_wgs84))in zip(current_plan.to_crs("EPSG:3857").iterrows(), current_plan.iterrows()):
                for ((start, end), (start_wgs_84, end_wgs84)) in zip(iter_segments(row.geometry), iter_segments(row_wgs84.geometry)):
                    segment = [start, end]
                    normal = segment_normal(start, end, direction)
                    midpoint = Point((end_wgs84[0] - start_wgs_84[0])/2 + start_wgs_84[0],(end_wgs84[1] - start_wgs_84[1])/2 + start_wgs_84[1])
                
                    tp = Point((end_wgs84[0] - start_wgs_84[0])/2 + start_wgs_84[0],(end_wgs84[1] - start_wgs_84[1])/2 + start_wgs_84[1])
                    width = m.width_at(tp)
                    scores = []
                    new_segs = []
                    for alternative_start, alternative_end in zip(alternatives, alternatives):
                        new_start = (start[0] + normal[0]*alternative_start*2*width, start[1] + normal[1]*alternative_start*2*width)
                        new_end = (end[0] + normal[0]*alternative_end*2*width, end[1] + normal[1]*alternative_end*2*width)
                        overlaps = []
                        for p1, p2 in paired_interpolation(start, end, new_start, new_end, n=10):
                            p1_wgs84 = transformer_localtowgs.transform(p1[0], p1[1]) 
                            p1_p = Point(p1_wgs84[0], p1_wgs84[1])
                            p1_width = m.width_at(p1_p)
                            p2_wgs84 = transformer_localtowgs.transform(p2[0], p2[1])
                            p2_p = Point(p2_wgs84[0], p2_wgs84[1])
                            p2_width = m.width_at(p2_p)
                            if p2_width <= 0:
                                continue
                            overlaps.append((p1_width/2 + p2_width/2 - local_distance(p2, p1))/ p2_width)
                        new_segs.append((new_start, new_end))
                        scores.append(np.linalg.norm(np.asarray(overlaps) - 0.15))
                    best_i = np.argmin(scores)
                    best_new_seg = new_segs[best_i]
                    new_segments.append((transformer_localtowgs.transform(best_new_seg[0][0], best_new_seg[0][1]) , transformer_localtowgs.transform(best_new_seg[1][0], best_new_seg[1][1])))
            connected_new_segments = connect_segments(new_segments)
            current_plan = gpd.GeoDataFrame(geometry = [connected_new_segments], crs = "EPSG:4326").to_crs("EPSG:3857")
            
            fixed_line = gpd.GeoDataFrame(geometry = [current_plan.union_all()], crs = current_plan.crs).explode()
            fixed_line = fixed_line.loc[~fixed_line.is_ring]
            try:
                current_plan['geometry'] =  gpd.GeoDataFrame(geometry = [linemerge(fixed_line.union_all())], crs = current_plan.crs)
            except:
                current_plan['geometry'] =  gpd.GeoDataFrame(geometry = [fixed_line.union_all()], crs = current_plan.crs)
            current_swath = m.simple_survey_line(current_plan)
            current_plan = current_plan.to_crs("EPSG:4326")
            plans.append(current_plan)
            swaths.append(current_swath.to_crs("EPSG:4326"))
            acc += 1
            print(direction, acc)    
            cds = gpd.GeoDataFrame( pd.concat(swaths))
            new_remaining = gpd.overlay(polygon.to_crs("EPSG:3857"), cds.to_crs("EPSG:3857"), "difference").area.values[0]/polygon.to_crs("EPSG:3857").area.values[0]
            print(prev_remaining_area, new_remaining)
            print("change in remaining area", prev_remaining_area - new_remaining)
            if prev_remaining_area - new_remaining < 0.001:
                prev_remaining_area = new_remaining
                break
            prev_remaining_area = new_remaining
            # print(segment, normal, width)
    cut_down_plans = []
    cut_down_swaths = []
    for plan in plans:
        pruned_plan = gpd.overlay(plan.to_crs("EPSG:3857"), polygon.to_crs("EPSG:3857"), "intersection")
        if len(pruned_plan.explode()) >1:
            pruned_plan['geometry'] = pruned_plan.explode().iloc[1:].union_all()
        # pruned_plan = pruned_plan.explode()
        if pruned_plan.to_crs("EPSG:3857").length.sum() > 10:
            cut_down_plans.append(pruned_plan)
            cut_down_swaths.append(m.simple_survey_line(pruned_plan.to_crs("EPSG:3857")))

    cds = gpd.GeoDataFrame( pd.concat(cut_down_swaths))
    print(cds.to_crs("epsg:4326").to_json(), file=open(args.output_swath, "w"))
    cdp = gpd.GeoDataFrame( pd.concat(cut_down_plans))
    print(cdp.to_crs("epsg:4326").to_json(), file=open(args.output_path, "w"))