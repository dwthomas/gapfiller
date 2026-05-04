import numpy as np
import geopandas as gpd
import shapely
import os
from pyproj import Geod, CRS
from shapely.geometry import Point, LineString, Polygon, MultiLineString, MultiPolygon, shape, box, mapping
from shapely.ops import unary_union, snap, linemerge
import time
import subprocess
import glob
import rasterio
from rasterio.plot import reshape_as_image
from rasterio.windows import from_bounds
from rasterio.features import shapes
from rasterio.transform import xy, rowcol
from pyproj import Transformer
from rasterio.merge import merge
from rasterio.mask import mask as rio_mask

import pandas as pd


wgs84 = CRS.from_epsg(4326)
web_mercator = CRS.from_epsg(3857)
metric_crs = web_mercator

land_buffer_width = 3000

def remove_holes(polygon):
    if polygon.is_valid:
        # Recreate the polygon with no holes (interiors)
        return Polygon(polygon.exterior)
    return polygon

def shrink_polygon(polygon, distance):
    """Shrinks a polygon by offsetting its exterior ring inward.
    
    - Uses `unary_union` to merge fragmented pieces after shrinking.
    - Handles MultiLineString by constructing a valid MultiPolygon.
    - Preserves interior holes where possible.
    """
    
    # Offset the exterior ring inward
    offset_ext = polygon.exterior.offset_curve(-distance)

    # If offsetting results in a MultiLineString, try to create multiple polygons
    if isinstance(offset_ext, MultiLineString):
        polygons = [Polygon(line) for line in offset_ext.geoms if line.is_ring]
    else:
        polygons = [Polygon(offset_ext)] if offset_ext.is_ring else []

    # Preserve interior holes from the original polygon
    for poly in polygons:
        poly = Polygon(poly.exterior, holes=[hole for hole in polygon.interiors if poly.contains(hole)])

    # Merge all resulting polygons into a single geometry
    result = unary_union(polygons)

    return result if not result.is_empty else polygon  # Return original if shrinking failed


def p2vgp(p):
    return vg.Point(p[0], p[1])

def isend(edge, point):
    return (edge.p1 == point) or (edge.p2 == point)
sea_wkt = """PROJCS["ProjWiz_Custom_Equidistant_Cylindrical",
 GEOGCS["GCS_WGS_1984",
  DATUM["D_WGS_1984",
   SPHEROID["WGS_1984",6378137.0,298.257223563]],
  PRIMEM["Greenwich",0.0],
  UNIT["Degree",0.0174532925199433]],
 PROJECTION["Equidistant_Cylindrical"],
 PARAMETER["False_Easting",0.0],
 PARAMETER["False_Northing",0.0],
 PARAMETER["Central_Meridian",130.078125],
 PARAMETER["Standard_Parallel_1",25.8090083],
 UNIT["Meter",1.0]]"""
sea_crs = CRS.from_wkt(sea_wkt)

#intersects = guc[guc.intersects(gdf_path_line.geometry.iloc[0]).any()]
def sci_utility(intersections):
    return intersections.to_crs('ESRI:54009').length.sum()/1000

def length(gdf_path_line):
    return gdf_path_line.to_crs('ESRI:54009').length.sum()/1000

def get_verts(point_on_edge, polygon):
    for i in range(len(polygon.exterior.coords) - 1):
        edge = LineString([polygon.exterior.coords[i], polygon.exterior.coords[i + 1]])
        if edge.distance(point_on_edge) < 10:
            return polygon.exterior.coords[i], polygon.exterior.coords[i + 1]
    return None

def get_edges(point_on_edge, polygon, stop_point):
    start = None
    for i in range(len(polygon.exterior.coords) - 1):
        edge = LineString([polygon.exterior.coords[i], polygon.exterior.coords[i + 1]])
        if edge.distance(point_on_edge) < 10:
            start = i
            break #= polygon.exterior.coords[i], polygon.exterior.coords[i + 1]
    forward_edges = []
    for i in range(len(polygon.exterior.coords) - 1):
        index = (start + i) % len(polygon.exterior.coords)
        next_i = (start + i + 1) % len(polygon.exterior.coords)
        edge = LineString([polygon.exterior.coords[index], polygon.exterior.coords[next_i]])
        if i == 0:
            edge = LineString([point_on_edge, polygon.exterior.coords[next_i]])
        if edge.distance(stop_point) < 10:
            forward_edges.append(LineString([polygon.exterior.coords[index], stop_point]))
            break
        else:
            forward_edges.append(edge)
    reverse_edges = []
    for i in range(len(polygon.exterior.coords) - 1, 0, -1):
        index = (start + i) % len(polygon.exterior.coords)
        next_i = (start + i - 1) % len(polygon.exterior.coords)
        edge = LineString([polygon.exterior.coords[index], polygon.exterior.coords[next_i]])
        if edge.distance(stop_point) < 10:
            reverse_edges.append(LineString([polygon.exterior.coords[index], stop_point]))
            break
        else:
            reverse_edges.append(edge)
    return forward_edges, reverse_edges

def combine_almost_continuous_lines(multi_line, tolerance=10):
    if isinstance(multi_line, LineString):
        return multi_line  # Already a LineString

    if not isinstance(multi_line, MultiLineString):
        raise ValueError("Input must be a LineString or MultiLineString")

    # Step 1: Snap lines to close small gaps
    snapped = snap(multi_line, multi_line, tolerance)

    # print(snapped)

    # Step 2: Merge snapped lines
    merged = unary_union(snapped)

    # print(merged)
    
    # Step 3: Ensure we produced a single LineString
    if isinstance(merged, LineString):
        return merged

    # Step 4: Attempt to order and merge segments
    if isinstance(merged, MultiLineString):
        lines = list(merged.geoms)
        ordered_lines = [lines.pop(0)]

        while lines:
            current = ordered_lines[-1]
            for i, line in enumerate(lines):
                # Check if the current line connects to any other line
                if current.coords[-1] == line.coords[0]:
                    ordered_lines.append(lines.pop(i))
                    break
                elif current.coords[-1] == line.coords[-1]:
                    ordered_lines.append(LineString(line.coords[::-1]))
                    lines.pop(i)
                    break
                elif current.coords[0] == line.coords[-1]:
                    ordered_lines.insert(0, lines.pop(i))
                    break
                elif current.coords[0] == line.coords[0]:
                    ordered_lines.insert(0, LineString(line.coords[::-1]))
                    lines.pop(i)
                    break
            else:
                break  # No more connections found

        # Combine ordered lines if all are connected
        if len(lines) == 0:
            return LineString([pt for line in ordered_lines for pt in line.coords])

    raise ValueError("Cannot combine lines into a single continuous LineString")

def get_pos(lat, lng):
    return lat, lng

def line_to_ellipse(line, width, metric_crs = web_mercator, resolution=64):
    """
    Constructs a GeoDataFrame of ellipses where each ellipse has the two points 
    from each segment of the line as the foci, and the sum of distances that 
    defines the ellipse is equal to the width. The ellipses are rotated to align 
    with the line segments.

    Parameters:
        line (LineString): A LineString representing the line.
        width (float): The sum of distances (major axis length) for the ellipses.

    Returns:
        GeoDataFrame: A GeoDataFrame containing the ellipses as geometries.
    """
    ellipses = []

    metric_line = line.to_crs(metric_crs)

    # Iterate through each segment of the line
    for i in range(len(metric_line.geometry[0].coords) - 1):
        p1 = Point(metric_line.geometry[0].coords[i])
        p2 = Point(metric_line.geometry[0].coords[i + 1])

        # Calculate the center of the ellipse
        center = Point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2)

        # Calculate the distance between the two foci
        foci_distance = p1.distance(p2)

        # Calculate the semi-major axis (half of the width)
        semi_major = (width + foci_distance)/ 2        

        # Calculate the semi-minor axis using the ellipse formula
        semi_minor = (semi_major**2 - (foci_distance / 2)**2)**0.5

        # Create a unit circle and scale it to the ellipse dimensions
        ellipse = shapely.affinity.scale(center.buffer(1, resolution=resolution), xfact=semi_major, yfact=semi_minor)

        # Calculate the angle of rotation (in degrees) to align the ellipse with the line segment
        angle = np.degrees(np.arctan2(p2.y - p1.y, p2.x - p1.x))

        # Rotate the ellipse around its center
        ellipse = shapely.affinity.rotate(ellipse, angle, origin='center')

        # Append the ellipse to the list
        ellipses.append(ellipse)
    # Create a GeoDataFrame from the ellipses
    gdf = gpd.GeoDataFrame(geometry=ellipses, crs=metric_line.crs)

    return gdf.to_crs(line.crs)

def load_gebco_region(tile_paths: list[str], polygon):
    """
    Args:
        tile_paths: List of paths to all GEBCO .tif files
        polygon: A Shapely geometry (Polygon or MultiPolygon) in WGS84
    """
    # print(type(polygon))
    polygon = polygon.geometry[0]
    bbox = polygon.bounds  # (min_lon, min_lat, max_lon, max_lat)
    min_lon, min_lat, max_lon, max_lat = bbox

    # Find overlapping tiles using bounding box (fast)
    overlapping = []
    for path in tile_paths:
        with rasterio.open(path) as src:
            b = src.bounds
            if b.left < max_lon and b.right > min_lon and b.bottom < max_lat and b.top > min_lat:
                overlapping.append(path)

    if not overlapping:
        raise ValueError("No GEBCO tiles overlap the requested polygon.")

    # print(f"Found {len(overlapping)} overlapping tile(s).")

    # Convert polygon to GeoJSON-like dict for rasterio
    geom = [mapping(polygon)]

    if len(overlapping) == 1:
        with rasterio.open(overlapping[0]) as src:
            data, transform = rio_mask(src, geom, crop=True, nodata=-1)
            crs = src.crs
    else:
        # Mosaic tiles first, then mask to polygon
        datasets = [rasterio.open(p) for p in overlapping]
        try:
            mosaic, mosaic_transform = merge(
                datasets,
                bounds=(min_lon, min_lat, max_lon, max_lat),
            )
            crs = datasets[0].crs
        finally:
            for ds in datasets:
                ds.close()

        # Write mosaic to a memory file, then apply polygon mask
        from rasterio.io import MemoryFile
        with MemoryFile() as memfile:
            with memfile.open(
                driver="GTiff",
                height=mosaic.shape[1],
                width=mosaic.shape[2],
                count=1,
                dtype=mosaic.dtype,
                crs=crs,
                transform=mosaic_transform,
            ) as mem_ds:
                mem_ds.write(mosaic)
                data, transform = rio_mask(mem_ds, geom, crop=True, nodata=-1)

    return data[0], transform, crs 


def load_raster(filename, bbox):
    with rasterio.open(filename) as ds:
        window = from_bounds(*bbox, transform=ds.transform)
        crs = ds.crs
        transform = rasterio.windows.transform(window, ds.transform)
        return ds.read(1, window=window), transform, crs

class Map:
    def __init__(self, mask, gebco_folder, extinction_file = "EM302nautilus.txt"):
        self.beam = None
        if os.path.isfile(extinction_file):
            self.beam = pd.read_csv(extinction_file, sep = r"\s+", header=None, names= ["depth", "extinction"], index_col = "depth")
        else:
            extinction = []
            for pair in extinction_file.split(","):
                depth, ext = pair.strip().split()
                extinction.append((float(depth), float(ext)))
            self.beam = pd.DataFrame(extinction, columns=["depth", "extinction"]).set_index("depth")
        
        tid_files = glob.glob(os.path.join(gebco_folder, "*_tid_*.tif"))
        self.tid_raster, self.tid_transform, self.tid_crs = load_gebco_region(tid_files, mask)
        depth_files = glob.glob(os.path.join(gebco_folder, "*_sub_ice_*.tif"))
        self.depth_raster, self.depth_transform, self.depth_crs = load_gebco_region(depth_files, mask)
        # self.land_raster = (255 * (self.depth_raster < 0)).astype(np.uint8)
        self.land_raster = (self.tid_raster == 0).astype(np.uint8)
        self.tid_raster = self.tid_raster.astype(np.uint8)
        self.unmapped_raster = (self.tid_raster != 11).astype(np.uint8)  # (self.tid_raster > 17) * (1 - self.land_raster)

        # print(self.depth_raster.shape)
        # bbox = (-180.0, 0, 180, 90)
        
        # left, bottom, right, top = bbox
        # Polygonize land        
        land_polygons = []
        for geom, value in shapes(
            self.land_raster,
            mask=self.land_raster,
            transform=self.tid_transform
        ):
            if value == 1:
                land_polygons.append(shape(geom))
        # Create GeoDataFrame
        self.land_polygons = gpd.GeoDataFrame(geometry=land_polygons, crs=self.tid_crs)
        self.grow_land_polygons()

        #polygonize unmapped
        unmapped_polygons = []
        for geom, value in shapes(
            self.unmapped_raster,
            mask=self.unmapped_raster,
            transform=self.tid_transform
        ):
            if value == 1:
                unmapped_polygons.append(shape(geom))
        # Create GeoDataFrame
        self.unmapped_polygons  = gpd.GeoDataFrame(geometry=unmapped_polygons, crs=self.tid_crs)
        merged = self.unmapped_polygons.geometry.union_all()
        self.unmapped_polygons = gpd.GeoDataFrame(
            geometry=[merged], crs=self.unmapped_polygons.crs
        ).explode(index_parts=False, ignore_index=True)
        self.shrink_unmapped_polygons()
        


    def plot_mapped_raster(self):
        plt.figure(figsize=(10, 5))
        plt.imshow(self.mapped_raster, cmap="Blues_r")
        plt.colorbar(label="TID")
        plt.title("GEBCO 2025 Bathymetry")
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")

    def plot_min_raster(self):
        plt.figure(figsize=(10, 5))
        plt.imshow(self.min_depth_raster, cmap="Blues_r")
        plt.colorbar(label="Elevation")
        plt.title("GEBCO 2008 Bathymetry")
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")


    def plot_land(self, m = None):
        if not m:
            m = folium.Map(location=[lat1, lon1],  zoom_start=8, tiles="Esri.OceanBasemap")

        # Add land polygons
        folium.GeoJson(
            self.land_polygons,
            name="Land",
            style_function=lambda x: {
                "fillColor": "#2ecc71",   # green land
                "color": "#145a32",       # border
                "weight": 1,
                "fillOpacity": 0.6,
            },
        ).add_to(m)
        
        folium.LayerControl().add_to(m)
        
        return m


        

    def plot_unmapped(self, m = None):
        if not m:
            m = folium.Map(location=[lat1, lon1], zoom_start=8, tiles="Esri.OceanBasemap")
        folium.GeoJson(
            self.unmapped_polygons,
            style_function=lambda feature: {
                "fillColor": feature["properties"]["color"],
                "color": feature["properties"]["color"],
                "weight": 1,
                "fillOpacity": 0.6,
            },
        ).add_to(m)
        return m

    def shrink_unmapped_polygons(self):
        shrunk = []
        for pgon in self.unmapped_polygons.geometry:
            centroid = pgon.centroid
            beam_width = self.width_at(centroid)
            pgon_merc = gpd.GeoSeries([pgon], crs=self.unmapped_polygons.crs)
            pgon_merc = pgon_merc.simplify(.005)

            pgon_merc = pgon_merc.to_crs(web_mercator).iloc[0]
            new_gon_merc = pgon_merc.buffer(-beam_width / 2, join_style=3)
            if new_gon_merc.area > 0:
                new_gon = gpd.GeoSeries([new_gon_merc], crs=web_mercator).to_crs(self.unmapped_polygons.crs).iloc[0]
                shrunk.append(new_gon)

        self.unmapped_polygons = gpd.GeoDataFrame(geometry=shrunk, crs=self.unmapped_polygons.crs)

    def grow_land_polygons(self):
        grown = []
        for pgon in self.land_polygons.geometry:
            centroid = pgon.centroid
            
            pgon_merc = gpd.GeoSeries([pgon], crs=self.land_polygons.crs)

            pgon_merc = pgon_merc.to_crs(web_mercator).iloc[0]
            new_gon_merc = pgon_merc.buffer(land_buffer_width, join_style=3)
            if new_gon_merc.area > 0:
                new_gon = gpd.GeoSeries([new_gon_merc], crs=web_mercator).to_crs(self.land_polygons.crs).iloc[0]
                new_gon = new_gon.simplify(.005)
                grown.append(new_gon)

        self.land_polygons = gpd.GeoDataFrame(geometry=grown, crs=self.land_polygons.crs)


    def index_of(self, point):
        transformer = Transformer.from_crs(
            "EPSG:4326",
            self.depth_crs,
            always_xy=True
        )

        x, y = transformer.transform(point.x, point.y)

        # Convert to row/col in window
        row, col = rowcol(self.depth_transform, x, y)
        return col, row

    def in_radius_of(self, point, radius):
        col, row = self.index_of(point)
        
        transformer = Transformer.from_crs("EPSG:4326", "EPSG:32634", always_xy=True)
        
        x1, y1 = rasterio.transform.xy(self.depth_transform, row, col)
        x1_m, y1_m = transformer.transform(x1, y1)

        points = [(col, row)]
        for i in [1, -1]: # up/down
            for j in [1, -1]: # left/right
                y = 1
                while True:
                    r = row + i*y
                    gotone = False
                    x = 1
                    while True:
                        c = col + j*x
                        x2, y2 = rasterio.transform.xy(self.depth_transform, r, c)
                        x2_m, y2_m = transformer.transform(x2, y2)

                        dist = np.sqrt((x2_m - x1_m)**2 + (y2_m - y1_m)**2)
                        if dist > radius:
                            break
                        gotone = True
                        points.append((c, r))
                        x+=1
                    if not gotone:
                        break
                    y += 1
        return points


        

    def depth_at(self, point):
        col, row = self.index_of(point)
        # print(row, col, self.depth_raster.shape)
        if row >= self.depth_raster.shape[0]:
            row = self.depth_raster.shape[0]-1
        if col >= self.depth_raster.shape[1]:
            col = self.depth_raster.shape[1]-1
        # print(col, row)
        return self.depth_raster[row, col]

    def width_at(self, point):
        s = self.beam["extinction"]
        depth = -self.depth_at(point)
        s2 = s.reindex(s.index.union([depth])).sort_index()   
        extinction = s2.interpolate(method='index').loc[depth]
        result = depth * extinction
        if not np.isfinite(result):
            return 0.0
            raise ValueError(f"Non-finite width calculated at point {point}. Depth: {depth}, Extinction: {extinction}")
        return result

    def survey_line(self, line):
        gdf_m = line.to_crs(metric_crs)
        line_m = gdf_m.geometry.iloc[0]
        to_wgs84 = Transformer.from_crs(metric_crs, "EPSG:4326", always_xy=True)

        step = 1000  # meters
        length = line_m.length
        distances = np.arange(0, length + step, step)

        def unit(v):
            return v / np.linalg.norm(v)
        
        left_pts = []
        right_pts = []
        
        for s in distances:
            p = line_m.interpolate(s)
        
            # tangent via finite difference
            eps = 1.0
            s1 = max(s - eps, 0)
            s2 = min(s + eps, length)
        
            p1 = line_m.interpolate(s1)
            p2 = line_m.interpolate(s2)
            # print(p, p1, p2)
            t = np.array([p2.x - p1.x, p2.y - p1.y])
            try:
                t = unit(t)
            except:
                continue
        
            # perpendicular normal
            n = np.array([-t[1], t[0]])
        
            # convert sample point → WGS84 for width()
            lon, lat = to_wgs84.transform(p.x, p.y)
            p_wgs84 = Point(lon, lat)
            try:
                w = self.width_at(p_wgs84) / 2.0   # meters
            except:
                continue
            left = Point(p.x + n[0]*w, p.y + n[1]*w)
            right = Point(p.x - n[0]*w, p.y - n[1]*w)

            if np.isfinite(left.x) and np.isfinite(left.y):
                left_pts.append(left)
            if np.isfinite(right.x) and np.isfinite(right.y):
                right_pts.append(right)
        
    
            # print(left, right)

        poly_m = Polygon(list(left_pts) + list(reversed(right_pts)))
        poly_wgs84 = gpd.GeoSeries([poly_m], crs=metric_crs).to_crs("EPSG:4326").iloc[0]
        poly_gdf = gpd.GeoDataFrame(geometry=[poly_wgs84], crs="EPSG:4326")

        lefts = gpd.GeoDataFrame(geometry = [LineString(list(left_pts))], crs = metric_crs).to_crs("EPSG:4326")
        rights = gpd.GeoDataFrame(geometry = [LineString(list(right_pts))], crs = metric_crs).to_crs("EPSG:4326")
        return poly_gdf, lefts, rights

def get_polys():
    if os.path.exists("glc_simp.feather") and os.path.exists("guc_simp.feather") and os.path.exists("glc_orig.feather"):
        glc_simp = gpd.read_feather("glc_simp.feather")
        guc_simp = gpd.read_feather("guc_simp.feather")
        glc_orig = gpd.read_feather("glc_orig.feather")
        return glc_simp, guc_simp, glc_orig, get_center(glc_simp.to_crs(epsg=4326))
    guc_simp = gpd.read_file("GebcoHICrop/guc.json")
    guc_simp = guc_simp.to_crs('ESRI:54009')
    ind = np.argsort(-guc_simp.geometry.area)
    x = guc_simp.iloc[ind]
    y_guc = guc_simp.iloc[ind]
    x.geometry = x.geometry.apply(remove_holes)
    x = x.buffer(-2000)
    guc_simp.geometry = x.simplify(1000)
    guc_simp = guc_simp[~guc_simp.geometry.is_empty]


    glc_orig = gpd.read_file("GebcoHICrop/merge-glc.json")
    glc_orig = glc_orig.to_crs('ESRI:54009')
    glc_simp = glc_orig.copy()
    ind = np.argsort(-glc_simp.geometry.area)
    x = glc_simp.iloc[ind]
    y = glc_simp.iloc[ind]
    x.geometry = x.geometry.apply(remove_holes)
    x.geometry = x.buffer(1000)
    x.geometry = x.simplify(1000)
    x.geometry = x.geometry.apply(lambda y: y.convex_hull)
    x = x.dissolve().explode(index_parts=False).reset_index(drop=True)
    x.geometry = x.simplify(1000)

    glc_simp = x
    glc_simp['geometry'] = glc_simp['geometry'].apply(remove_holes)
    glc_simp = glc_simp.loc[glc_simp.geometry.area > 50000000] # get rid of small areas
    glc_simp.to_feather("glc_simp.feather")
    glc_orig.to_feather("glc_orig.feather")
    guc_simp.to_feather("guc_simp.feather")
    return glc_simp.to_crs(epsg=4326), guc_simp.to_crs(epsg=4326), glc_orig.to_crs(epsg=4326), get_center(glc_simp.to_crs(epsg=4326))

def get_center(_glc_simp):
    center = _glc_simp.union_all().centroid
    return (np.asarray(center.coords[0])[::-1])
