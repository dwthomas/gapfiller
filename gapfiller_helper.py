#! /usr/bin/env python
import argparse
import sys
from pathlib import Path
import pandas as pd
import geopandas as gpd
from pyproj import CRS
from shapely.geometry import LineString

import tempfile
import subprocess

import utils

wgs84 = CRS.from_epsg(4326)

def existing_dir(path_str: str) -> str:
    path = Path(path_str)
    if not path.is_dir():
        raise argparse.ArgumentTypeError(f"Directory does not exist: {path_str}")
    return str(path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Helper for gapfiller")
    parser.add_argument("--source-lat", type=str, required=True, help="Source latitude (any GeoPandas-compatible format)")
    parser.add_argument("--source-lon", type=str, required=True, help="Source longitude (any GeoPandas-compatible format)")
    parser.add_argument("--dest-lat", type=str, required=True, help="Destination latitude (any GeoPandas-compatible format)")
    parser.add_argument("--dest-lon", type=str, required=True, help="Destination longitude (any GeoPandas-compatible format)")
    parser.add_argument(
        "--budget",
        type=float,
        required=True,
        help="Budget in meters.",
    )
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
    parser.add_argument("--bin-path", type=str, default="src/release", help="Location of local_search")
    args = parser.parse_args()

    source_lat = float(args.source_lat)
    source_lon = float(args.source_lon)
    dest_lat = float(args.dest_lat)
    dest_lon = float(args.dest_lon)

    swath = args.swath

    command = "{bin_path}/local_search --unmapped {unmapped} --land {land} --dst_srs ESRI:54009 --budget {budget} --plan {plan}"

    line = LineString([(source_lon, source_lat), (dest_lon, dest_lat)])

    line_gdf = gpd.GeoDataFrame(
        geometry=[line],
        crs=wgs84,
    )
    metric_line = line_gdf.to_crs(utils.metric_crs)
    budget = float(args.budget) + float(metric_line.geometry[0].length)

    gebco_folder = args.gebco_dir

    envelope = utils.line_to_ellipse(line_gdf, width=budget, resolution = 4)  # Example width of 100 km+

    m = utils.Map(envelope, gebco_folder, extinction_file=args.extinction)
    with tempfile.TemporaryDirectory(delete=True) as tmpdir:
        unmapped_output_path = Path(tmpdir) / "unmapped_polygons.json"
        unmapped_output_path.write_text(m.unmapped_polygons.to_json())
        land_output_path = Path(tmpdir) / "land_polygons.json"
        land_output_path.write_text(m.land_polygons.to_json())
        plan_output_path = Path(tmpdir) / "plan.json"
        plan_output_path.write_text(line_gdf.to_json())
        cmd = command.format(
                bin_path=args.bin_path,
                unmapped=unmapped_output_path,
                land=land_output_path,
                budget=budget,
                plan=plan_output_path,
            )
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
        )
        output_gdf = gpd.GeoDataFrame(
            geometry=gpd.GeoSeries.from_wkt([result.stdout.strip()]),
            crs=wgs84,
        )
        if swath:
            swath_gdf = m.survey_line(output_gdf)
            output_gdf = gpd.GeoDataFrame(pd.concat([output_gdf, swath_gdf[0]], ignore_index=True), crs = output_gdf.crs)
        print(output_gdf.to_json())    
