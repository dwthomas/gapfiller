#! /usr/bin/env python
import argparse
import json
import antimeridian
import sys
import faulthandler
from pathlib import Path
import pandas as pd
import geopandas as gpd
from pyproj import CRS
from shapely.geometry import LineString

import tempfile
import subprocess

from beam import utils

# Enable faulthandler to print Python tracebacks on fatal errors (SIGSEGV, etc.)
faulthandler.enable(file=sys.stderr, all_threads=True)

def existing_dir(path_str: str) -> str:
    path = Path(path_str)
    if not path.is_dir():
        raise argparse.ArgumentTypeError(f"Directory does not exist: {path_str}")
    return str(path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Helper for gapfiller")
    parser.add_argument("--source-lat", type=str, required=True, help="Source latitude (any GeoPandas-compatible format)")
    parser.add_argument("--source-lon", type=str, required=True, help="Source longitude (any GeoPandas-compatible format)")
    parser.add_argument("--dest-lat", "--end-lat", type=str, required=True, help="Destination latitude (any GeoPandas-compatible format)")
    parser.add_argument("--dest-lon", "--end-lon", type=str, required=True, help="Destination longitude (any GeoPandas-compatible format)")
    parser.add_argument("--tmpdir", type=str, default="/tmp", help="Temporary directory for intermediate files. Default: /tmp")
    parser.add_argument("--keep-tmp", action="store_true", help="Keep temporary files after execution. Default: False")
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
    parser.add_argument("--unmapped_file", type=str, help="Path to output unmapped raster file")
    args = parser.parse_args()

    source_lat = float(args.source_lat)
    source_lon = float(args.source_lon)
    dest_lat = float(args.dest_lat)
    dest_lon = float(args.dest_lon)

    swath = args.swath

    command = "{bin_path}/local_search --unmapped {unmapped} --land {land} --budget {budget} --plan {plan}"

    line = LineString([(source_lon, source_lat), (dest_lon, dest_lat)])

    line_gdf = gpd.GeoDataFrame(
        geometry=[line],
        crs=utils.wgs84,
    )
    centroid = line_gdf.geometry.iloc[0].centroid
    metric_crs = utils.metric_crs(centroid.x)
    # print("initial length:", line_gdf.to_crs( metric_crs).length[0])
    # print("budget:", args.budget, "total:", float(args.budget) + line_gdf.to_crs( metric_crs).length[0])
    budget = float(args.budget) + line_gdf.to_crs( metric_crs).length[0]
    metadata = {
        "initial_length_m": line_gdf.to_crs( metric_crs).length[0],
        "budget_m": float(args.budget),
    }
    gebco_folder = args.gebco_dir
    print("Calculating ellipse...", file=sys.stderr)
    envelope = utils.line_to_ellipse(line_gdf, width=float(args.budget), resolution = 4)  # Example width of 100 km+
    # print("envelope", envelope.to_crs(utils.wgs84).to_json())
    m = utils.Map(envelope, gebco_folder, extinction_file=args.extinction, unmapped_raster_path=args.unmapped_file)
    with tempfile.TemporaryDirectory(delete=not args.keep_tmp, dir=args.tmpdir) as tmpdir:
        # print(tmpdir)
        unmapped_output_path = Path(tmpdir) / "unmapped_polygons.json"
        unmapped_output_path.write_text(m.unmapped_polygons.to_json())
        # print(m.unmapped_polygons.to_json(), flush=True)
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
        # print(cmd)
        # Run the command, piping only stderr (capture stderr, let stdout go to the console)
        print("Searching for plans...", file=sys.stderr)
        result = subprocess.run(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        print("Processing Plan...", file=sys.stderr)
        # Prefer stderr if the external tool writes WKT there, otherwise use stdout.
        output_wkt = result.stdout.strip()
        # print("C++ stderr:", result.stderr.strip(), flush=True)
        output_gdf = gpd.GeoDataFrame(
            geometry=gpd.GeoSeries.from_wkt([output_wkt]),
            crs=utils.wgs84,
        )
        # print("output length:", output_gdf.to_crs( metric_crs).length[0])
        output_gdf = output_gdf.to_crs(metric_crs)
        output_gdf['geometry'] = output_gdf['geometry'].simplify(2000, preserve_topology=True)
        output_gdf = output_gdf.to_crs(utils.wgs84)
        geojson_dict = json.loads(output_gdf.to_json())
        fixed = antimeridian.fix_geojson(geojson_dict)
        output_gdf = gpd.GeoDataFrame.from_features(fixed["features"], crs=utils.wgs84)
        metadata["final_length_m"] = output_gdf.to_crs(metric_crs).length[0]
        metadata["remaining budget_m"] = metadata["initial_length_m"] + metadata["budget_m"] - metadata["final_length_m"]
        if swath:
            swath_gdf = m.simple_survey_line(output_gdf.to_crs(metric_crs))
            swath_gdf = swath_gdf.to_crs( metric_crs)
                        # print(output_gdf.area[1]/output_gdf.length[0])
            buf = 1000
            swath_gdf['geometry'] = swath_gdf['geometry'].buffer(buf)
            swath_gdf['geometry'] = swath_gdf['geometry'].simplify(buf, preserve_topology=True)
            swath_gdf['geometry'] = swath_gdf['geometry'].union_all()
            swath_gdf['geometry'] = swath_gdf['geometry'].buffer(-buf)
            # print(output_gdf.crs, swath_gdf.crs)
            output_gdf = gpd.GeoDataFrame(pd.concat([output_gdf, swath_gdf.to_crs(utils.wgs84)], ignore_index=True), crs =  utils.wgs84)

            initial_area = m.unmapped_polygons.to_crs(metric_crs).area.sum()
            new_unmapped = m.unmapped_polygons.to_crs(metric_crs).overlay(swath_gdf.to_crs(metric_crs), how='difference')
            metadata["newly_unmapped_area_m2"] = initial_area - new_unmapped.area.sum()

        geojson_dict = json.loads(output_gdf.to_json())
        fixed = antimeridian.fix_geojson(geojson_dict)
        output_gdf = gpd.GeoDataFrame.from_features(fixed["features"], crs=utils.wgs84)
       
        # fixed = geojson_dict
        geojson_dict = json.loads(output_gdf.to_json())
        geojson_dict["properties"] = metadata
        print(json.dumps(geojson_dict))