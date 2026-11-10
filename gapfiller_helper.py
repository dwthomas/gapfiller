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
        "--gap-bonus-km",
        type=float,
        default=10.0,
        help="Gap completion bonus in km of polygon-equivalent mapping length. Default: 10.",
    )
    parser.add_argument(
        "--consider-gaps",
        type=int,
        default=5,
        help="Maximum number of nearest unmapped lines to consider per iteration. Default: 5.",
    )
    parser.add_argument(
        "--n-agents",
        type=int,
        default=1,
        help="Number of agents. Default: 1.",
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

    command = (
        "{bin_path}/local_search --unmapped-points {unmapped_points} "
        "--unmapped-lines {unmapped_lines} --unmapped-polygons {unmapped_polygons} "
        "--land {land} --budget {budget} --gap-bonus-km {gap_bonus_km} "
        "--consider-gaps {consider_gaps} --plan {plan}"
    )

    line = LineString([(source_lon, source_lat), (dest_lon, dest_lat)])

    line_gdf = gpd.GeoDataFrame(
        geometry=[line],
        crs=utils.wgs84,
    )

    centroid = line_gdf.geometry.iloc[0].centroid
    metric_crs = utils.metric_crs(centroid.x)
    # print("initial length:", line_gdf.to_crs( metric_crs).length[0])
    # print("budget:", args.budget, "total:", float(args.budget) + line_gdf.to_crs( metric_crs).length[0])
    budget = float(args.budget) + line_gdf.to_crs(metric_crs).length[0]
    
    metadata = {
        "initial_length_m": line_gdf.to_crs( metric_crs).length[0],
        "budget_m": float(args.budget),
    }

    gebco_folder = args.gebco_dir
    print("Calculating ellipse...", file=sys.stderr)
    envelope = utils.line_to_ellipse(line_gdf, width=float(args.budget), resolution = 4)  # Example width of 100 km+
    # print("envelope", envelope.to_crs(utils.wgs84).to_json())
    m = utils.Map(envelope, gebco_folder, extinction_file=args.extinction, unmapped_raster_path=args.unmapped_file, n_agents = args.n_agents)
    with tempfile.TemporaryDirectory(delete=not args.keep_tmp, dir=args.tmpdir) as tmpdir:
        # print(tmpdir)
        unmapped_output_path = Path(tmpdir) / "unmapped_polygons.json"
        unmapped_json = json.loads(m.unmapped_polygons.to_json())
        unmapped_json["properties"] = {"unmapped_scores": m.unmapped_scores.values.tolist()}
        unmapped_output_path.write_text(json.dumps(unmapped_json))

        unmapped_lines_output_path = Path(tmpdir) / "unmapped_lines.json"
        unmapped_lines_json = json.loads(m.unmapped_lines.to_json())
        unmapped_lines_json["properties"] = {
            "unmapped_scores": [
                m.width_at(line.centroid)
                for line in m.unmapped_lines.geometry
            ]
        }
        unmapped_lines_output_path.write_text(json.dumps(unmapped_lines_json))

        unmapped_points_output_path = Path(tmpdir) / "unmapped_points.json"
        unmapped_points_json = json.loads(m.unmapped_points.to_json())
        unmapped_points_json["properties"] = {
            "unmapped_scores": [
                m.width_at(point)
                for point in m.unmapped_points.geometry
            ]
        }
        unmapped_points_output_path.write_text(json.dumps(unmapped_points_json))

        # print(m.unmapped_polygons.to_json(), flush=True)
        land_output_path = Path(tmpdir) / "land_polygons.json"
        land_output_path.write_text(m.land_polygons.to_json())
        plan_output_path = Path(tmpdir) / "plan.json"
        plan_output_path.write_text(line_gdf.to_json())
        cmd = command.format(
            bin_path=args.bin_path,
            unmapped_points=unmapped_points_output_path,
            unmapped_lines=unmapped_lines_output_path,
            unmapped_polygons=unmapped_output_path,
            land=land_output_path,
            budget=budget,
            gap_bonus_km=args.gap_bonus_km,
            consider_gaps=args.consider_gaps,
            plan=plan_output_path,
        )
        print(cmd)
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
        # print("Output WKT:", output_wkt, flush=True)
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
        centerline_gdf = output_gdf
        if args.n_agents == 2:
            split_plan = utils.split_plan_between_agents(m, centerline_gdf)
            output_gdf = gpd.GeoDataFrame.from_features(
                split_plan["features"], crs=utils.wgs84
            )
        if swath:
            if args.n_agents == 2:
                agent_swaths = [
                    m.simple_survey_line(
                        gpd.GeoDataFrame(
                            geometry=[agent_line], crs=utils.wgs84
                        ).to_crs(metric_crs),
                        width_multiplier=1 / args.n_agents,
                    )
                    for agent_line in output_gdf.geometry
                ]
                swath_gdf = gpd.GeoDataFrame(
                    pd.concat(agent_swaths, ignore_index=True), crs=utils.wgs84
                )
            else:
                swath_gdf = m.simple_survey_line(
                    centerline_gdf.to_crs(metric_crs)
                )
            swath_gdf = swath_gdf.to_crs(metric_crs)
                        # print(output_gdf.area[1]/output_gdf.length[0])
            buf = 2500
            swath_gdf['geometry'] = swath_gdf['geometry'].buffer(buf, cap_style = 'round', join_style = 'round')
            swath_gdf['geometry'] = swath_gdf['geometry'].simplify(buf, preserve_topology=True)
            if args.n_agents != 2:
                swath_gdf['geometry'] = swath_gdf['geometry'].union_all()
            swath_gdf['geometry'] = swath_gdf['geometry'].buffer(-buf, cap_style = 'round', join_style = 'round')
            swath_gdf['geometry'] = swath_gdf['geometry'].simplify(.5*buf, preserve_topology=True)
            # print(output_gdf.crs, swath_gdf.crs)
            output_gdf = gpd.GeoDataFrame(pd.concat([output_gdf, swath_gdf.to_crs(utils.wgs84)], ignore_index=True), crs =  utils.wgs84)

            initial_area = m.unmapped_polygons.to_crs(metric_crs).area.sum()
            swath_union = gpd.GeoDataFrame(
                geometry=[swath_gdf.geometry.union_all()], crs=metric_crs
            )
            new_unmapped = m.unmapped_polygons.to_crs(metric_crs).overlay(swath_union, how='difference')
            metadata["newly_unmapped_area_m2"] = initial_area - new_unmapped.area.sum()

        geojson_dict = json.loads(output_gdf.to_json())
        fixed = antimeridian.fix_geojson(geojson_dict, fix_winding=True)
        output_gdf = gpd.GeoDataFrame.from_features(fixed["features"], crs=utils.wgs84)
       
        # fixed = geojson_dict
        geojson_dict = json.loads(output_gdf.to_json())
        geojson_dict["properties"] = metadata
        print(json.dumps(geojson_dict))