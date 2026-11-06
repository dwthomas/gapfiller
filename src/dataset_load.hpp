#pragma once

#include "s2/mutable_s2shape_index.h" // MutableS2ShapeIndex
#include "s2/s2latlng.h"  // S2LatLng (for lat/lon -> S2Point conversion)
#include "s2/s2point.h"   // S2Point
#include "s2/s2polygon.h" // S2Polygon (includes S2Polygon::Shape)
#include "s2/s2shape.h" // S2Shape (usually pulled in transitively, but explicit is safer)
#include <gdal_priv.h>
#include <ogrsf_frmts.h>
std::tuple<std::vector<OGRPolygon *>, std::vector<double>, GDALDataset *>
load_unmapped(const std::string &unmapped_file);
std::tuple<std::vector<OGRPoint *>, std::vector<double>, GDALDataset *>
load_unmapped_points(const std::string &unmapped_points_file);
std::tuple<std::vector<OGRLineString *>, std::vector<double>, GDALDataset *>
load_unmapped_lines(const std::string &unmapped_lines_file);
std::tuple<std::vector<OGRPolygon *>, GDALDataset *>
load_land(const std::string &land_file);
std::tuple<std::vector<OGRPoint>, OGRLineString *, OGRSpatialReference>
load_plan(const std::string &plan_file);
