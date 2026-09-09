#pragma once

#include "s2/s2polygon.h"           // S2Polygon (includes S2Polygon::Shape)
#include "s2/s2point.h"             // S2Point
#include "s2/s2latlng.h"            // S2LatLng (for lat/lon -> S2Point conversion)
#include "s2/mutable_s2shape_index.h"  // MutableS2ShapeIndex
#include "s2/s2shape.h"             // S2Shape (usually pulled in transitively, but explicit is safer)

std::tuple<S2ShapeIndex, std::vector<double>> load_unmapped(const std::string& unmapped_file);

S2ShapeIndex load_land(const std::string& land_file);

S2Polyline load_plan(const std::string& plan_file);

