#pragma once

#include <ogr_geometry.h>
#include <ogrsf_frmts.h>


std::tuple<std::vector<OGRPolygon *>, std::vector<double>, GDALDataset *> load_unmapped(const std::string& unmapped_file);

std::tuple<std::vector<OGRPolygon *>, GDALDataset *> load_land(const std::string& land_file);

std::tuple<std::vector<OGRPoint>, OGRLineString *, OGRSpatialReference> load_plan(const std::string& plan_file);

