#pragma once

#include <fstream>
#include <iostream>

#include "json.hpp"

#include "dataset_load.hpp" 

std::tuple<std::vector<OGRPolygon *>, std::vector<double>, GDALDataset *> load_unmapped(const std::string& unmapped_file){
  GDALDataset *ds_unmapped = (GDALDataset *)GDALOpenEx(
      unmapped_file.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr);

  if (!ds_unmapped) {
    std::cerr << "Failed to open unmapped file: " << unmapped_file << std::endl;
    exit(1);
  }

  std::vector<OGRPolygon *> unmapped_polygons;
  std::vector<double> unmapped_polygon_beam_widths;


    // Extract polygons from unmapped dataset
  for (int i = 0; i < ds_unmapped->GetLayerCount(); ++i) {
    OGRLayer *layer = ds_unmapped->GetLayer(i);
    layer->ResetReading();
    OGRFeature *feature = nullptr;

    while ((feature = layer->GetNextFeature()) != nullptr) {
      OGRGeometry *geom = feature->GetGeometryRef();
      if (geom != nullptr &&
          wkbFlatten(geom->getGeometryType()) == wkbPolygon) {
        unmapped_polygons.push_back((OGRPolygon *)geom->clone());
      }
      OGRFeature::DestroyFeature(feature);
    }
  }
  
  // Read unmapped_polygon_beam_widths from JSON properties
  std::ifstream unmapped_json_file(unmapped_file);
  if (unmapped_json_file.is_open()) {
    try {
      nlohmann::json j = nlohmann::json::parse(unmapped_json_file);
      
      // std::cerr << "[DEBUG] Top-level JSON keys:\n";
      // for (auto& [key, value] : j.items()) {
      //   std::cerr << "  - " << key << "\n";
      // }

      if (j.contains("features") && j["features"].is_array() && j["features"].size() > 0) {
        // std::cerr << "[DEBUG] First feature keys:\n";
        for (auto& [key, value] : j["features"][0].items()) {
          // std::cerr << "  - " << key << "\n";
        }
        if (j["features"][0].contains("properties") && j["features"][0]["properties"].is_object()) {
          // std::cerr << "[DEBUG] First feature properties keys:\n";
          for (auto& [key, value] : j["features"][0]["properties"].items()) {
            // std::cerr << "  - " << key << "\n";
          }
        }
      }

      if (j.contains("properties") && j["properties"].is_object()) {
        auto& props = j["properties"];
        if (props.contains("unmapped_scores")) {
          
          if (props["unmapped_scores"].is_object()) {
            for (auto& [key, value] : props["unmapped_scores"].items()) {
              unmapped_polygon_beam_widths.push_back(value.get<double>());
            }
          } else if (props["unmapped_scores"].is_array()) {
            for (auto& value : props["unmapped_scores"]) {
              unmapped_polygon_beam_widths.push_back(value.get<double>());
            }
          } else {
            std::cerr << "[DEBUG] unmapped_scores is neither an object nor an array!\n";
          }
        }
      }
    } catch (const std::exception& e) {
      std::cerr << "Error reading unmapped_scores from JSON: " << e.what() << std::endl;
    }
    unmapped_json_file.close();
  }
  std::cerr << "Unmapped polygon beam widths: ";
  for (auto bw: unmapped_polygon_beam_widths) {
     std::cerr << bw << " ";
  }
  std::cerr << std::endl;
  return std::make_tuple(unmapped_polygons, unmapped_polygon_beam_widths, ds_unmapped);
}

std::tuple<std::vector<OGRPolygon *>, GDALDataset *> load_land(const std::string& land_file){
  GDALDataset *ds_land = (GDALDataset *)GDALOpenEx(
      land_file.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr);
  if (!ds_land) {
    std::cerr << "Failed to open land file: " << land_file << std::endl;
    GDALClose(ds_land);
    exit(1);
  }
  
  std::vector<OGRPolygon *> land_polygons;

  // Extract polygons from land dataset
  for (int i = 0; i < ds_land->GetLayerCount(); ++i) {
    OGRLayer *layer = ds_land->GetLayer(i);
    layer->ResetReading();
    OGRFeature *feature = nullptr;
    while ((feature = layer->GetNextFeature()) != nullptr) {
      OGRGeometry *geom = feature->GetGeometryRef();
      if (geom != nullptr &&
          wkbFlatten(geom->getGeometryType()) == wkbPolygon) {
        land_polygons.push_back((OGRPolygon *)geom->clone());
      }
      OGRFeature::DestroyFeature(feature);
    }
  }
  return std::make_tuple(land_polygons, ds_land);
}