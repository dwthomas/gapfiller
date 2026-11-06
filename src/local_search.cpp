#include <iostream>
#include <limits>
#include <ogr_geometry.h>
#include <ogrsf_frmts.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "dataset_load.hpp"
#include "local_search_options.hpp"
#include "local_search_utils.hpp"

using local_search_utils::add_line_gap;
using local_search_utils::calculateCentroid;
using local_search_utils::closestPointOnPolygon;
using local_search_utils::find_alternative_paths;
using local_search_utils::IntersectionInfo;
using local_search_utils::score;
using local_search_utils::simplify;
using local_search_utils::transform_line;


int main(int argc, char *argv[]) {
  const auto parsed_options =
      local_search_cli::parse_command_line(argc, argv);
  if (!parsed_options) {
    return 1;
  }
  const auto &cli_options = *parsed_options;

  GDALAllRegister();

  auto [unmapped_polygons, unmapped_polygon_beam_widths, ds_unmapped] =
      load_unmapped(cli_options.unmapped_polygons_file);
  auto [unmapped_lines, unmapped_line_beam_widths, ds_unmapped_lines] =
      load_unmapped_lines(cli_options.unmapped_lines_file);
  auto [unmapped_points, unmapped_point_beam_widths, ds_unmapped_points] =
      load_unmapped_points(cli_options.unmapped_points_file);
  auto [land_polygons, ds_land] = load_land(cli_options.land_file);

  const double default_beam_width =
      unmapped_polygon_beam_widths.empty()
          ? 1.0
          : std::accumulate(unmapped_polygon_beam_widths.begin(),
                            unmapped_polygon_beam_widths.end(), 0.0) /
                unmapped_polygon_beam_widths.size();
  auto prepare_gap_beam_widths = [default_beam_width](
                                     std::vector<double> &beam_widths,
                                     std::size_t feature_count,
                                     const char *feature_type) {
    if (beam_widths.empty() && feature_count > 0) {
      std::cerr << "Unmapped " << feature_type
                << " have no beam-width scores; using " << default_beam_width
                << " for gap bonuses." << std::endl;
      beam_widths.assign(feature_count, default_beam_width);
    }
    if (beam_widths.size() != feature_count) {
      std::cerr << "Unmapped " << feature_type
                << " count does not match beam-width scores." << std::endl;
      return false;
    }
    return true;
  };
  if (!prepare_gap_beam_widths(unmapped_line_beam_widths, unmapped_lines.size(),
                               "lines") ||
      !prepare_gap_beam_widths(unmapped_point_beam_widths,
                               unmapped_points.size(), "points")) {
    GDALClose(ds_unmapped);
    GDALClose(ds_unmapped_lines);
    GDALClose(ds_unmapped_points);
    GDALClose(ds_land);
    return 1;
  }

  auto [initial_plan, line, planSrcSRS] = load_plan(cli_options.plan_file);


  // Define the source (WGS84) and target (Mollweide) spatial references
  OGRLayer *unmapped_layer = ds_unmapped->GetLayer(0);
  if (!unmapped_layer) {
    std::cerr << "Unmapped file does not contain any layers." << std::endl;
    GDALClose(ds_unmapped);
    GDALClose(ds_unmapped_lines);
    GDALClose(ds_unmapped_points);
    GDALClose(ds_land);
    return 1;
  }
  const OGRSpatialReference *spatialRef = unmapped_layer->GetSpatialRef();
  if (!spatialRef) {
    std::cerr << "Unmapped file does not have a valid spatial reference."
              << std::endl;
    GDALClose(ds_unmapped);
    GDALClose(ds_unmapped_lines);
    GDALClose(ds_unmapped_points);
    GDALClose(ds_land);
    return 1;
  }
  // Calculate centroid of the plan points
  auto [center_lat, center_lon] = calculateCentroid(initial_plan);
  
  std::cerr << "Plan centroid: (" << center_lat << ", " << center_lon << ")" << std::endl;
  
  // Create local azimuthal equidistant projection centered at plan centroid
  OGRSpatialReference dstSRS;
  // dstSRS = pyproj.CRS.from_proj4("+proj=merc +lon_0=-157.5 +datum=WGS84 +units=m")
  // std::string aeqd_proj4 = "+proj=aeqd +lat_0=" + std::to_string(center_lat) + 
  //                           " +lon_0=" + std::to_string(center_lon) + 
  //                           " +datum=WGS84 +units=m";
  
  // if (dstSRS.SetFromUserInput(aeqd_proj4.c_str()) != OGRERR_NONE) {
  //   std::cerr << "Failed to set local AEQD spatial reference: " << aeqd_proj4 << std::endl;
  //   GDALClose(ds_unmapped);
  //   GDALClose(ds_land);
  //   return 1;
  // }
  
  // std::cerr << "Using local AEQD projection: " << aeqd_proj4 << std::endl;
  std::string merc_proj4 = "+proj=merc +a=6378137 +b=6378137 +lon_0=" + std::to_string(center_lon) +
                          " +datum=WGS84 +units=m";
  if (dstSRS.SetFromUserInput(merc_proj4.c_str()) != OGRERR_NONE) {
    std::cerr << "Failed to set Mercator spatial reference: " << merc_proj4 << std::endl;
    GDALClose(ds_unmapped);
    GDALClose(ds_unmapped_lines);
    GDALClose(ds_unmapped_points);
    GDALClose(ds_land);
    return 1;
  }
  std::cerr << "Using Mercator projection: " << merc_proj4 << std::endl;

  // Create a coordinate transformation
  OGRCoordinateTransformation *transform =
      OGRCreateCoordinateTransformation(&planSrcSRS, &dstSRS);
  OGRCoordinateTransformation *inverse_transform =
      OGRCreateCoordinateTransformation(&dstSRS, &planSrcSRS);
  if (!transform) {
    std::cerr << "Failed to create coordinate transformation to AEQD." << std::endl;
    GDALClose(ds_unmapped);
    GDALClose(ds_unmapped_lines);
    GDALClose(ds_unmapped_points);
    GDALClose(ds_land);
    return 1;
  }

  // Transform unmapped polygons
  std::vector<OGRPolygon*> metric_unmapped_polygons;
  metric_unmapped_polygons.reserve(unmapped_polygons.size());

  for (auto* polygon : unmapped_polygons) {
    metric_unmapped_polygons.push_back(static_cast<OGRPolygon*>(polygon->clone()));
  }
  for (auto &polygon : metric_unmapped_polygons) {
    if (polygon->transform(transform) != OGRERR_NONE) {
      std::cerr << "Failed to transform an unmapped polygon to AEQD." << std::endl;
    }
  }

  std::vector<OGRLineString *> metric_unmapped_lines;
  metric_unmapped_lines.reserve(unmapped_lines.size());
  for (auto *feature_line : unmapped_lines) {
    auto *metric_line =
        static_cast<OGRLineString *>(feature_line->clone());
    if (metric_line->transform(transform) != OGRERR_NONE) {
      std::cerr << "Failed to transform an unmapped line to AEQD."
                << std::endl;
    }
    metric_unmapped_lines.push_back(metric_line);
  }

  std::vector<OGRPoint *> metric_unmapped_points;
  metric_unmapped_points.reserve(unmapped_points.size());
  for (auto *feature_point : unmapped_points) {
    auto *metric_point = static_cast<OGRPoint *>(feature_point->clone());
    if (metric_point->transform(transform) != OGRERR_NONE) {
      std::cerr << "Failed to transform an unmapped point to AEQD."
                << std::endl;
    }
    metric_unmapped_points.push_back(metric_point);
  }

  // order unmapped polygons by area descending
  std::sort(metric_unmapped_polygons.begin(), metric_unmapped_polygons.end(),
            [](OGRPolygon *a, OGRPolygon *b) {
              return a->get_Area() > b->get_Area();
            });

  // Transform land polygons
  for (auto &polygon : land_polygons) {
    if (polygon->transform(transform) != OGRERR_NONE) {
      std::cerr << "Failed to transform a land polygon to AEQD." << std::endl;
    }
  }

  OGRCoordinateTransformation *planTransform =
      OGRCreateCoordinateTransformation(&planSrcSRS, &dstSRS);
  if (!planTransform) {
    std::cerr << "Failed to create plan coordinate transformation to AEQD." << std::endl;
    GDALClose(ds_unmapped);
    GDALClose(ds_unmapped_lines);
    GDALClose(ds_unmapped_points);
    GDALClose(ds_land);
    return 1;
  }

  // // Transform plan points from plan CRS (WGS84 by default) to AEQD.
  // for (auto &point : initial_plan) {
  //   double x = point.getX();
  //   double y = point.getY();
  //   if (!planTransform->Transform(1, &x, &y)) {
  //     std::cerr << "Failed to transform a plan point to AEQD." << std::endl;
  //     continue;
  //   }
  //   point.setX(x);
  //   point.setY(y);
  // }
  auto metric_plan = transform_line(initial_plan, planTransform);

  auto initial_score = score(metric_plan, metric_unmapped_polygons, unmapped_polygon_beam_widths, inverse_transform);
  std::cerr << "Initial plan score: " << initial_score << std::endl;
  int iternum = 0;

  std::vector<OGRPolygon> used_polys;
  std::vector<bool> used_line_gaps(metric_unmapped_lines.size(), false);
  std::vector<bool> used_point_gaps(metric_unmapped_points.size(), false);
  double accumulated_gap_bonus = 0.0;

  double initial_length = -1;

  while (true) {
    ++iternum;
    // Find intersections between the path and unmapped polygons
    std::vector<IntersectionInfo> intersection_infos;

    OGRLineString line;
    for (const auto& p : simplify(metric_plan)) {
      line.addPoint(&p);    
    }

    if (initial_length < 0) {
      initial_length = line.get_Length();
    }

    double remainder = cli_options.budget - line.get_Length();
    std::cerr << "C++ line length: " << line.get_Length() << " budget: "
              << cli_options.budget << " remainder: " << remainder
              << std::endl;
    int non_intersecting = 5;

    for (auto polygon : metric_unmapped_polygons) {
      if (std::find(used_polys.begin(), used_polys.end(), *polygon) != used_polys.end()) {
        continue;
      }
      for (size_t i = 0; i < metric_plan.size() - 1; i++) {
        OGRLineString segment;
        segment.addPoint(&metric_plan[i]);
        segment.addPoint(&metric_plan[i + 1]);

        if (polygon->Intersects(&segment)) {

          std::cerr << "Found intersection with polygon at segment "
                    << i << std::endl;

          OGRGeometry *intersection = polygon->Intersection(&segment);
          if (intersection != nullptr &&
              wkbFlatten(intersection->getGeometryType()) == wkbLineString) {
            IntersectionInfo info;
            info.intersection = (OGRLineString *)intersection->clone();
            info.path_segment = (OGRLineString *)segment.clone();
            info.unmapped_region = (OGRPolygon *)polygon->clone();
            info.index = i;
            intersection_infos.push_back(info);
          }
          OGRGeometryFactory::destroyGeometry(intersection);
        }
        if (non_intersecting > 0 && segment.Distance(polygon) <= remainder) {
            non_intersecting--;
          std::cerr << "Distance to polygon: " << segment.Distance(polygon) << "Remaining budget: " << remainder << std::endl;
            OGRPoint start, end;
            segment.StartPoint(&start);
            segment.EndPoint(&end);

            OGRPoint * closestToStart = closestPointOnPolygon(*polygon, start);
            OGRPoint * closestToEnd = closestPointOnPolygon(*polygon, end);
            
            OGRLineString * intersection = new OGRLineString();
            intersection->addPoint(closestToStart);
            intersection->addPoint(closestToEnd);

            IntersectionInfo info;
            info.intersection = intersection;
            info.path_segment = (OGRLineString *)segment.clone();
            info.unmapped_region = (OGRPolygon *)polygon->clone();
            info.index = i;
            intersection_infos.push_back(info);
        }
      }
    }

    // Store the intersection info in a vector for later use
    for (const auto &info : intersection_infos) {
      // std::cout << "Intersection LineString: ";
      // for (int i = 0; i < info.intersection->getNumPoints(); ++i) {
      //   OGRPoint point;
      //   info.intersection->getPoint(i, &point);
      //   std::cout << "(" << point.getX() << ", " << point.getY() << ") ";
      // }
      //  std::cout << std::endl;
      // call find alternative paths
    }

    std::vector<std::vector<OGRPoint>> options;
    std::vector<OGRPolygon *> used_polys_vector;
    std::vector<int> option_line_indices;
    std::vector<int> option_point_indices;

    struct NearbyGapLine {
      std::size_t line_index;
      std::size_t insertion_index;
      double centroid_distance;
    };
    std::vector<NearbyGapLine> nearby_gap_lines;
    if (cli_options.consider_gaps > 0) {
      for (std::size_t line_index = 0;
           line_index < metric_unmapped_lines.size(); ++line_index) {
        if (used_line_gaps[line_index] ||
            metric_unmapped_lines[line_index]->getNumPoints() < 2) {
          continue;
        }

        OGRPoint centroid;
        if (metric_unmapped_lines[line_index]->Centroid(&centroid) !=
            OGRERR_NONE) {
          std::cerr << "Failed to calculate centroid for unmapped line "
                    << line_index << "." << std::endl;
          continue;
        }

        double closest_distance = std::numeric_limits<double>::infinity();
        std::size_t closest_vertex = 0;
        for (std::size_t vertex = 0; vertex + 1 < metric_plan.size(); ++vertex) {
          const double distance = centroid.Distance(&metric_plan[vertex]);
          if (distance < closest_distance) {
            closest_distance = distance;
            closest_vertex = vertex;
          }
        }
        if (closest_distance <= remainder) {
          nearby_gap_lines.push_back(
              {line_index, closest_vertex, closest_distance});
        }
      }

      std::sort(nearby_gap_lines.begin(), nearby_gap_lines.end(),
                [](const NearbyGapLine &a, const NearbyGapLine &b) {
                  return a.centroid_distance < b.centroid_distance;
                });
      if (nearby_gap_lines.size() >
          static_cast<std::size_t>(cli_options.consider_gaps)) {
        nearby_gap_lines.resize(
            static_cast<std::size_t>(cli_options.consider_gaps));
      }

      for (const auto &nearby_line : nearby_gap_lines) {
        auto candidate =
            add_line_gap(metric_plan, nearby_line.insertion_index,
                         *metric_unmapped_lines[nearby_line.line_index]);
        OGRLineString candidate_line;
        for (const auto &point : simplify(candidate)) {
          candidate_line.addPoint(&point);
        }
        if (candidate_line.get_Length() < cli_options.budget) {
          options.push_back(std::move(candidate));
          used_polys_vector.push_back(nullptr);
          option_line_indices.push_back(
              static_cast<int>(nearby_line.line_index));
          option_point_indices.push_back(-1);
        }
      }
    }

    struct NearbyGapPoint {
      std::size_t point_index;
      std::size_t insertion_index;
      double distance;
    };
    std::vector<NearbyGapPoint> nearby_gap_points;
    if (cli_options.consider_gaps > 0) {
      for (std::size_t point_index = 0;
           point_index < metric_unmapped_points.size(); ++point_index) {
        if (used_point_gaps[point_index]) {
          continue;
        }

        double closest_distance = std::numeric_limits<double>::infinity();
        std::size_t closest_vertex = 0;
        for (std::size_t vertex = 0; vertex < metric_plan.size(); ++vertex) {
          const double distance =
              metric_unmapped_points[point_index]->Distance(&metric_plan[vertex]);
          if (distance < closest_distance) {
            closest_distance = distance;
            closest_vertex = vertex;
          }
        }
        if (closest_distance <= remainder) {
          nearby_gap_points.push_back(
              {point_index, closest_vertex, closest_distance});
        }
      }

      std::sort(nearby_gap_points.begin(), nearby_gap_points.end(),
                [](const NearbyGapPoint &a, const NearbyGapPoint &b) {
                  return a.distance < b.distance;
                });
      if (nearby_gap_points.size() >
          static_cast<std::size_t>(cli_options.consider_gaps)) {
        nearby_gap_points.resize(
            static_cast<std::size_t>(cli_options.consider_gaps));
      }

      constexpr double point_gap_half_length_m = 1000.0;
      for (const auto &nearby_point : nearby_gap_points) {
        const auto &point = *metric_unmapped_points[nearby_point.point_index];
        const auto insertion_index = nearby_point.insertion_index;
        double dx = metric_plan[insertion_index].getX() - point.getX();
        double dy = metric_plan[insertion_index].getY() - point.getY();
        double direction_length = std::hypot(dx, dy);
        if (direction_length == 0.0 && insertion_index > 0) {
          dx = metric_plan[insertion_index - 1].getX() - point.getX();
          dy = metric_plan[insertion_index - 1].getY() - point.getY();
          direction_length = std::hypot(dx, dy);
        }
        if (direction_length == 0.0 &&
            insertion_index + 1 < metric_plan.size()) {
          dx = metric_plan[insertion_index + 1].getX() - point.getX();
          dy = metric_plan[insertion_index + 1].getY() - point.getY();
          direction_length = std::hypot(dx, dy);
        }
        if (direction_length == 0.0) {
          continue;
        }

        OGRPoint toward_plan(
            point.getX() + point_gap_half_length_m * dx / direction_length,
            point.getY() + point_gap_half_length_m * dy / direction_length);
        OGRPoint away_from_plan(
            point.getX() - point_gap_half_length_m * dx / direction_length,
            point.getY() - point_gap_half_length_m * dy / direction_length);
        OGRLineString pseudo_line;
        pseudo_line.addPoint(&toward_plan);
        pseudo_line.addPoint(&away_from_plan);

        auto candidate =
            add_line_gap(metric_plan, insertion_index, pseudo_line);
        OGRLineString candidate_line;
        for (const auto &candidate_point : simplify(candidate)) {
          candidate_line.addPoint(&candidate_point);
        }
        if (candidate_line.get_Length() < cli_options.budget) {
          options.push_back(std::move(candidate));
          used_polys_vector.push_back(nullptr);
          option_line_indices.push_back(-1);
          option_point_indices.push_back(
              static_cast<int>(nearby_point.point_index));
        }
      }
    }


    for (std::size_t i = 0; i < intersection_infos.size(); i++) {
        OGRPoint source;
        
        if (intersection_infos[i].intersection->getNumPoints() < 1) {
            // std::cerr << "Intersection has no points." << std::endl;
            continue;
        }
        intersection_infos[i].intersection->getPoint(0, &source);
        OGRPoint target;
        intersection_infos[i].intersection->getPoint(
            intersection_infos[i].intersection->getNumPoints() - 1, &target);

        OGRPoint *next_point = &metric_plan[std::min(
            (std::size_t)intersection_infos[i].index + 1, metric_plan.size() - 1)];
        // std::cerr << "target point: (" << target.getX() << ", " << target.getY()
        //         << ")" << std::endl;

        // std::cerr << "polygon: ";
        OGRLinearRing *ring =
            intersection_infos[i].unmapped_region->getExteriorRing();

        for (int i = 0; i < ring->getNumPoints(); i++) {
        OGRPoint p;
        ring->getPoint(i, &p);
        // std::cerr << i << ": (" << p.getX() << ", " << p.getY() << ")\n";
        }
        // std::cerr << std::endl;
        // for (std::size_t i = 0; i < plan.size() - 1; ++i) {
        //     std::cerr << "checking plan point: (" << plan[i].getX() << ",
        //     " << plan[i].getY() << ")" << std::endl; if
        //     (plan[i].Equals(&target)) {
        //         next_point = &plan[i + 1];
        //         break;
        //     }
        // }
        // if(next_point == nullptr){
        //     std::cerr << "Next point not found after target point." << std::endl;
        //     break;
        // }

        std::pair<OGRLineString *, OGRLineString *> alt_paths =
            find_alternative_paths(*intersection_infos[i].unmapped_region, source,
                                target, *next_point);

        std::vector<OGRPoint> new_plan1;

        for (std::size_t j = 0; j <= intersection_infos[i].index; j++) {
            new_plan1.push_back(metric_plan[j]);
        }
        for (int k = 0; k < alt_paths.first->getNumPoints(); k++) {
            OGRPoint point;
            alt_paths.first->getPoint(k, &point);
            new_plan1.push_back(point);
        }
        for (std::size_t j = intersection_infos[i].index + 2; j < metric_plan.size(); j++) {
            new_plan1.push_back(metric_plan[j]);
        }

        std::vector<OGRPoint> new_plan2;

        for (std::size_t j = 0; j <= intersection_infos[i].index; j++) {
            new_plan2.push_back(metric_plan[j]);
        }
        for (int k = 0; k < alt_paths.second->getNumPoints(); k++) {
            OGRPoint point;
            alt_paths.second->getPoint(k, &point);
            new_plan2.push_back(point);
        }
        for (std::size_t j = intersection_infos[i].index + 2; j < metric_plan.size(); j++) {
            new_plan2.push_back(metric_plan[j]);
        }

        // std::cout << "New Plan 1:" << std::endl;
        // for (const auto &point : new_plan1) {
        //     std::cout << "Point: (" << point.getX() << ", " << point.getY() << ")"
        //             << std::endl;
        // }
        // std::cout << "New Plan 2:" << std::endl;
        // for (const auto &point : new_plan2) {           
        //     std::cout << "Point: (" << point.getX() << ", " << point.getY() << ")"
        //             << std::endl;
        // }

        for (const auto& path : {new_plan1, new_plan2}) {
            OGRLineString line;
            for (const auto& p : simplify(path)) {
                line.addPoint(&p);    
            }
            double length = line.get_Length();
            std::cerr << "Path length: " << length
                      << " budget: " << cli_options.budget << std::endl;
            if(length < cli_options.budget){
                options.push_back(path);
                used_polys_vector.push_back(
                    intersection_infos[i].unmapped_region);
                option_line_indices.push_back(-1);
                option_point_indices.push_back(-1);
            }
        }
    }
    if (options.empty()) {
      std::cerr << "No more intersections found. Exiting local improvement."
                << std::endl;
      OGRLineString final_line;
      for (const auto &p : simplify(metric_plan)) {
          final_line.addPoint(&p);    
      }
      double dur = final_line.get_Length();
      // convert to wgs84
      OGRCoordinateTransformation *inv_transform =
          OGRCreateCoordinateTransformation(&dstSRS, &planSrcSRS);
      final_line.transform(inv_transform);
      char *wkt = nullptr;
      final_line.exportToWkt(&wkt);
      std::cout << wkt << std::endl;
      std::cerr << "Final plan duration: " << dur << std::endl;
      std::cerr << "Initial plan duration: " << initial_length << std::endl; 
      std::cerr << "Improvement: " << dur - initial_length << std::endl;
      std::cerr << "Budget: " << cli_options.budget << std::endl << std::endl;
      CPLFree(wkt);
      break;
    } 
    // std::cerr << "options:\n";
    double best_score = 0.0;
    long best_i = -1;
    long i = 0;
    for (std::size_t option_index = 0; option_index < options.size();
         ++option_index) {
      const auto &option = options[option_index];
      double s = score(option, metric_unmapped_polygons, unmapped_polygon_beam_widths, inverse_transform);
      s += accumulated_gap_bonus;
      if (option_line_indices[option_index] >= 0) {
        const auto line_index =
            static_cast<std::size_t>(option_line_indices[option_index]);
        s += cli_options.gap_bonus_km * 1000.0 *
             unmapped_line_beam_widths[line_index];
      }
      if (option_point_indices[option_index] >= 0) {
        const auto point_index =
            static_cast<std::size_t>(option_point_indices[option_index]);
        s += cli_options.gap_bonus_km * 1000.0 *
             unmapped_point_beam_widths[point_index];
      }
        // std::cerr << "Score: " << s;
        // for (const auto &point : option) {
        //     std::cerr << " (" << point.getX() << ", " << point.getY() << ")\t";
        // }
        // std::cerr << std::endl;
        if( s > best_score){
            best_score = s;
            best_i = i;
        }
        ++i;
    }
    if (best_score > score(metric_plan, metric_unmapped_polygons, unmapped_polygon_beam_widths, inverse_transform) + accumulated_gap_bonus) {
      metric_plan = options[best_i];
      if (option_line_indices[best_i] >= 0) {
        const auto line_index =
            static_cast<std::size_t>(option_line_indices[best_i]);
        used_line_gaps[line_index] = true;
        accumulated_gap_bonus +=
            cli_options.gap_bonus_km * 1000.0 *
            unmapped_line_beam_widths[line_index];
      }
      if (option_point_indices[best_i] >= 0) {
        const auto point_index =
            static_cast<std::size_t>(option_point_indices[best_i]);
        used_point_gaps[point_index] = true;
        accumulated_gap_bonus +=
            cli_options.gap_bonus_km * 1000.0 *
            unmapped_point_beam_widths[point_index];
      }
      if (used_polys_vector[best_i]) {
        char *wkt = nullptr;
        used_polys_vector[best_i]->exportToWkt(&wkt);
        std::cerr << "Using polygon WKT: " << (wkt ? wkt : "") << std::endl;
        CPLFree(wkt);
      } else {
        std::cerr << "Using polygon: nullptr" << std::endl;
      }
      if (used_polys_vector[best_i]) {
        used_polys.push_back(*used_polys_vector[best_i]);
      }
      std::cerr << "Updated plan with improved score: " << best_score
                 << std::endl;
       // print geojson of plan
      OGRLineString final_line;
      for (const auto &p : metric_plan) {
          final_line.addPoint(&p);    
      }
      // convert to wgs84
      std::cerr << "Iteration number: " << iternum << std::endl;
      if(iternum >= 50){ // need to not allow doubling back
        OGRCoordinateTransformation *inv_transform =
            OGRCreateCoordinateTransformation(&dstSRS, &planSrcSRS);
        final_line.transform(inv_transform);
        char *wkt = nullptr;
        final_line.exportToWkt(&wkt);
        std::cout << wkt << std::endl;
        CPLFree(wkt);
        break;     
      }
    } 
    else {
      std::cerr << "No improvement found. Exiting local improvement."
                << std::endl;
      // print geojson of plan
      OGRLineString final_line;
      for (const auto &p : metric_plan) {
          final_line.addPoint(&p);    
      }
      // convert to wgs84
      OGRCoordinateTransformation *inv_transform =
          OGRCreateCoordinateTransformation(&dstSRS, &planSrcSRS);
      final_line.transform(inv_transform);
      char *wkt = nullptr;
      final_line.exportToWkt(&wkt);
      std::cout << wkt << std::endl;
      CPLFree(wkt);
      break;                
    }

    // break;
  }

  OCTDestroyCoordinateTransformation(transform);

 
  for (auto *feature_line : unmapped_lines) {
    OGRGeometryFactory::destroyGeometry(feature_line);
  }
  for (auto *feature_line : metric_unmapped_lines) {
    OGRGeometryFactory::destroyGeometry(feature_line);
  }
  for (auto *feature_point : unmapped_points) {
    OGRGeometryFactory::destroyGeometry(feature_point);
  }
  for (auto *feature_point : metric_unmapped_points) {
    OGRGeometryFactory::destroyGeometry(feature_point);
  }
  GDALClose(ds_unmapped);
  GDALClose(ds_unmapped_lines);
  GDALClose(ds_unmapped_points);
  GDALClose(ds_land);

  return 0;
}