#include <iostream>
#include <limits>
#include <ogr_geometry.h>
#include <ogrsf_frmts.h>

#include <algorithm>
#include <boost/program_options.hpp>
#include <cmath>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
#include <unordered_map>
#include <functional>

namespace po = boost::program_options;

// Calculate centroid of points, handling date line wrapping
std::pair<double, double> calculateCentroid(const std::vector<OGRPoint>& points) {
    if (points.empty()) {
        return {0.0, 0.0};
    }
    
    double sum_lat = 0.0;
    double sum_lon_x = 0.0;  // sum of cos(lon)
    double sum_lon_y = 0.0;  // sum of sin(lon)
    
    for (const auto& point : points) {
        double lat = point.getY();
        double lon = point.getX();
        
        sum_lat += lat;
        
        // Convert to radians and use unit vector approach for longitude
        // This handles date line wrapping correctly
        double lon_rad = lon * M_PI / 180.0;
        sum_lon_x += std::cos(lon_rad);
        sum_lon_y += std::sin(lon_rad);
    }
    
    double n = static_cast<double>(points.size());
    double avg_lat = sum_lat / n;
    
    // Calculate average longitude from the unit vectors
    double avg_lon_rad = std::atan2(sum_lon_y / n, sum_lon_x / n);
    double avg_lon = avg_lon_rad * 180.0 / M_PI;
    
    return {avg_lat, avg_lon};
}

OGRPoint * closestPointOnPolygon(const OGRPolygon &polygon, const OGRPoint &point) {
    OGRPoint * closestPoint = new OGRPoint();
    double minDist = std::numeric_limits<double>::infinity();

    const OGRLinearRing * ring = polygon.getExteriorRing();
    if (!ring) {
        std::cerr << "Polygon does not have an exterior ring." << std::endl;
        return nullptr;
    }

    int numPoints = ring->getNumPoints();
    for (int i = 0; i < numPoints - 1; ++i) {
        OGRPoint p1, p2;
        ring->getPoint(i, &p1);
        ring->getPoint(i + 1, &p2);

        OGRLineString edge;
        edge.addPoint(&p1);
        edge.addPoint(&p2);

        OGRPoint projectedPoint;
        // edge.Project(&point);

        double dist_along = edge.Project(&point); // distance from start of edge to projection
        double edgeLen = edge.get_Length();
        if (edgeLen > 0.0) {
            double t = dist_along / edgeLen;
            t = std::max(0.0, std::min(1.0, t)); // clamp just in case
            projectedPoint.setX(p1.getX() + t * (p2.getX() - p1.getX()));
            projectedPoint.setY(p1.getY() + t * (p2.getY() - p1.getY()));
        } else {
            // degenerate edge: use one endpoint
            projectedPoint.setX(p1.getX());
            projectedPoint.setY(p1.getY());
        }

        double dist = projectedPoint.Distance(&point);


        if (dist < minDist) {
            minDist = dist;
            closestPoint->setX(projectedPoint.getX());
            closestPoint->setY(projectedPoint.getY());
        }
    }

    return closestPoint;
}

std::pair<OGRPolygon, int> add_to_polygon(const OGRPolygon &polygon,
                                          const OGRPoint &point) {
  std::vector<OGRPoint> pgon_points;
  int min_i = 0;
  double min_dist = std::numeric_limits<double>::infinity();
  int n = polygon.getExteriorRing()->getNumPoints();
  // std::cerr << "Polygon has " << n << " points." << std::endl;
  for (int i = 0; i < n; i++) {
    OGRPoint p1, p2;
    polygon.getExteriorRing()->getPoint(i, &p1);
    pgon_points.push_back(p1);
    polygon.getExteriorRing()->getPoint((i + 1) % n, &p2);
    OGRLineString edge;
    edge.addPoint(&p1);
    edge.addPoint(&p2);
    double dist = edge.Distance(&point);
    if (dist < min_dist) {
      min_dist = dist;
      min_i = i + 1;
    }
  }
  pgon_points.insert(pgon_points.begin() + min_i, point);
  OGRLinearRing ring;
  for (const auto &p : pgon_points) {
    ring.addPoint(&p);
  }
  ring.closeRings();
  OGRPolygon new_polygon;
  new_polygon.addRing(&ring);
  return {new_polygon, min_i};
}

double score(const std::vector<OGRPoint>& path, const std::vector<OGRPolygon *> &unmapped_polygons) {
    OGRLineString line;
    for (const auto& p : path) {
        line.addPoint(&p);    
    }
    double length = line.get_Length();
    return length;
}


std::vector<OGRPoint> simplify(const std::vector<OGRPoint>& path){
  // remove loops in path
  std::unordered_map<double, std::unordered_map<double, std::size_t>> last_seen;

  for(std::size_t i = 0; i < path.size(); i++){
    double x = path[i].getX();
    double y = path[i].getY();
    last_seen[x][y] = i;
  }

  std::vector<OGRPoint> simplified_path;

  std::size_t i = 0;
  while (i < path.size()){
    double x = path[i].getX();
    double y = path[i].getY();
    simplified_path.push_back(path[i]);
    i = last_seen[x][y] + 1;
  }
  return simplified_path;
}


std::pair<OGRLineString *, OGRLineString *>
find_alternative_paths(const OGRPolygon &polygon, const OGRPoint &source,
                       const OGRPoint &target, const OGRPoint &next_point) {
  int source_index;
  // std::cerr << "pgonsize " << polygon.getExteriorRing()->getNumPoints()
  //           << std::endl;
  auto x = add_to_polygon(polygon, target);
  // std::cerr << "pgonsize " << x.first.getExteriorRing()->getNumPoints()
  //           << std::endl;
  x = add_to_polygon(x.first, source);
  // std::cerr << "pgonsize " << x.first.getExteriorRing()->getNumPoints()
  //           << std::endl;
  source_index = x.second;
  int min_index = 0;
  double min_dist = std::numeric_limits<double>::infinity();
  for (int i = 0; i < x.first.getExteriorRing()->getNumPoints(); i++) {
    OGRPoint p;
    x.first.getExteriorRing()->getPoint(i, &p);
    double dist = p.Distance(&next_point);
    if (dist < min_dist) {
      min_dist = dist;
      min_index = i;
    }
  }

  std::vector<OGRPoint> line1;
  int i = source_index;
  while (i != min_index) {
    OGRPoint p;
    x.first.getExteriorRing()->getPoint(i, &p);
    line1.push_back(p);
    i = (i + 1) % x.first.getExteriorRing()->getNumPoints();
  }
  line1.push_back(next_point);

  std::vector<OGRPoint> line2;
  i = source_index;
  // std::cerr << "size = " << x.first.getExteriorRing()->getNumPoints()
  //           << std::endl;
  while (i != min_index) {
    OGRPoint p;
    x.first.getExteriorRing()->getPoint(i, &p);
    line2.push_back(p);
    // std::cerr << "l2 point: (" << i << " " << p.getX() << ", " << p.getY()
    //           << ")" << std::endl;
    i = (i - 1);
    if (i < 0) {
      i = x.first.getExteriorRing()->getNumPoints() - 1;
    }
  }
  line2.push_back(next_point);
  OGRLineString *ls1 = new OGRLineString();
  for (const auto &p : line1){
    ls1->addPoint(&p);
  }

  OGRLineString *ls2 = new OGRLineString();
  for (const auto &p : line2){
    ls2->addPoint(&p);
  }

  return std::make_pair(ls1, ls2);
}

struct IntersectionInfo {
  // the intersection
  OGRLineString *intersection;
  // the segment of the path that produced the intersection
  OGRLineString *path_segment;
  // the unmapped region that produced the intersection
  OGRPolygon * unmapped_region;
  std::size_t index;
};

std::pair<double, OGRLineString *>
local_improvement(const IntersectionInfo &info) {
  // Find the endpoints of the intersection on the path segment
  OGRPoint start_point, end_point;
  info.intersection->StartPoint(&start_point);
  info.intersection->EndPoint(&end_point);

  // Find the edge where the start_point lies and traverse the edges to form two
  // lines
  OGRLineString *line1 = new OGRLineString();
  OGRLineString *line2 = new OGRLineString();

  OGRLinearRing *ring = info.unmapped_region->getExteriorRing();
  if (!ring) {
    std::cerr << "Polygon does not have an exterior ring." << std::endl;
    return {0.0, nullptr};
  }

  std::cerr << "Polygon: ";
  for (int i = 0; i < ring->getNumPoints(); ++i) {
    OGRPoint p;
    ring->getPoint(i, &p);
    std::cerr << "(" << p.getX() << ", " << p.getY() << ") ";
  }
  std::cerr << std::endl;

  int num_points = ring->getNumPoints();
  bool start_found = false;

  // Traverse the ring to find the edge containing start_point
  for (int i = 0; i < num_points - 1; ++i) {
    OGRPoint p1, p2;
    ring->getPoint(i, &p1);
    ring->getPoint(i + 1, &p2);

    OGRLineString edge;
    edge.addPoint(&p1);
    edge.addPoint(&p2);

    if (start_point.Within(&edge)) {
      start_found = true;

      // Build line1 from start_point to end_point
      line1->addPoint(&start_point);
      for (int j = i + 1; j < num_points; j++) {
        OGRPoint next_point;
        ring->getPoint(j, &next_point);
        line1->addPoint(&next_point);
        if (next_point.Equals(&end_point)) {
          break;
        }
      }

      // Build line2 in reverse direction from start_point to end_point
      line2->addPoint(&start_point);
      for (int j = i; j >= 0; --j) {
        OGRPoint prev_point;
        ring->getPoint(j, &prev_point);
        line2->addPoint(&prev_point);
        if (prev_point.Equals(&end_point)) {
          break;
        }
      }

      break;
    }
  }

  if (!start_found) {
    std::cerr << "Start point is not on the polygon's edge." << std::endl;
    delete line1;
    delete line2;
    return {0.0, nullptr};
  }

  return {0.0, nullptr};
}

int main(int argc, char *argv[]) {
  std::string unmapped_file, land_file;

  po::options_description desc("Allowed options");
desc.add_options()
    ("help,h", "produce help message")
    ("unmapped", po::value<std::string>(&unmapped_file)->required(), "Unmapped geojson file")
    ("land", po::value<std::string>(&land_file)->required(), "Land geojson file")
    ("plan", po::value<std::string>()->required(), "Plan geojson file")
    // ("dst_crs", po::value<std::string>()->required(), "Destination spatial reference system (CRS)")
    ("budget", po::value<double>()->default_value(100000.0), "Budget for path length");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
      std::cout << desc << "\n";
      return 1;
    }

    po::notify(vm);
  } catch (const po::error &e) {
    std::cerr << "Error: " << e.what() << "\n";
    std::cout << desc << "\n";
    return 1;
  }

  GDALAllRegister();

  GDALDataset *ds_unmapped = (GDALDataset *)GDALOpenEx(
      unmapped_file.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr);
  GDALDataset *ds_land = (GDALDataset *)GDALOpenEx(
      land_file.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr);

  if (!ds_unmapped) {
    std::cerr << "Failed to open unmapped file: " << unmapped_file << std::endl;
    return 1;
  }
  if (!ds_land) {
    std::cerr << "Failed to open land file: " << land_file << std::endl;
    GDALClose(ds_unmapped);
    return 1;
  }
  std::vector<OGRPolygon *> unmapped_polygons;
  std::vector<OGRPolygon *> land_polygons;

  double budget = vm["budget"].as<double>();

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
  std::vector<OGRPoint> initial_plan;
  GDALDataset *ds_plan =
      (GDALDataset *)GDALOpenEx(vm["plan"].as<std::string>().c_str(),
                                GDAL_OF_VECTOR, nullptr, nullptr, nullptr);

  if (!ds_plan) {
    std::cerr << "Failed to open plan file: " << vm["plan"].as<std::string>()
              << std::endl;
    return 1;
  }

  OGRLayer *plan_layer = ds_plan->GetLayer(0);
  if (!plan_layer) {
    std::cerr << "Plan file does not contain any layers." << std::endl;
    GDALClose(ds_plan);
    return 1;
  }

  plan_layer->ResetReading();
  OGRFeature *feature = plan_layer->GetNextFeature();
  if (!feature) {
    std::cerr << "Plan file does not contain any features." << std::endl;
    GDALClose(ds_plan);
    return 1;
  }

  OGRGeometry *geom = feature->GetGeometryRef();
  if (!geom || wkbFlatten(geom->getGeometryType()) != wkbLineString) {
    std::cerr << "Plan file does not contain a LineString feature."
              << std::endl;
    OGRFeature::DestroyFeature(feature);
    GDALClose(ds_plan);
    return 1;
  }

  OGRLineString *line = (OGRLineString *)geom;
  for (int i = 0; i < line->getNumPoints(); ++i) {
    OGRPoint point;
    line->getPoint(i, &point);
    initial_plan.push_back(point);
  }

  OGRSpatialReference planSrcSRS;
  const OGRSpatialReference *planSpatialRef = plan_layer->GetSpatialRef();
  if (planSpatialRef) {
    planSrcSRS = *planSpatialRef;
  } else {
    // If the plan has no CRS metadata, assume WGS84 input coordinates.
    planSrcSRS.SetFromUserInput("EPSG:4326");
  }

  // std::cout << "Plan points:" << std::endl;
  // for (const auto &point : initial_plan) {
  //   std::cout << "Point: (" << point.getX() << ", " << point.getY() << ")"
  //             << std::endl;
  // }

  // Define the source (WGS84) and target (Mollweide) spatial references
  OGRLayer *unmapped_layer = ds_unmapped->GetLayer(0);
  if (!unmapped_layer) {
    std::cerr << "Unmapped file does not contain any layers." << std::endl;
    GDALClose(ds_unmapped);
    GDALClose(ds_land);
    return 1;
  }
  const OGRSpatialReference *spatialRef = unmapped_layer->GetSpatialRef();
  if (!spatialRef) {
    std::cerr << "Unmapped file does not have a valid spatial reference."
              << std::endl;
    GDALClose(ds_unmapped);
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
  std::string merc_proj4 = "+proj=webmerc +lon_0=" + std::to_string(center_lon) +
                          " +datum=WGS84 +units=m";
  if (dstSRS.SetFromUserInput(merc_proj4.c_str()) != OGRERR_NONE) {
    std::cerr << "Failed to set Mercator spatial reference: " << merc_proj4 << std::endl;
    GDALClose(ds_unmapped);
    GDALClose(ds_land);
    return 1;
  }
  std::cerr << "Using Mercator projection: " << merc_proj4 << std::endl;

  // Create a coordinate transformation
  OGRCoordinateTransformation *transform =
      OGRCreateCoordinateTransformation(&planSrcSRS, &dstSRS);
  if (!transform) {
    std::cerr << "Failed to create coordinate transformation to AEQD." << std::endl;
    return 1;
  }

  // Transform unmapped polygons
  for (auto &polygon : unmapped_polygons) {
    if (polygon->transform(transform) != OGRERR_NONE) {
      std::cerr << "Failed to transform an unmapped polygon to AEQD." << std::endl;
    }
  }

  // order unmapped polygons by area descending
  std::sort(unmapped_polygons.begin(), unmapped_polygons.end(),
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
    return 1;
  }

  // Transform plan points from plan CRS (WGS84 by default) to AEQD.
  for (auto &point : initial_plan) {
    double x = point.getX();
    double y = point.getY();
    if (!planTransform->Transform(1, &x, &y)) {
      std::cerr << "Failed to transform a plan point to AEQD." << std::endl;
      continue;
    }
    point.setX(x);
    point.setY(y);
  }
  auto plan = initial_plan;

  auto initial_score = score(plan, unmapped_polygons);
  std::cerr << "Initial plan score: " << initial_score << std::endl;
  int iternum = 0;

  std::vector<OGRPolygon> used_polys;

  double initial_length = -1;

  while (true) {
    ++iternum;
    // Find intersections between the path and unmapped polygons
    std::vector<IntersectionInfo> intersection_infos;

    OGRLineString line;
    for (const auto& p : simplify(plan)) {
      line.addPoint(&p);    
    }

    if (initial_length < 0) {
      initial_length = line.get_Length();
    }

    double remainder = budget - line.get_Length();
    std::cerr << "C++ line length: " << line.get_Length() << " budget: " << budget << " remainder: " << remainder << std::endl;
    int non_intersecting = 5;

    for (auto polygon : unmapped_polygons) {
      if (std::find(used_polys.begin(), used_polys.end(), *polygon) != used_polys.end()) {
        continue;
      }
      for (size_t i = 0; i < plan.size() - 1; i++) {
        OGRLineString segment;
        segment.addPoint(&plan[i]);
        segment.addPoint(&plan[i + 1]);

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

        OGRPoint *next_point = &plan[std::min(
            (std::size_t)intersection_infos[i].index + 1, plan.size() - 1)];
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
            new_plan1.push_back(plan[j]);
        }
        for (int k = 0; k < alt_paths.first->getNumPoints(); k++) {
            OGRPoint point;
            alt_paths.first->getPoint(k, &point);
            new_plan1.push_back(point);
        }
        for (std::size_t j = intersection_infos[i].index + 2; j < plan.size(); j++) {
            new_plan1.push_back(plan[j]);
        }

        std::vector<OGRPoint> new_plan2;

        for (std::size_t j = 0; j <= intersection_infos[i].index; j++) {
            new_plan2.push_back(plan[j]);
        }
        for (int k = 0; k < alt_paths.second->getNumPoints(); k++) {
            OGRPoint point;
            alt_paths.second->getPoint(k, &point);
            new_plan2.push_back(point);
        }
        for (std::size_t j = intersection_infos[i].index + 2; j < plan.size(); j++) {
            new_plan2.push_back(plan[j]);
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
            std::cerr << "Path length: " << length << " budget: " << budget << std::endl;
            if(length < budget){
                options.push_back(path);
                used_polys_vector.push_back(
                    intersection_infos[i].unmapped_region);
            }
        }
    }
    if (options.empty()) {
      std::cerr << "No more intersections found. Exiting local improvement."
                << std::endl;
      OGRLineString final_line;
      for (const auto &p : simplify(plan)) {
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
      std::cerr << "Budget: " << budget << std::endl << std::endl;      
      CPLFree(wkt);
      break;
    } 
    // std::cerr << "options:\n";
    double best_score = 0.0;
    long best_i = -1;
    long i = 0;
    for (const auto& option : options) {
      double s = score(option, unmapped_polygons);
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
    if (best_score > score(plan, unmapped_polygons)) {
      plan = options[best_i];
      if (used_polys_vector[best_i]) {
        char *wkt = nullptr;
        used_polys_vector[best_i]->exportToWkt(&wkt);
        std::cerr << "Using polygon WKT: " << (wkt ? wkt : "") << std::endl;
        CPLFree(wkt);
      } else {
        std::cerr << "Using polygon: nullptr" << std::endl;
      }
      used_polys.push_back(*used_polys_vector[best_i]);
      std::cerr << "Updated plan with improved score: " << best_score
                 << std::endl;
       // print geojson of plan
      OGRLineString final_line;
      for (const auto &p : plan) {
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
      for (const auto &p : plan) {
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

  OGRFeature::DestroyFeature(feature);
  GDALClose(ds_plan);

  GDALClose(ds_unmapped);
  GDALClose(ds_land);

  return 0;
}