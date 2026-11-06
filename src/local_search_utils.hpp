#pragma once

#include <ogr_geometry.h>
#include <ogrsf_frmts.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

#include <GeographicLib/Geodesic.hpp>

namespace local_search_utils {

struct IntersectionInfo {
  OGRLineString *intersection;
  OGRLineString *path_segment;
  OGRPolygon *unmapped_region;
  std::size_t index;
};

inline const GeographicLib::Geodesic &geod =
    GeographicLib::Geodesic::WGS84();

inline std::pair<double, double>
calculateCentroid(const std::vector<OGRPoint> &points) {
  if (points.empty()) {
    return {0.0, 0.0};
  }

  double sum_lat = 0.0;
  double sum_lon_x = 0.0;
  double sum_lon_y = 0.0;

  for (const auto &point : points) {
    const double lat = point.getY();
    const double lon = point.getX();
    sum_lat += lat;

    const double lon_rad = lon * M_PI / 180.0;
    sum_lon_x += std::cos(lon_rad);
    sum_lon_y += std::sin(lon_rad);
  }

  const double n = static_cast<double>(points.size());
  const double avg_lat = sum_lat / n;
  const double avg_lon_rad =
      std::atan2(sum_lon_y / n, sum_lon_x / n);
  const double avg_lon = avg_lon_rad * 180.0 / M_PI;

  return {avg_lat, avg_lon};
}

inline OGRPoint *closestPointOnPolygon(const OGRPolygon &polygon,
                                      const OGRPoint &point) {
  OGRPoint *closestPoint = new OGRPoint();
  double minDist = std::numeric_limits<double>::infinity();

  const OGRLinearRing *ring = polygon.getExteriorRing();
  if (!ring) {
    std::cerr << "Polygon does not have an exterior ring." << std::endl;
    return nullptr;
  }

  const int numPoints = ring->getNumPoints();
  for (int i = 0; i < numPoints - 1; ++i) {
    OGRPoint p1, p2;
    ring->getPoint(i, &p1);
    ring->getPoint(i + 1, &p2);

    OGRLineString edge;
    edge.addPoint(&p1);
    edge.addPoint(&p2);

    OGRPoint projectedPoint;
    const double dist_along = edge.Project(&point);
    const double edgeLen = edge.get_Length();
    if (edgeLen > 0.0) {
      double t = dist_along / edgeLen;
      t = std::max(0.0, std::min(1.0, t));
      projectedPoint.setX(p1.getX() + t * (p2.getX() - p1.getX()));
      projectedPoint.setY(p1.getY() + t * (p2.getY() - p1.getY()));
    } else {
      projectedPoint.setX(p1.getX());
      projectedPoint.setY(p1.getY());
    }

    const double dist = projectedPoint.Distance(&point);
    if (dist < minDist) {
      minDist = dist;
      closestPoint->setX(projectedPoint.getX());
      closestPoint->setY(projectedPoint.getY());
    }
  }

  return closestPoint;
}

inline std::pair<OGRPolygon, int> add_to_polygon(const OGRPolygon &polygon,
                                                const OGRPoint &point) {
  std::vector<OGRPoint> pgon_points;
  int min_i = 0;
  double min_dist = std::numeric_limits<double>::infinity();
  const int n = polygon.getExteriorRing()->getNumPoints();
  for (int i = 0; i < n; i++) {
    OGRPoint p1, p2;
    polygon.getExteriorRing()->getPoint(i, &p1);
    pgon_points.push_back(p1);
    polygon.getExteriorRing()->getPoint((i + 1) % n, &p2);
    OGRLineString edge;
    edge.addPoint(&p1);
    edge.addPoint(&p2);
    const double dist = edge.Distance(&point);
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

inline std::vector<OGRPoint>
transform_line(const std::vector<OGRPoint> &line,
               OGRCoordinateTransformation *transform) {
  std::vector<OGRPoint> transformed_line;
  for (const auto &point : line) {
    double x = point.getX();
    double y = point.getY();
    if (!transform->Transform(1, &x, &y)) {
      std::cerr << "Failed to transform a plan point to AEQD." << std::endl;
      continue;
    }
    transformed_line.emplace_back(x, y);
  }
  return transformed_line;
}

inline std::vector<OGRPoint>
add_line_gap(const std::vector<OGRPoint> &plan, std::size_t insertion_index,
             const OGRLineString &gap_line) {
  std::vector<OGRPoint> successor(plan.begin(),
                                  plan.begin() + insertion_index + 1);
  OGRPoint start, end;
  gap_line.StartPoint(&start);
  gap_line.EndPoint(&end);

  const bool reverse = start.Distance(&plan[insertion_index]) >
                       end.Distance(&plan[insertion_index]);
  if (reverse) {
    for (int i = gap_line.getNumPoints() - 1; i >= 0; --i) {
      OGRPoint point;
      gap_line.getPoint(i, &point);
      successor.push_back(point);
    }
  } else {
    for (int i = 0; i < gap_line.getNumPoints(); ++i) {
      OGRPoint point;
      gap_line.getPoint(i, &point);
      successor.push_back(point);
    }
  }
  successor.insert(successor.end(), plan.begin() + insertion_index + 1,
                   plan.end());
  return successor;
}

inline double true_length(
    OGRLineString *line,
    OGRCoordinateTransformation *inverse_transform) {
  std::vector<OGRPoint> line_points;
  for (int i = 0; i < line->getNumPoints(); i++) {
    OGRPoint p1;
    line->getPoint(i, &p1);
    line_points.push_back(p1);
  }

  const auto wgs84_line = transform_line(line_points, inverse_transform);

  double length = 0.0;
  if (line->getNumPoints() < 2) {
    return length;
  }
  for (std::size_t i = 0; i < wgs84_line.size() - 1; ++i) {
    const OGRPoint &p1 = wgs84_line[i];
    const OGRPoint &p2 = wgs84_line[i + 1];
    const double lat1 = p1.getY();
    const double lon1 = p1.getX();
    const double lat2 = p2.getY();
    const double lon2 = p2.getX();

    double distance_m;
    geod.Inverse(lat1, lon1, lat2, lon2, distance_m);
    length += distance_m;
  }
  return length;
}

inline double score(
    const std::vector<OGRPoint> &metric_path,
    const std::vector<OGRPolygon *> &metric_unmapped_polygons,
    const std::vector<double> &unmapped_polygon_beam_widths,
    OGRCoordinateTransformation *inverse_transform) {
  double sc = 0.0;
  for (std::size_t pi = 0; pi < metric_unmapped_polygons.size(); pi++) {
    auto polygon = metric_unmapped_polygons[pi];
    auto pwidth = unmapped_polygon_beam_widths[pi];
    for (std::size_t i = 0; i < metric_path.size() - 1; i++) {
      OGRLineString segment;
      segment.addPoint(&metric_path[i]);
      segment.addPoint(&metric_path[i + 1]);

      if (polygon->Intersects(&segment)) {
        OGRGeometry *intersection = polygon->Intersection(&segment);
        if (intersection != nullptr &&
            wkbFlatten(intersection->getGeometryType()) == wkbLineString) {
          auto *intersection_line =
              static_cast<OGRLineString *>(intersection);
          sc += true_length(intersection_line, inverse_transform) * pwidth;
        }
        OGRGeometryFactory::destroyGeometry(intersection);
      }
    }
  }
  return sc;
}

inline std::vector<OGRPoint> simplify(const std::vector<OGRPoint> &path) {
  std::unordered_map<double, std::unordered_map<double, std::size_t>> last_seen;

  for (std::size_t i = 0; i < path.size(); i++) {
    const double x = path[i].getX();
    const double y = path[i].getY();
    last_seen[x][y] = i;
  }

  std::vector<OGRPoint> simplified_path;
  std::size_t i = 0;
  while (i < path.size()) {
    const double x = path[i].getX();
    const double y = path[i].getY();
    simplified_path.push_back(path[i]);
    i = last_seen[x][y] + 1;
  }
  return simplified_path;
}

inline std::pair<OGRLineString *, OGRLineString *>
find_alternative_paths(const OGRPolygon &polygon, const OGRPoint &source,
                       const OGRPoint &target, const OGRPoint &next_point) {
  auto x = add_to_polygon(polygon, target);
  x = add_to_polygon(x.first, source);
  const int source_index = x.second;
  int min_index = 0;
  double min_dist = std::numeric_limits<double>::infinity();
  for (int i = 0; i < x.first.getExteriorRing()->getNumPoints(); i++) {
    OGRPoint p;
    x.first.getExteriorRing()->getPoint(i, &p);
    const double dist = p.Distance(&next_point);
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
  while (i != min_index) {
    OGRPoint p;
    x.first.getExteriorRing()->getPoint(i, &p);
    line2.push_back(p);
    i = (i - 1);
    if (i < 0) {
      i = x.first.getExteriorRing()->getNumPoints() - 1;
    }
  }
  line2.push_back(next_point);
  auto *ls1 = new OGRLineString();
  for (const auto &p : line1) {
    ls1->addPoint(&p);
  }

  auto *ls2 = new OGRLineString();
  for (const auto &p : line2) {
    ls2->addPoint(&p);
  }

  return std::make_pair(ls1, ls2);
}

inline std::pair<double, OGRLineString *>
local_improvement(const IntersectionInfo &info) {
  OGRPoint start_point, end_point;
  info.intersection->StartPoint(&start_point);
  info.intersection->EndPoint(&end_point);

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

  const int num_points = ring->getNumPoints();
  bool start_found = false;
  for (int i = 0; i < num_points - 1; ++i) {
    OGRPoint p1, p2;
    ring->getPoint(i, &p1);
    ring->getPoint(i + 1, &p2);

    OGRLineString edge;
    edge.addPoint(&p1);
    edge.addPoint(&p2);

    if (start_point.Within(&edge)) {
      start_found = true;
      line1->addPoint(&start_point);
      for (int j = i + 1; j < num_points; j++) {
        OGRPoint next_point;
        ring->getPoint(j, &next_point);
        line1->addPoint(&next_point);
        if (next_point.Equals(&end_point)) {
          break;
        }
      }

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

}  // namespace local_search_utils
