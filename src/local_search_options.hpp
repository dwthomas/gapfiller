#pragma once

#include <boost/program_options.hpp>

#include <iostream>
#include <optional>
#include <string>

namespace local_search_cli {

namespace po = boost::program_options;

struct Options {
  std::string unmapped_points_file;
  std::string unmapped_lines_file;
  std::string unmapped_polygons_file;
  std::string land_file;
  std::string plan_file;
  double budget = 100000.0;
  double gap_bonus_km = 10.0;
  int consider_gaps = 5;
};

inline std::optional<Options> parse_command_line(int argc, char *argv[]) {
  Options options;
  po::options_description description("Allowed options");
  description.add_options()
      ("help,h", "produce help message")
      ("unmapped-points",
       po::value<std::string>(&options.unmapped_points_file)->required(),
       "Unmapped points GeoJSON file")
      ("unmapped-lines",
       po::value<std::string>(&options.unmapped_lines_file)->required(),
       "Unmapped lines GeoJSON file")
      ("unmapped-polygons",
       po::value<std::string>(&options.unmapped_polygons_file)->required(),
       "Unmapped polygons GeoJSON file")
      ("land", po::value<std::string>(&options.land_file)->required(),
       "Land geojson file")
      ("plan", po::value<std::string>(&options.plan_file)->required(),
       "Plan geojson file")
      ("budget", po::value<double>(&options.budget)->default_value(100000.0),
       "Budget for path length")
      ("gap-bonus-km",
       po::value<double>(&options.gap_bonus_km)->default_value(10.0),
       "Gap completion bonus, in km of polygon-equivalent mapping length")
      ("consider-gaps",
       po::value<int>(&options.consider_gaps)->default_value(5),
       "Maximum number of nearest unmapped lines and points to consider per "
       "iteration");

  po::variables_map variables;
  try {
    po::store(po::parse_command_line(argc, argv, description), variables);
    if (variables.count("help")) {
      std::cout << description << "\n";
      return std::nullopt;
    }
    po::notify(variables);
  } catch (const po::error &error) {
    std::cerr << "Error: " << error.what() << "\n";
    std::cout << description << "\n";
    return std::nullopt;
  }

  return options;
}

}  // namespace local_search_cli
