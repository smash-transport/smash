/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/inputfunctions.h"

#include <sstream>
#include <vector>

#include "smash/logging.h"
#include "smash/stringfunctions.h"

namespace smash {

std::vector<Line> line_parser(const std::string &input) {
  const auto &log = logger<LogArea::InputParser>();
  log.trace() << source_location << input;
  std::istringstream input_stream(input);
  std::vector<Line> lines;
  lines.reserve(50);

  std::string line;
  int line_number = 0;
  while (std::getline(input_stream, line)) {
    ++line_number;
    const auto hash_pos = line.find('#');
    if (hash_pos != std::string::npos) {
      // Found a comment, remove it from the line and look further
      line = line.substr(0, hash_pos);
    }
    if (line.find_first_not_of(" \t") == std::string::npos) {
      // Only whitespace (or nothing) on this line. Next, please.
      continue;
    }
    line = trim(line);
    lines.emplace_back(line_number, std::move(line));
    line = std::string();
  }
  return lines;
}

}  // namespace smash
