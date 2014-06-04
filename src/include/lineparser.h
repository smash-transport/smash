/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_LINEPARSER_H_
#define SRC_INCLUDE_LINEPARSER_H_

#include "particles.h"
#include <sstream>
#include <string>
#include <vector>

namespace Smash {

namespace {/*{{{*/
std::string trim(const std::string &s) {
  const auto begin = s.find_first_not_of(" \t\n\r");
  if (begin == std::string::npos) {
    return {};
  }
  const auto end = s.find_last_not_of(" \t\n\r");
  return s.substr(begin, end - begin + 1);
}
struct Line {/*{{{*/
  Line() = default;
  Line(int n, std::string &&t) : number(n), text(std::move(t)) {
  }
  int number;
  std::string text;
};/*}}}*/

std::string build_error_string(std::string message, const Line &line) {/*{{{*/
  return message + " (on line " + std::to_string(line.number) + ": \"" +
         line.text + "\")";
}/*}}}*/

/**
 * Helper function for parsing particles.txt and decaymodes.txt.
 *
 * This function goes through an input stream line by line and removes
 * comments and empty lines. The remaining lines will be returned as a vector
 * of strings and linenumber pairs (Line).
 *
 * \param input an lvalue reference to an input stream
 */
std::vector<Line> line_parser(const std::string &input) {/*{{{*/
  std::istringstream input_stream(input);
  std::vector<Line> lines;
  lines.reserve(50);

  std::string line;
  int line_number = 0;
  while (std::getline(input_stream, line)) {
    ++line_number;
    const auto hash_pos = line.find('#');
    if (hash_pos != std::string::npos) {
      // found a comment, remove it from the line and look further
      line = line.substr(0, hash_pos);
    }
    if (line.find_first_not_of(" \t") == std::string::npos) {
      // only whitespace (or nothing) on this line. Next, please.
      continue;
    }
    lines.emplace_back(line_number, std::move(line));
    line = std::string();
  }
  return std::move(lines);
}/*}}}*/

void ensure_all_read(std::istream &input, const Line &line) {/*{{{*/
  std::string tmp;
  input >> tmp;
  if (!input.eof()) {
    throw Particles::LoadFailure(
        build_error_string("While loading the Particle data:\nGarbage (" + tmp +
                               ") at the remainder of the line.",
                           line));
  }
}/*}}}*/
}  // unnamed namespace/*}}}*/

}  // namespace Smash

#endif  // SRC_INCLUDE_LINEPARSER_H_
