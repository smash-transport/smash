/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_INPUTFUNCTIONS_H_
#define SRC_INCLUDE_INPUTFUNCTIONS_H_

#include <sstream>
#include <string>
#include <vector>

#include "particles.h"
#include "logging.h"

namespace Smash {

namespace {/*{{{*/
/// takes a string and strips leading and trailing whitespaces.
std::string trim(const std::string &s) {
  const auto begin = s.find_first_not_of(" \t\n\r");
  if (begin == std::string::npos) {
    return {};
  }
  const auto end = s.find_last_not_of(" \t\n\r");
  return s.substr(begin, end - begin + 1);
}
/// Line consists of a line number and the contents of that line
struct Line {/*{{{*/
  /// initialize line with empty string and number
  Line() = default;
  /// initialize a line with line number \p n and text \p t
  Line(int n, std::string &&t) : number(n), text(std::move(t)) {
  }
  /// line number
  int number;
  /// line content.
  std::string text;
};/*}}}*/

/** builds a meaningful error message
 *
 * Takes the message and quotes the Line where the error occurs
 *
 * \param[in] message Error message
 * \param[in] line Line object containing line number and line content.
 */
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

/// makes sure that nothing is left to read from this line.
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

/**
 * Utility function to read a complete input stream (e.g. file) into one string.
 *
 * \param input The input stream. Since it reads until EOF und thus "uses up the
 * whole input stream" the function takes an rvalue reference to the stream
 * object (just pass a temporary).
 *
 * \note There's no slicing here: the actual istream object is a temporary that
 * is not destroyed until read_all returns.
 */
std::string read_all(std::istream &&input) {
  return {std::istreambuf_iterator<char>{input},
          std::istreambuf_iterator<char>{}};
}

}  // unnamed namespace/*}}}*/

}  // namespace Smash

#endif  // SRC_INCLUDE_INPUTFUNCTIONS_H_
