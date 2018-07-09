/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/logging.h"

#include <array>

#include "smash/configuration.h"
#include "smash/stringfunctions.h"

namespace smash {

static einhard::LogLevel global_default_loglevel = einhard::ALL;

einhard::LogLevel default_loglevel() { return global_default_loglevel; }

void set_default_loglevel(einhard::LogLevel level) {
  global_default_loglevel = level;
}

/**
 * An array that stores all pre-configured Logger objects. The objects can be
 * accessed via the logger function.
 */
static std::array<einhard::Logger<>, std::tuple_size<LogArea::AreaTuple>::value>
    global_logger_collection;

einhard::Logger<> &retrieve_logger_impl(int id) {
  return global_logger_collection[id];
}

template <int index, int stop = 0>
constexpr typename std::enable_if<(index == stop), int>::type
find_longest_logger_name() {
  using LogAreaTag = typename std::remove_reference<decltype(
      std::get<index>(std::declval<LogArea::AreaTuple &>()))>::type;
  return LogAreaTag::textual_length();
}
template <int index, int stop = 0, int mid = (index + stop) / 2>
constexpr typename std::enable_if<(index > stop), int>::type
find_longest_logger_name() {
  return find_longest_logger_name<index, mid + 1>() >
                 find_longest_logger_name<mid, stop>()
             ? find_longest_logger_name<index, mid + 1>()
             : find_longest_logger_name<mid, stop>();
}

template <std::size_t index, int>
inline typename std::enable_if<(index == 0)>::type create_all_loggers_impl(
    Configuration &) {}  // do nothing to end the recursion

/**
 * \internal
 * Recurse over the log areas in the LogArea::AreaTuple type. (The recursion is
 * ended via the overload above.)
 *
 * For every entry in the list the corresponding Logger object in
 * global_logger_collection is set up with area name and verbosity.
 *
 * \tparam index Recursion index.
 * \tparam longest_name Length of longest log area name.
 * \param[inout] config Configuration object.
 */
template <std::size_t index,
          int longest_name = find_longest_logger_name<index - 1>()>
inline typename std::enable_if<(index != 0)>::type create_all_loggers_impl(
    Configuration &config) {
  using LogAreaTag = typename std::remove_reference<decltype(
      std::get<index - 1>(std::declval<LogArea::AreaTuple &>()))>::type;
  static_assert(LogAreaTag::id == index - 1,
                "The order of types in LogArea::AreaTuple does not match the "
                "id values in the LogArea types. Please fix! (see top of "
                "'include/logging.h')");
  auto &logger = global_logger_collection[LogAreaTag::id];
  const auto tmp = utf8::fill_both(LogAreaTag::textual(), longest_name);
  logger.setAreaName(tmp);
  logger.setVerbosity(
      config.take({LogAreaTag::textual()}, global_default_loglevel));
  create_all_loggers_impl<index - 1, longest_name>(config);
}

void create_all_loggers(Configuration config) {
  create_all_loggers_impl<std::tuple_size<LogArea::AreaTuple>::value>(config);
}

}  // namespace smash
