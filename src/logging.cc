/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/logging.h"

#include <array>

#include "include/configuration.h"
#include "include/stringfunctions.h"

namespace Smash {

static einhard::LogLevel global_default_loglevel = einhard::ALL;

einhard::LogLevel default_loglevel() {
  return global_default_loglevel;
}

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

template <std::size_t index>
constexpr typename std::enable_if<(index == 0), int>::type
find_longest_logger_name() {
  return 0;
}
template <std::size_t index>
constexpr typename std::enable_if<(index != 0), int>::type
find_longest_logger_name() {
  using LogAreaTag = typename std::remove_reference<decltype(
      std::get<index - 1>(std::declval<LogArea::AreaTuple &>()))>::type;
  return LogAreaTag::textual_length() > find_longest_logger_name<index - 1>()
             ? LogAreaTag::textual_length()
             : find_longest_logger_name<index - 1>();
}

template <std::size_t index, int>
inline typename std::enable_if<(index == 0)>::type create_all_loggers_impl(
    Configuration &) {}  // do nothing to end the recursion

/*!
 * \internal
 * Recurse over the log areas in the LogArea::AreaTuple type. (The recursion is
 * ended via the overload above.)
 *
 * For every entry in the list the corresponding Logger object in
 * global_logger_collection is set up with area name and verbosity.
 */
template <std::size_t index,
          int longest_name = find_longest_logger_name<index>()>
inline typename std::enable_if<(index != 0)>::type create_all_loggers_impl(
    Configuration &config) {
  using LogAreaTag = typename std::remove_reference<decltype(
      std::get<index - 1>(std::declval<LogArea::AreaTuple &>()))>::type;
  static_assert(LogAreaTag::id == index - 1,
                "The order of types in LogArea::AreaTuple does not match the "
                "id values in the LogArea types. Please fix! (see top of "
                "'include/logging.h')");
  auto &logger = global_logger_collection[LogAreaTag::id];
  const auto tmp = fill_both(LogAreaTag::textual(), longest_name);
  logger.setAreaName(tmp);
  logger.setVerbosity(
      config.take({LogAreaTag::textual()}, global_default_loglevel));
  create_all_loggers_impl<index - 1, longest_name>(config);
}

void create_all_loggers(Configuration config) {
  create_all_loggers_impl<std::tuple_size<LogArea::AreaTuple>::value>(config);
}

}  // namespace Smash
