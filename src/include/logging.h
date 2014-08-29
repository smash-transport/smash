/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_LOGGING_H_
#define SRC_INCLUDE_LOGGING_H_

#include <einhard.hpp>
#include <tuple>
#include <yaml-cpp/yaml.h>
#include <stdexcept>

#include "macros.h"

namespace Smash {
class Configuration;

/**
 * Declares the necessary interface to identify a new log area.
 */
#define DECLARE_LOGAREA(id__, name__)                                     \
  struct name__ {                                                         \
    static constexpr int id = id__;                                       \
    static constexpr const char *textual() { return #name__; }            \
    static constexpr int textual_length() { return sizeof(#name__) - 1; } \
  }

/**
 * The namespace where log areas are declared.
 *
 * To add a new area add one more line with DECLARE_LOGAREA at the bottom: Pick
 * the next number for the id and a name to identify it in the log and source
 * code. Then add the name to the end of the AreaTuple.
 */
namespace LogArea {
DECLARE_LOGAREA( 0, Main);
DECLARE_LOGAREA( 1, Experiment);
DECLARE_LOGAREA( 2, Box);
DECLARE_LOGAREA( 3, Collider);
DECLARE_LOGAREA( 4, Nucleus);
DECLARE_LOGAREA( 5, Sphere);
DECLARE_LOGAREA( 6, Action);
DECLARE_LOGAREA( 7, InputParser);
DECLARE_LOGAREA( 8, ParticleType);
DECLARE_LOGAREA( 9, FindScatter);
DECLARE_LOGAREA(10, Legacy);
DECLARE_LOGAREA(11, Clock);
DECLARE_LOGAREA(12, DecayModes);
DECLARE_LOGAREA(13, Resonances);

/// This type collects all existing log areas so they will be created with the
/// correct log level automatically.
using AreaTuple = std::tuple<Main, Experiment, Box, Collider, Nucleus, Sphere,
                             Action, InputParser, ParticleType, FindScatter,
                             Legacy, Clock, DecayModes, Resonances>;
}  // namespace LogArea

/**
 * Called from main() right after the Configuration object is fully set up to
 * create all logger objects (as defined by LogArea::AreaTuple) with the correct
 * area names and log levels.
 *
 * \param config A configuration object with the log area names as toplevel
 *               keys.
 */
void create_all_loggers(Configuration config);

/** \internal
 * Returns the einhard::Logger object created for the area with the associated
 * index \p id.
 */
einhard::Logger<> &retrieve_logger_impl(int id);

/**
 * Returns the einhard::Logger object created for the named area (see the LogArea types).
 */
template <typename LogAreaTag>
inline einhard::Logger<> &logger() {
  static_assert(LogAreaTag::id < std::tuple_size<LogArea::AreaTuple>::value &&
                    LogAreaTag::id >= 0,
                "The LogArea::AreaTuple is out of sync with the declared log "
                "areas. Please fix! (see top of 'include/logging.h')");
  return retrieve_logger_impl(LogAreaTag::id);
}

/**
 * Hackery that is required to output the location in the source code where the
 * log statement occurs.
 */
#define source_location \
  __FILE__ ":" + std::to_string(__LINE__) + " (" + __func__ + ')'

/**
 * Return the default log level to use if no specific level is configured.
 */
einhard::LogLevel default_loglevel();

/**
 * Set the default log level (what will be returned from subsequent
 * default_loglevel calls).
 *
 * \param level The new log level. See einhard::LogLevel.
 */
void set_default_loglevel(einhard::LogLevel level);

/**
 * Formatting helper
 */
template <typename T>
struct FormattingHelper {
  const T &value;
  const int width;
  const int precision;
  const char *const unit;
  friend std::ostream &operator<<(std::ostream &out,
                                  const FormattingHelper &h) {
    if (h.width > 0) {
      out << std::setfill(' ') << std::setw(h.width);
    }
    if (h.precision >= 0) {
        out << std::setprecision(h.precision);
    }
    out << h.value;
    if (h.unit) {
      out << ' ' << h.unit;
    }
    return out;
  }
};
template <typename T>
FormattingHelper<T> format(const T &value, const char *unit = 0, int width = -1,
                           int precision = -1) {
  return {value, width, precision, unit};
}
}  // namespace Smash

namespace YAML {
/** \internal
 * Enables YAML-cpp to auto-convert a YAML Node to and from an einhard::LogLevel.
 */
template <>
struct convert<einhard::LogLevel> {
  static Node encode(const einhard::LogLevel &x) {
    return Node{einhard::getLogLevelString(x)};
  }
  static bool decode(const Node &node, einhard::LogLevel &x) {
    if (!node.IsScalar()) {
      return false;
    } else {
      x = einhard::getLogLevel(node.Scalar());
      return true;
    }
  }
};
}  // namespace YAML

#endif  // SRC_INCLUDE_LOGGING_H_
