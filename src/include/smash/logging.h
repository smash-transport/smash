/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_LOGGING_H_
#define SRC_INCLUDE_LOGGING_H_

#include <stdexcept>
#include <tuple>

#include <yaml-cpp/yaml.h>  // NOLINT(build/include_order)
#include <einhard.hpp>

#include "macros.h"

namespace smash {

/// Configuration object to set the verbosity of each area independently.
class Configuration;

/*!\Userguide
 * \page input_logging_ Logging
 * This section controls the logging levels for different areas of the code.
 * It contains the following keys, all of which have string values and are
 * optional:
 * \li \key default (determines the default logging level for all areas)
 * \li \key Legacy
 * \li \key Main
 * \li \key Experiment
 * \li \key Box
 * \li \key Collider
 * \li \key Sphere
 * \li \key Action
 * \li \key InputParser
 * \li \key ParticleType
 * \li \key FindScatter
 * \li \key Clock
 * \li \key DecayModes
 * \li \key Resonances
 * \li \key ScatterAction
 * \li \key Distributions
 * \li \key Propagation
 * \li \key Grid
 * \li \key List
 * \li \key Nucleus
 * \li \key Density
 * \li \key PauliBlocking
 * \li \key Tmn
 * \li \key Fpe
 * \li \key Lattice
 * \li \key AdaptiveTS
 * \li \key Pythia
 * \li \key GrandcanThermalizer
 * \li \key CrossSections
 *
 * Each key can have one of the following values:
 * \li ALL   - Log all messages (default)
 * \li TRACE - The lowest severity for messages describing the program flow
 * \li DEBUG - Debug messages
 * \li INFO  - Messages of informational nature
 * \li WARN  - Warning messages
 * \li ERROR - Non-fatal errors
 * \li FATAL - Messages that indicate terminal application failure
 * \li OFF   - If selected no messages will be output
 *
 * Note that the logging levels TRACE and DEBUG are only available in
 * debugging builds (i.e. with -DCMAKE_BUILD_TYPE=Debug). \n
 *
 * \n
 * Example: Configuring the Logging Area
 * --------------
 * To activate all informational messages (\key INFO) for all logging areas,
 * include the following in the configuration file:
 *\verbatim
 Logging:
     default: INFO
     Legacy: INFO
     Main: INFO
     Experiment: INFO
     Box: INFO
     Collider: INFO
     Sphere: INFO
     ParticleType: INFO
     Action: INFO
     FindScatter: INFO
     DecayModes: INFO
     ScatterAction: INFO
     Pythia: INFO
     InputParser: INFO
     Clock: INFO
     Resonances: INFO
     Distributions: INFO
     Propagation: INFO
     Grid: INFO
     List: INFO
     Nucleus: INFO
     Density: INFO
     PauliBlocking: INFO
     Tmn: INFO
     Fpe: INFO
     Lattice: INFO
     AdaptiveTS: INFO
     GrandcanThermalizer: INFO
     CrossSections: INFO
 \endverbatim
 */

/** \addtogroup logging
 * @{
 *
 * \brief The interfaces in this group are used for debug and informational
 *        console output.
 *
 * SMASH uses the \ref einhard logging library for debug and info/warn/error
 * output to stdout. This library builds upon the C++ <a
 * href="http://en.cppreference.com/w/cpp/io/basic_ostream">ostream</a> classes
 * and thus uses stream operators for converting objects into a text
 * representation.
 *
 * The \ref einhard library supports named output streams (which simply means
 * they automatically add the name to the prefix string). We use this feature to
 * define log areas in SMASH that can be configured independently. The \ref
 * einhard::Logger class supports two options: colorization and verbosity. For
 * colorization we stay with the default of auto-detecting a color-terminal. For
 * verbosity we use a Configuration object to set the verbosity of each area
 * independently. This way we have control over the amount of debug output at
 * runtime and without the need to touch the code/recompile.
 *
 * To output something from your code do the following:
 * \code
 * const auto &log = logger<LogArea::Action>();
 * log.trace(source_location);
 * log.debug("particle: ", p);
 * log.warn("something is fishy");
 * \endcode
 *
 * The einhard::Logger class supports two ways of writing to an output stream:
 * Use stream operators or pass the list of objects for output as parameters.
 * Thus \code
 * log.debug("particle: ", p);
 * \endcode and \code
 * log.debug() << "particle: " << p;
 * \endcode are equivalent (except for a small optimization opportunity in the
 * former variant, that could make it slightly more efficient). You can see,
 * though, that the former variant is more concise and often much easier to type
 * than the stream operators.
 */

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
DECLARE_LOGAREA(0, Main);
DECLARE_LOGAREA(1, Experiment);
DECLARE_LOGAREA(2, Box);
DECLARE_LOGAREA(3, Collider);
DECLARE_LOGAREA(4, Sphere);
DECLARE_LOGAREA(5, Action);
DECLARE_LOGAREA(6, InputParser);
DECLARE_LOGAREA(7, ParticleType);
DECLARE_LOGAREA(8, FindScatter);
DECLARE_LOGAREA(9, Clock);
DECLARE_LOGAREA(10, DecayModes);
DECLARE_LOGAREA(11, Resonances);
DECLARE_LOGAREA(12, ScatterAction);
DECLARE_LOGAREA(13, Distributions);
DECLARE_LOGAREA(14, Propagation);
DECLARE_LOGAREA(15, Grid);
DECLARE_LOGAREA(16, List);  // ListModus
DECLARE_LOGAREA(17, Nucleus);
DECLARE_LOGAREA(18, Density);
DECLARE_LOGAREA(19, PauliBlocking);
DECLARE_LOGAREA(20, Tmn);
DECLARE_LOGAREA(21, Fpe);
DECLARE_LOGAREA(22, Lattice);
DECLARE_LOGAREA(23, Sampling);
DECLARE_LOGAREA(24, AdaptiveTS);
DECLARE_LOGAREA(25, Pythia);
DECLARE_LOGAREA(26, GrandcanThermalizer);
DECLARE_LOGAREA(27, CrossSections);

/**
 * This type collects all existing log areas so they will be created with the
 * correct log level automatically.
 */
using AreaTuple =
    std::tuple<Main, Experiment, Box, Collider, Sphere, Action, InputParser,
               ParticleType, FindScatter, Clock, DecayModes, Resonances,
               ScatterAction, Distributions, Propagation, Grid, List, Nucleus,
               Density, PauliBlocking, Tmn, Fpe, Lattice, Sampling, AdaptiveTS,
               Pythia, GrandcanThermalizer, CrossSections>;
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
 * \return The einhard::Logger object created for the area with the associated
 * index \p id.
 */
einhard::Logger<> &retrieve_logger_impl(int id);

/**
 * \return The einhard::Logger object created for the named area (see the
 * LogArea types).
 * \tparam LogAreaTag Type determining the logging area.
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
 * \return The default log level to use if no specific level is configured.
 */
einhard::LogLevel default_loglevel();

/**
 * Set the default log level (what will be returned from subsequent
 * default_loglevel calls).
 *
 * \param level The new log level. See einhard::LogLevel.
 */
void set_default_loglevel(einhard::LogLevel level);

/** \internal
 * Formatting helper.
 *
 * \tparam T Value that is being formatted.
 */
template <typename T>
struct FormattingHelper {
  /// Value that is being formatted
  const T &value;
  /// Output width
  const int width;
  /// Precision that value is being formatted with
  const int precision;
  /// Unit that is attached at the end of value
  const char *const unit;
  /**
   * Nicely formatted output.
   * \param out Output stream
   * \param h FormattingHelper with given output parameters.
   */
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

/**
 * Acts as a stream modifier for std::ostream to output an object with an
 * optional suffix string and with a given field width and precision.
 *
 * \tparam T Value that is being formatted.
 * \param value The object to be written to the stream.
 * \param unit  An optional suffix string, typically used for a unit. May be
 *              nullptr.
 * \param width The field width to use for \p value.
 * \param precision The precision to use for \p value.
 */
template <typename T>
FormattingHelper<T> format(const T &value, const char *unit, int width = -1,
                           int precision = -1) {
  return {value, width, precision, unit};
}
}  // namespace smash

namespace YAML {
/** \internal
 * Enables YAML-cpp to auto-convert a YAML Node to and from an
 * einhard::LogLevel.
 */
template <>
struct convert<einhard::LogLevel> {
  /**
   * Convert from einhard::LogLevel to YAML::Node.
   *
   * \param x Log level.
   * \return Corresponding YAML node.
   */
  static Node encode(const einhard::LogLevel &x) {
    return Node{einhard::getLogLevelString(x)};
  }

  /**
   * Convert from YAML::Node to einhard::LogLevel.
   *
   * \param[in] node YAML node.
   * \param[out] x Where the corresponding log level will be stored if the
   * conversion was successful.
   * \return Whether the conversion was successful.
   */
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

// @}

#endif  // SRC_INCLUDE_LOGGING_H_
