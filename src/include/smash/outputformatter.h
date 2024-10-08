/*
 *
 *    Copyright (c) 2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_OUTPUTFORMATTER_H_
#define SRC_INCLUDE_SMASH_OUTPUTFORMATTER_H_

#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "smash/particledata.h"

namespace smash {

/**
 * Structure to convert a given value into ASCII format, such that all methods
 * return a \c std::string.
 */
struct ToASCII {
  /// Return type of this converter.
  using type = std::string;

  /**
   * Converts an integer.
   *
   * \param[in] value number to convert
   */
  type as_integer(const int value) const { return std::to_string(value); }

  /**
   * Converts a double with 6 digits of precision.
   *
   * \param[in] value number to convert
   */
  type as_double(const double value) const {
    std::ostringstream os{};
    os << std::setprecision(6) << value;
    return os.str();
  }

  /**
   * Converts a double with 9 digits of precision.
   *
   * \param[in] value number to convert
   */
  type as_precise_double(const double value) const {
    std::ostringstream os{};
    os << std::setprecision(9) << value;
    return os.str();
  }

  /**
   * Because %ToASCII converts into strings, this simply returns the string
   * itself.
   *
   * \param[inout] str string to be written
   */
  type as_string(const std::string& str) const { return str; }
};

/**
 * \tparam Converter format desired for the output.
 * A general formatter used for output purposes, which currently only works for
 * ASCII-based formats. In the future, a Binary converter will be implemented.
 *
 * In case further quantities need to be made outputable, one needs to add them
 * in the constructor with the appropriate getter from ParticleData, as well as
 * insert the proper pair in the units map.
 */
template <typename Converter,
          std::enable_if_t<std::is_same_v<Converter, ToASCII>, bool> = true>
class OutputFormatter {
 public:
  /**
   * Creates the formatter. This reduces the number of literal strings flying
   * around in the codebase, and since this is called only once per output file,
   * in the beginning of the run, there is almost no efficiency lost compared to
   * having fixed strings.
   *
   * \param[in] in_quantities list of quantities to be output.
   */
  explicit OutputFormatter(const std::vector<std::string>& in_quantities)
      : quantities_(in_quantities) {
    validate_quantities();
    for (const std::string& quantity : quantities_) {
      if (quantity == "t") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_double(in.position()[0]);
        });
      } else if (quantity == "x") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_double(in.position()[1]);
        });
      } else if (quantity == "y") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_double(in.position()[2]);
        });
      } else if (quantity == "z") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_double(in.position()[3]);
        });
      } else if (quantity == "mass") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_double(in.effective_mass());
        });
      } else if (quantity == "p0") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_precise_double(in.momentum()[0]);
        });
      } else if (quantity == "px") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_precise_double(in.momentum()[1]);
        });
      } else if (quantity == "py") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_precise_double(in.momentum()[2]);
        });
      } else if (quantity == "pz") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_precise_double(in.momentum()[3]);
        });
      } else if (quantity == "pdg") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_string(in.pdgcode().string().c_str());
        });
      } else if (quantity == "ID" || quantity == "id") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_integer(in.id());
        });
      } else if (quantity == "charge") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_integer(in.type().charge());
        });
      } else if (quantity == "ncoll") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_integer(
              in.get_history().collisions_per_particle);
        });
      } else if (quantity == "form_time") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_double(in.formation_time());
        });
      } else if (quantity == "xsecfac") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_double(in.xsec_scaling_factor());
        });
      } else if (quantity == "proc_id_origin") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_integer(in.get_history().id_process);
        });
      } else if (quantity == "proc_type_origin") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_integer(
              static_cast<int>(in.get_history().process_type));
        });
      } else if (quantity == "time_last_coll") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_double(
              in.get_history().time_last_collision);
        });
      } else if (quantity == "pdg_mother1") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_string(
              in.get_history().p1.string().c_str());
        });
      } else if (quantity == "pdg_mother2") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_string(
              in.get_history().p2.string().c_str());
        });
      } else if (quantity == "baryon_number") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_integer(in.pdgcode().baryon_number());
        });
      } else if (quantity == "strangeness") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_integer(in.pdgcode().strangeness());
        });
      } else if (quantity == "spin_projection") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_integer(in.spin_projection());
        });
      } else if (quantity == "0") {  // for OSCAR1999
        getters_.push_back([this]([[maybe_unused]] const ParticleData& in) {
          return this->converter_.as_integer(0);
        });
      }
    }
  }

  /**
   * Produces the line with formatted data for the body of the output file.
   *
   * \param[in] p particle whose information is to be written.
   * \return string with formatted data separated by a space.
   */
  typename Converter::type data_line(const ParticleData& p) const {
    return std::accumulate(
        std::begin(getters_), std::end(getters_), std::string{},
        [&](const std::string& ss,
            const std::function<typename Converter::type(const ParticleData&)>
                getter) {
          return ss.empty() ? getter(p) : ss + " " + getter(p);
        });
  }

  /**
   * Produces the line with quantities for the header of the output file.
   * \return string with name of quantities separated by a space.
   */
  typename Converter::type quantities_line() const {
    return std::accumulate(
        std::begin(quantities_), std::end(quantities_), std::string{},
        [&](const std::string& ss, const std::string& s) {
          return ss.empty() ? converter_.as_string(s)
                            : ss + " " + converter_.as_string(s);
        });
  }

  /**
   * Produces the line with units for the header of the output file.
   * \return string with units separated by a space.
   */
  typename Converter::type unit_line() const {
    return std::accumulate(
        std::begin(quantities_), std::end(quantities_), std::string{},
        [&](const std::string& ss, const std::string& s) {
          return ss.empty() ? converter_.as_string(units_.at(s))
                            : ss + " " + converter_.as_string(units_.at(s));
        });
  }

 private:
  /// Desired format for data conversion. Currently only ASCII is available.
  Converter converter_{};

  /// List of quantities to be written.
  std::vector<std::string> quantities_{};

  /// List of getters for the corresponding quantities.
  std::vector<std::function<typename Converter::type(const ParticleData&)>>
      getters_{};

  ///  Map with known quantities and corresponding units.
  const std::map<std::string, std::string> units_ = {
      {"t", "fm"},
      {"x", "fm"},
      {"y", "fm"},
      {"z", "fm"},
      {"mass", "GeV"},
      {"p0", "GeV"},
      {"px", "GeV"},
      {"py", "GeV"},
      {"pz", "GeV"},
      {"pdg", "none"},
      {"ID", "none"},
      {"id", "none"},  // used in OSCAR1999
      {"charge", "e"},
      {"ncoll", "none"},
      {"form_time", "fm"},
      {"xsecfac", "none"},
      {"proc_id_origin", "none"},
      {"proc_type_origin", "none"},
      {"time_last_coll", "fm"},
      {"pdg_mother1", "none"},
      {"pdg_mother2", "none"},
      {"baryon_number", "none"},
      {"strangeness", "none"},
      {"spin_projection", "none"},
      {"0", "0"}};  // for OSCAR1999

  /// Checks whether the quantities requested are known and unique
  void validate_quantities() {
    if (quantities_.empty()) {
      throw std::invalid_argument(
          "OutputFormatter: Quantities not given, "
          "please fix the configuration file.");
    }
    std::string error_message{};
    std::string repeated{};
    for (const std::string& quantity : quantities_) {
      if (std::count(quantities_.begin(), quantities_.end(), quantity) > 1) {
        repeated += "'" + quantity + "',";
      }
    }
    if (!repeated.empty()) {
      error_message += "OutputFormatter: Repeated quantities: " + repeated +
                       " please fix the configuration file.\n";
    }
    std::string unknown{};
    for (const std::string& quantity : quantities_) {
      if (units_.count(quantity) == 0) {
        unknown += "'" + quantity + "',";
      }
    }
    if (!unknown.empty()) {
      error_message += "OutputFormatter: Unknown quantities: " + unknown +
                       " please fix the configuration file.\n";
    }
    if (!repeated.empty() || !unknown.empty())
      throw std::invalid_argument(error_message);
  }
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_OUTPUTFORMATTER_H_
