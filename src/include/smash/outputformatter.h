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
#include <string>
#include <vector>

#include "smash/particledata.h"

namespace smash {

/// Structure to convert a given value into ASCII
struct ASCII {
  /// Return type of this converter.
  using ReturnType = std::string;

  /**
   * Converts a value into a string with the default format specifier for the
   * value type.
   *
   * \param[in] value Value to be written.
   * \return Formatted string.
   */
  template <typename T>
  std::string operator()(const T& value) const {
    return std::to_string(value);
  }
  /**
   * Overload of the above for when a string is given.
   *
   * \param[in] value string to be written.
   * \return the same string.
   */
  std::string operator()(const std::string& value) const { return value; }
  /**
   * Overload of the above for when a literal string (sequence of chars) is
   * given.
   *
   * \param[in] value string to be written.
   * \return the same string.
   */
  std::string operator()(const char* value) const { return value; }

  /**
   * Converts a value into a string with a format specifier, occupying a maximum
   * of 20 chars.
   *
   * \param[in] value Value to be written.
   * \param[in] format Format specifier to use.
   * \return Formatted string.
   */
  template <typename T>
  std::string operator()(const T& value, const char* format) const {
    // Maximum size of a string is 20.
    constexpr std::size_t kMaxSize = 20;
    char buffer[kMaxSize];
    std::snprintf(buffer, kMaxSize, format, value);
    return buffer;
  }
};

/// Structure to convert a given value into Binary
struct Binary {
  /// Return type of this converter.
  using ReturnType = const void*;

  /**
   * Gives the pointer to the stored data with no format specified.
   *
   * \param[in] value Value to be written.
   * \return a pointer to the value.
   */
  template <typename T>
  const void* operator()(const T& value) const {
    return static_cast<const void*>(&value);
  }
};

/**
 * \tparam Converter format desired for the output.
 * A general formatter used for output purposes, which currently only works for
 * ASCII-based formats.
 */
template <typename Converter>
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
          return this->converter_(in.position()[0], "%g");
        });
      } else if (quantity == "x") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.position()[1], "%g");
        });
      } else if (quantity == "y") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.position()[2], "%g");
        });
      } else if (quantity == "z") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.position()[3], "%g");
        });
      } else if (quantity == "mass") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.effective_mass(), "%g");
        });
      } else if (quantity == "p0") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.momentum()[0], "%.9g");
        });
      } else if (quantity == "px") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.momentum()[1], "%.9g");
        });
      } else if (quantity == "py") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.momentum()[2], "%.9g");
        });
      } else if (quantity == "pz") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.momentum()[3], "%.9g");
        });
      } else if (quantity == "pdg") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.pdgcode().string().c_str(), "%s");
        });
      } else if (quantity == "ID" || quantity == "id") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.id(), "%i");
        });
      } else if (quantity == "charge") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.type().charge(), "%i");
        });
      } else if (quantity == "ncoll") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.get_history().collisions_per_particle,
                                  "%i");
        });
      } else if (quantity == "form_time") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.formation_time(), "%g");
        });
      } else if (quantity == "xsecfac") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.xsec_scaling_factor(), "%g");
        });
      } else if (quantity == "proc_id_origin") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.get_history().id_process, "%i");
        });
      } else if (quantity == "proc_type_origin") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(
              static_cast<int>(in.get_history().process_type), "%i");
        });
      } else if (quantity == "time_last_coll") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.get_history().time_last_collision, "%g");
        });
      } else if (quantity == "pdg_mother1") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.get_history().p1.string().c_str(), "%s");
        });
      } else if (quantity == "pdg_mother2") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.get_history().p2.string().c_str(), "%s");
        });
      } else if (quantity == "baryon_number") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.pdgcode().baryon_number(), "%i");
        });
      } else if (quantity == "strangeness") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.pdgcode().strangeness(), "%i");
        });
      } else if (quantity == "spin_projection") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_(in.spin_projection(), "%i");
        });
      } else if (quantity == "0") {  // for OSCAR1999
        getters_.push_back([this]([[maybe_unused]] const ParticleData& in) {
          return this->converter_(0, "%i");
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
  const std::string data_line(const ParticleData& p) {
    std::string line;
    for (const auto& getter : getters_) {
      line += getter(p) + " ";
    }
    if (!line.empty()) {
      line.pop_back();
      line += "\n";
    }
    return line;
  }

  /**
   * Produces the line with quantities for the header of the output file.
   * \return string with name of quantities separated by a space.
   */
  const std::string quantities_line() {
    std::string header;
    for (const std::string& q : quantities_) {
      header += this->converter_(q) + " ";
    }
    header.pop_back();
    return header;
  }

  /**
   * Produces the line with units for the header of the output file.
   * \return string with units separated by a space.
   */
  std::string unit_line() {
    std::string line;
    for (const std::string& q : quantities_) {
      line += this->converter_(units_.at(q)) + " ";
    }
    line.pop_back();
    return line;
  }

 private:
  /// Desired format for data conversion. Currently only ASCII is available.
  Converter converter_;

  /// List of quantities to be written.
  std::vector<std::string> quantities_;

  /// List of getters for the corresponding quantities.
  std::vector<
      std::function<typename Converter::ReturnType(const ParticleData&)>>
      getters_;

  /**
   *  Map with known quantities and corresponding units.
   */
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
      {"id", "none"},
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
      {"0", "0"}};

  /// Checks whether the quantities requested are known
  void validate_quantities() {
    for (auto& quantity : quantities_) {
      if (units_.count(quantity) == 0) {
        throw std::invalid_argument("OutputFormatter: Unknown quantity: " +
                                    quantity);
      }
    }
  }
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_OUTPUTFORMATTER_H_
