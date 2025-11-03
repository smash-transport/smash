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
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "smash/particledata.h"
#include "smash/particles.h"

namespace smash {

/**
 * Structure to convert a given value into ASCII format, such that all methods
 * return a \c std::string.
 */
struct ToASCII {
  /// Return type of this converter.
  using type = std::string;

  /// Character used to separate fields in ASCII output.
  static constexpr std::optional<char> separator{' '};

  /// Character used to mark the end of a record in ASCII output.
  static constexpr std::optional<char> end_of_line{'\n'};

  /**
   * Converts an integer.
   *
   * \param[in] value number to convert
   */
  type as_integer(int value) const { return std::to_string(value); }

  /**
   * Converts a double with 6 digits of precision.
   *
   * \note The usage of \c std::snprintf over \c std::ostringstream is because
   * of performance reasons (in C++20 this will be replaced with \c std::format
   * which is even better). The returned string is constructed from the buffer
   * in a way to exclude the terminating null character from the buffer. The
   * hard-coded buffer size should fit any number, but a couple of assert are
   * used to possibly investigate unexpected behaviour.
   *
   * \warning Since the buffer size is needed twice, it makes sense to store it
   * in a variable. However, cpplint complains if the variable name is not
   * starting with \c k followed by CamelCase.
   *
   * \param[in] value number to convert
   */
  type as_double(double value) {
    constexpr size_t kBufferSize = 13;
    char buffer[kBufferSize];
    const auto length = std::snprintf(buffer, kBufferSize, "%g", value);
    assert(static_cast<size_t>(length) < kBufferSize);
    assert(length > 0);
    return std::string{buffer, buffer + length};
  }

  /**
   * Converts a double with 9 digits of precision.
   *
   * \see \c as_double for further information.
   *
   * \note The duplication of the code of the \c as_double method is done on
   * purpose as naively extracting a function passing the \c std::snprintf
   * format string as parameter would trigger a warning in compilation (the
   * format has to be a literal in order to be checked by the compiler at
   * compile time) and the effort to avoid this is not worth now, especially
   * since this code will be changed anyhow when using C++20.
   *
   * \param[in] value number to convert
   */
  type as_precise_double(double value) {
    constexpr size_t kBufferSize = 16;
    char buffer[kBufferSize];
    const auto length = std::snprintf(buffer, kBufferSize, "%.9g", value);
    assert(static_cast<size_t>(length) < kBufferSize);
    assert(length > 0);
    return std::string{buffer, buffer + length};
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
 * Structure to convert a given value into binary format, such that all methods
 * return a \c std::vector<char>.
 */
class ToBinary {
 public:
  /// Return type of this converter.
  using type = std::vector<char>;

  /// separator char is not needed for Binary format
  static constexpr std::optional<char> separator = std::nullopt;

  /// Endline char is not needed for Binary format
  static constexpr std::optional<char> end_of_line = std::nullopt;

  /**
   * Converts an integer to binary format.
   *
   * \param[in] value number to convert
   * \return a vector of char representing the binary format of the integer
   */
  type as_integer(int value) const { return as_binary_data(value); }

  /**
   * Converts a double to binary format.
   *
   * \param[in] value number to convert
   * \return a vector of char representing the binary format of the double
   */
  type as_double(double value) const { return as_binary_data(value); }

  /**
   * Converts a double to binary format, intended for precise representation.
   * Note that for binary output there is no difference between this and the
   * \c ToBinary::as_double method, but this method has still to be introduced
   * to allow other classes to get the converter class as a template parameter.
   *
   * \param[in] value number to convert
   * \return a vector of char representing the binary format of the double
   */
  type as_precise_double(double value) const { return as_double(value); }

  /**
   * Converts a string to binary format.
   *
   * \param[in] str string to convert
   * \return a vector of char representing the binary format of the string
   */
  type as_string(const std::string& str) const {
    return type(str.begin(), str.end());
  }

 private:
  /**
   * Template method to convert numbers into binary format.
   *
   * \param[in] value number to convert
   * \return a vector of char representing the binary format of the number
   */
  template <typename T>
  type as_binary_data(T value) const {
    type binary_data(sizeof(T));
    std::memcpy(binary_data.data(), &value, sizeof(T));
    return binary_data;
  }
};

/**
 * A general-purpose formatter for output, supporting both ASCII and binary
 * formats.
 *
 * This class allows the output of particle data in a flexible and configurable
 * manner, either in human-readable ASCII format or compact binary format. It
 * uses a template parameter `Converter` to determine the desired output format,
 * which must conform to the interface of either \c ToASCII or \c ToBinary.
 *
 * New quantities can be added for output by:
 * 1. Adding their corresponding getter to the constructor, which extracts the
 *    value from a \c ParticleData instance.
 * 2. Adding the proper key-value pair to the \c units_ map, specifying the unit
 *    of the quantity.
 *
 * \tparam Converter The desired output format. At the moment it must be either
 *                   \c ToASCII or `ToBinary`.
 */
template <typename Converter,
          std::enable_if_t<std::is_same_v<Converter, ToASCII> ||
                               std::is_same_v<Converter, ToBinary>,
                           bool> = true>
class OutputFormatter {
 public:
  /**
   * Creates the formatter. This reduces the number of literal strings flying
   * around in the codebase, and since this is called only once per output file,
   * in the beginning of the run, there is almost no efficiency lost compared to
   * having fixed strings.
   *
   * \param[in] in_quantities list of quantities to be output.
   *
   * \throw std::invalid_argument if the list of quantities is empty
   * \throw std::invalid_argument if unknown quantities exist in the list
   * \throw std::invalid_argument if incompatible quantities exist in the list
   * \throw std::invalid_argument if there are repeated quantities
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
          return this->converter_.as_integer(in.pdgcode().get_decimal());
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
          return this->converter_.as_integer(in.get_history().p1.get_decimal());
        });
      } else if (quantity == "pdg_mother2") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_integer(in.get_history().p2.get_decimal());
        });
      } else if (quantity == "baryon_number") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_integer(in.pdgcode().baryon_number());
        });
      } else if (quantity == "strangeness") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_integer(in.pdgcode().strangeness());
        });
      } else if (quantity == "spin0") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_double(in.spin_vector()[0]);
        });
      } else if (quantity == "spinx") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_double(in.spin_vector()[1]);
        });
      } else if (quantity == "spiny") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_double(in.spin_vector()[2]);
        });
      } else if (quantity == "spinz") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_double(in.spin_vector()[3]);
        });
      } else if (quantity == "0") {  // for OSCAR1999
        getters_.push_back([this]([[maybe_unused]] const ParticleData& in) {
          return this->converter_.as_integer(0);
        });
      } else if (quantity == "tau") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_double(in.hyperbolic_time());
        });
      } else if (quantity == "eta" || quantity == "eta_s") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_double(in.spatial_rapidity());
        });
      } else if (quantity == "mt") {
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_double(in.transverse_mass());
        });
      } else if (quantity == "Rap" || quantity == "y_rap") {
        // "Rap" is used for compatibility with vHLLE
        getters_.push_back([this](const ParticleData& in) {
          return this->converter_.as_double(in.rapidity());
        });
      }
    }
  }

  /**
   * Computes and returns the total size of the formatted representation
   * of a single particle using all registered getters.
   *
   * The total size consists of:
   *  - the sum of the sizes of all field strings (from each getter),
   *  - (n - 1) separator characters if a separator is defined,
   *    since separators appear only between fields, not after the last one,
   *  - one end-of-line character if an end-of-line marker is defined.
   *
   * \param[in] sample  A representative particle used to compute the size.
   * \return Total size of the formatted representation of one particle.
   */
  std::size_t compute_single_size(const ParticleData& sample) const {
    const std::size_t n = getters_.size();
    if (n == 0) {
      return 0;
    }
    std::size_t size = 0;
    for (const auto& getter : getters_) {
      const typename Converter::type tmp = getter(sample);
      size += tmp.size();
    }

    if (Converter::separator.has_value()) {
      size += (n - 1);
    }

    if (Converter::end_of_line.has_value()) {
      size += 1;
    }

    return size;
  }
  /**
   * Produces a line of data representing a single particle
   * suitable for writing to an output file.
   *
   * \param[in] p Particle whose information is to be written.
   * \return A buffer containing the formatted data.
   *
   * \see append_to_buffer
   */
  typename Converter::type single_particle_data(const ParticleData& p) const {
    typename Converter::type chunk{};
    append_to_buffer(p, chunk);
    return chunk;
  }

  /**
   * Builds multiple particle lines for a range of particles by appending the
   * Converter::type representation of each particle into a single buffer.
   *
   * The resulting buffer reflects exactly what the Converter defines(
   * per-record separators and end-of-line markers, if any). No additional
   * formatting is applied beyond what the Converter specifies.
   *
   * \tparam Range Container type — constrained to `Particles` or
   * `ParticleList`.
   * \param[in] particles Container of particles whose data should be
   * formatted and added to the buffer
   * \return Converter::type A buffer containing the concatenated
   * representation of all particles in the block.
   *
   * \note Capacity is reserved using a size estimate based on the first
   * element.
   *
   * \see append_to_buffer
   *
   * \return Converter::type
   */
  template <class Range,
            std::enable_if_t<std::is_same_v<Range, Particles> ||
                                 std::is_same_v<Range, ParticleList>,
                             bool> = true>
  typename Converter::type particles_data_chunk(const Range& particles) const {
    typename Converter::type chunk{};
    if (particles.size() == 0)
      return chunk;

    chunk.reserve(particles.size() * compute_single_size(particles.front()));

    for (const ParticleData& p : particles) {
      append_to_buffer(p, chunk);
    }
    return chunk;
  }

  /**
   * Produces the line with quantities for the header of the output file.
   *
   * Appends each quantity name to the output buffer, separating them
   * by the separator character if it is defined in the Converter. End of line
   * character is added last if defined.
   *
   * \return Converter::type
   */
  typename Converter::type quantities_line() const {
    typename Converter::type out{};

    // Educated size guess to reduce the number of reallocations in push_backs:
    out.reserve(quantities_.size() * 5);

    for (const auto& string_name : quantities_) {
      if constexpr (Converter::separator) {
        if (!out.empty()) {
          out.push_back(*Converter::separator);
        }
      }
      typename Converter::type name = converter_.as_string(string_name);
      out.insert(out.end(), name.begin(), name.end());
    }
    if constexpr (Converter::end_of_line) {
      out.push_back(*Converter::end_of_line);
    }
    return out;
  }

  /**
   * Produces the line with units for the header of the output file.
   *
   * Looks up each quantity’s unit in the units_ map and appends it to the
   * output buffer. Each unit is separated by the separator character if
   * defined. End of line character is added last if defined.
   *
   * \return Converter::type
   */
  typename Converter::type unit_line() const {
    typename Converter::type out{};

    // Educated size guess to reduce the number of reallocations in push_backs:
    out.reserve(quantities_.size() * 5);

    for (const auto& key : quantities_) {
      if constexpr (Converter::separator) {
        if (!out.empty()) {
          out.push_back(*Converter::separator);
        }
      }
      const auto& unit = converter_.as_string(units_.at(key));
      out.insert(out.end(), unit.begin(), unit.end());
    }
    if constexpr (Converter::end_of_line) {
      out.push_back(*Converter::end_of_line);
    }
    return out;
  }

 private:
  /// Member to convert data into the correct output format.
  Converter converter_{};

  /// List of quantities to be written.
  std::vector<std::string> quantities_{};

  /// List of getters to extract output data from the `ParticleData` object.
  std::vector<std::function<typename Converter::type(const ParticleData&)>>
      getters_{};

  /// Map with known quantities and corresponding units.
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
      {"0", "0"},  // for OSCAR1999;
      {"tau", "fm"},
      {"eta", "none"},
      {"eta_s", "none"},
      {"mt", "GeV"},
      {"Rap", "none"},
      {"y_rap", "none"},
      {"spin0", "none"},
      {"spinx", "none"},
      {"spiny", "none"},
      {"spinz", "none"}};

  /// Checks whether the quantities requested are known and unique
  void validate_quantities() {
    if (quantities_.empty()) {
      throw std::invalid_argument(
          "OutputFormatter: Empty quantities handed over to the class.");
    }
    std::string error_message{};
    std::string repeated{};
    for (const std::string& quantity : quantities_) {
      if (std::count(quantities_.begin(), quantities_.end(), quantity) > 1) {
        repeated += "'" + quantity + "',";
      }
    }
    if (!repeated.empty()) {
      error_message += "Repeated \"Quantities\": " + repeated +
                       " please fix the configuration file.\n";
    }
    std::string unknown{};
    for (const std::string& quantity : quantities_) {
      if (units_.count(quantity) == 0) {
        unknown += "'" + quantity + "',";
      }
    }
    if (!unknown.empty()) {
      error_message += "Unknown \"Quantities\": " + unknown +
                       " please fix the configuration file.\n";
    }
    if (!repeated.empty() || !unknown.empty())
      throw std::invalid_argument(error_message);

    const bool oscar1999_id_is_given =
        std::find(quantities_.begin(), quantities_.end(), "id") !=
        quantities_.end();
    const bool oscar2013_id_is_given =
        std::find(quantities_.begin(), quantities_.end(), "ID") !=
        quantities_.end();
    if (oscar1999_id_is_given && oscar2013_id_is_given) {
      throw std::invalid_argument(
          "Both 'id' and 'ID' cannot be provided in the \"Quantities\" key "
          "together. Please, fix the configuration file.");
    }
  }

  /**
   * Appends the Converter::type representation of a single particle to the
   * provided buffer.
   *
   * The buffer is used exactly as provided by the caller; no assumptions or
   * modifications are made to its initial state or capacity.
   *
   * \param[in] particle Particle whose information is to be appended.
   * \param[out] buffer Destination buffer to which the converted data is
   * appended.
   */

  void append_to_buffer(const ParticleData& particle,
                        typename Converter::type& buffer) const {
    for (std::size_t i = 0; i < getters_.size(); ++i) {
      if constexpr (Converter::separator) {
        if (i > 0)
          buffer.push_back(*Converter::separator);
      }
      auto data = getters_[i](particle);
      buffer.insert(buffer.end(), std::make_move_iterator(data.begin()),
                    std::make_move_iterator(data.end()));
    }
    if constexpr (Converter::end_of_line) {
      buffer.push_back(*Converter::end_of_line);
    }
  }
};
namespace details {
/**
 * \brief Writes particle data in multiple chunks if the total buffer size
 *        exceeds a predefined maximum.
 *
 * This method avoids creating a single excessively large binary buffer when
 * writing many particles at once. Instead, it splits the write into several
 * smaller chunks. This can prevent excessive memory usage and improve
 * stability on systems or filesystems that may have trouble with very large
 * write calls.
 *
 * The maximum buffer size is currently set to 1 GB (10^9 bytes), but this can
 * be adapted in the future if needed. If the total data size of the particle
 * block is below this threshold, the method simply delegates to the
 * `write` function which should perform a single write call.
 *
 * Otherwise, the data is accumulated particle by particle until the buffer
 * reaches the threshold. The buffer is then given to the write function which
 * should flush to disk.
 *
 * \note This utility does not strictly belong to this file. At the moment,
 *       the objects using it do not have a clean hierarchy that would allow
 *       both OscarOutput and BinaryOutput to inherit a shared implementation,
 *       hence it lives here temporarily.
 *
 * \todo Once the hierarchy is cleaned up, move this into a common base class
 *       that OscarOutput and BinaryOutput inherit from.
 *
 * \tparam Converter Converter used by OutputFormatter to produce the
 *         buffer type (must define Converter::type).
 * \tparam Range Container type — enforced to be either `Particles` or
 *         `ParticleList`.
 * \param[in] particles Container of particles whose particle_line
 *            representation is to be written.
 * \param[in] formatter Formatter responsible for converting particles into
 *            the corresponding Converter::type buffer representation.
 * \param[in] write Callable that receives each filled buffer chunk and
 *            performs the actual write to the underlying output (e.g. file).
 * \param[in] max_buffer_bytes Maximum buffer size in bytes before the
 *            accumulated data is flushed via \p write (default: 1'000'000'000).
 */

template <typename Converter, class Range,
          std::enable_if_t<std::is_same_v<Range, Particles> ||
                               std::is_same_v<Range, ParticleList>,
                           bool> = true>
void write_in_chunk_impl(
    const Range& particles, const OutputFormatter<Converter>& formatter,
    std::function<void(const typename Converter::type&)> write,
    std::size_t max_buffer_bytes = 1'000'000'000) {
  if (particles.size() == 0)
    return;

  const std::size_t bytes_per_particle =
      formatter.compute_single_size(particles.front());

  if (2.0 * bytes_per_particle > max_buffer_bytes) {
    throw std::runtime_error(
        "write_in_chunk_impl: the estimated size of a single particle line "
        "exceeds half of the configured max_buffer_bytes.\n"
        "This effectively means only one particle would fit per chunk, "
        "which defeats the purpose of chunked writing.\n"
        "Increase max_buffer_bytes to at least twice the particle line size "
        "to use this function correctly.");
  }

  if (particles.size() * bytes_per_particle <= max_buffer_bytes) {
    write(formatter.particles_data_chunk(particles));
    return;
  }

  using Buffer = typename Converter::type;
  Buffer buffer;
  buffer.reserve(max_buffer_bytes);
  std::size_t current_size = 0;

  for (const auto& particle : particles) {
    Buffer line = formatter.single_particle_data(particle);
    const std::size_t line_size = line.size();
    if (current_size + line_size > max_buffer_bytes) {
      write(buffer);
      buffer.clear();
      current_size = 0;
    }
    buffer.insert(buffer.end(), std::make_move_iterator(line.begin()),
                  std::make_move_iterator(line.end()));
    current_size += line_size;
  }

  if (!buffer.empty()) {
    write(buffer);
  }
}
}  // namespace details

/**
 * \brief User-facing wrapper for chunked particle writing.
 *
 * Forwards the call to the internal implementation.
 *
 * \see details::write_in_chunk_impl
 */
template <typename Converter, class Range,
          std::enable_if_t<std::is_same_v<Range, Particles> ||
                               std::is_same_v<Range, ParticleList>,
                           bool> = true>
void write_in_chunk(
    const Range& particles, const OutputFormatter<Converter>& formatter,
    std::function<void(const typename Converter::type&)> write) {
  details::write_in_chunk_impl(particles, formatter, write);
}

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_OUTPUTFORMATTER_H_
