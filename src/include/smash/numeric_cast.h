/*
 *
 *    Copyright (c) 2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_NUMERIC_CAST_H_
#define SRC_INCLUDE_SMASH_NUMERIC_CAST_H_

#include <limits>

namespace smash {

/**
 * Function template to perform a safe numeric conversion between types.
 * \tparam D Destination type
 * \tparam S Source type
 * \param[in] value Input value to be converted
 * \return The input value converted to the new type \c D .
 * \throw std::overflow_error If the value to be converted cannot be represented
 * by any value of the destination type.
 *
 * \b NOTE: This template is only declared and not implemented! Only
 * specializations are implemented so that only given conversions do not result
 * in compilation errors.
 */
template <typename D, typename S>
constexpr D numeric_cast(const S value);

/**
 * Function template specialization to perform a safe numeric conversion from
 * \c size_t to \c uint32_t .
 * \param[in] value Input value to be converted
 * \return The input value converted to the \c uint32_t type.
 * \throw std::overflow_error
 * If the value to be converted is larger than the maximum of \c uint32_t type.
 */
template <>
constexpr uint32_t numeric_cast(const size_t value) {
  // Conversion from unsigned to unsigned
  //  => the only possible problem is a positive overflow.
  if (value > std::numeric_limits<uint32_t>::max()) {
    throw std::overflow_error(
        "Input value overflows the target uint32_t type!");
  }
  return static_cast<uint32_t>(value);
}

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_NUMERIC_CAST_H_
