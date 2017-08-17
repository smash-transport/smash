/*  TODO(mkretz): license?

    Copyright (C) 2011-2014 Matthias Kretz <kretz@kde.org>

    Vc is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    Vc is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with Vc.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef SRC_TESTS_ULP_H_
#define SRC_TESTS_ULP_H_

#include <cmath>

#include <limits>

#ifdef _MSC_VER
namespace std {
static inline bool isnan(float x) { return _isnan(x); }
static inline bool isnan(double x) { return _isnan(x); }
}  // namespace std
#endif

template <typename T>
static T ulpDiffToReference(T val, T ref) {
  if (val == ref || (std::isnan(val) && std::isnan(ref))) {
    return 0;
  }
  if (ref == T(0)) {
    return 1 + ulpDiffToReference(std::abs(val), std::numeric_limits<T>::min());
  }
  if (val == T(0)) {
    return 1 + ulpDiffToReference(std::numeric_limits<T>::min(), std::abs(ref));
  }

  int exp;
  /*tmp = */ frexp(ref, &exp);  // ref == tmp * 2 ^ exp => tmp == ref * 2 ^ -exp
  // tmp is now in the range [0.5, 1.0[
  // now we want to know how many times we can fit 2^-numeric_limits<T>::digits
  // between tmp and
  // val * 2 ^ -exp
  return ldexp(std::abs(ref - val), std::numeric_limits<T>::digits - exp);
}

template <typename T>
static T ulpDiffToReferenceSigned(T val, T ref) {
  return ulpDiffToReference(val, ref) * (val < ref ? -1 : 1);
}

#endif  // SRC_TESTS_ULP_H_

// vim: foldmethod=marker
