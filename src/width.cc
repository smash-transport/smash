/*
 *
 *    Copyright (c) 2012-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "include/width.h"

#include <cstdio>
#include <stdexcept>
#include <istream>
#include <gsl/gsl_integration.h>

#include "include/resonances.h"

namespace Smash {

float BlattWeisskopf(const float x, const int L)
#ifdef NDEBUG
    noexcept
#endif
{
  const auto x2 = x * x;
  switch (L) {
    case 0:
      return 1.f;
    case 1:
      return x2 / (1.f + x2);
    case 2: {
      const auto x4 = x2 * x2;
      return x4 / (9.f + 3.f * x2 + x4);
    }
    /* The following lines should be correct. But since nothing in SMASH uses
     * L > 2, this code is untested and dead. Therefore we only keep it as a
     * reference for later.
     * See also input sanitization in load_decaymodes in decaymodes.cc.
    case 3:
      return x4 * x2 / (225.f + 45.f * x2 + 6.f * x4 + x4 * x2);
    case 4:
      return x4 * x4 /
             (11025.f + 1575.f * x2 + 135.f * x4 + 10.f * x2 * x4 + x4 * x4);
    */
#ifndef NDEBUG
    default:
      throw std::invalid_argument(
          std::string("Wrong angular momentum in BlattWeisskopf: ") +
          std::to_string(L));
#endif
  }
  return 0.f;
}


double Post_FF_sqr(double m, double M0, double s0, double L) {
  double FF = (L*L*L*L + (s0-M0*M0)*(s0-M0*M0)/4.) /
              (L*L*L*L + (m*m-(s0+M0*M0)/2.) * (m*m-(s0+M0*M0)/2.));
  return FF*FF;
}


}  // namespace Smash
