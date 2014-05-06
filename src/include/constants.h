/*
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_CONSTANTS_H_
#define SRC_INCLUDE_CONSTANTS_H_

#include <cmath>

namespace Smash {

/* GeV <-> fm conversion factor */
const double hbarc = 0.197327053;
/* mb <-> fm^2 conversion factor */
const double fm2_mb = 0.1;
/* Numerical error tolerance */
const double really_small = 1.0e-6;
const double twopi = 2. * M_PI;

}  // namespace Smash

#endif  // SRC_INCLUDE_CONSTANTS_H_
