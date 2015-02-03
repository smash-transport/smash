/*
 *    Copyright (c) 2013-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_CONSTANTS_H_
#define SRC_INCLUDE_CONSTANTS_H_

#include <cmath>

namespace Smash {

/* GeV <-> fm conversion factor */
constexpr double hbarc = 0.197327053;
/* mb <-> fm^2 conversion factor */
constexpr double fm2_mb = 0.1;
/* Numerical error tolerance */
constexpr double really_small = 1.0e-6;
constexpr double twopi = 2. * M_PI;
/* Ground state density of symmetric nuclear matter, fm^-3*/
constexpr double rho0 = 0.168;

}  // namespace Smash

#endif  // SRC_INCLUDE_CONSTANTS_H_
