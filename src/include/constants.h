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

/**
 * GeV <-> fm conversion factor
 *
 * \fpPrecision  This is \c double to make sure to always have sufficient precision.
 */
constexpr double hbarc = 0.197327053;
/**
 * mb <-> fm^2 conversion factor
 */
constexpr float fm2_mb = 0.1;
/**
 * Numerical error tolerance
 */
constexpr float really_small = 1.0e-6;
/**
 * \f$ 2\pi \f$
 *
 * \fpPrecision  This is \c double to make sure to always have sufficient precision.
 */
constexpr double twopi = 2. * M_PI;
/**
 * Ground state density of symmetric nuclear matter, fm^-3
 */
constexpr float nuclear_density = 0.168;
/**
 * nucleon mass in GeV
 *
 * Note that this should be the same as in particles.txt.
 */
constexpr float nucleon_mass = 0.938;
/**
 * kaon mass in GeV
 *
 * Note that this should be the same as in particles.txt.
 */
constexpr float kaon_mass = 0.494;

/**
 * Fine-struture constant, approximately 1/137
 */
constexpr double alpha = 7.2973525698e-3;

/**
 * The maximal cross section (in mb) for which it is guaranteed that all
 * collisions with this cross section will be found.
 *
 * This means that all particle pairs, where the transverse distance is smaller
 * or equal to \f$ \sqrt{200mb/\pi} \f$, will be checked for collions.
 *
 * This maximum occurs in the Delta peak of the pi+p cross section.
 * The only exception of physical cross sections going above 200 mb are the
 * elastic NN cross sections, which diverge at threshold.
 */
constexpr float maximum_cross_section = 200.f;  // mb

}  // namespace Smash

#endif  // SRC_INCLUDE_CONSTANTS_H_
