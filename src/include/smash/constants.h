/*
 *    Copyright (c) 2013-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_CONSTANTS_H_
#define SRC_INCLUDE_CONSTANTS_H_

#include <cmath>
#include <cstdint>
#include <limits>

/**
 * \file
 *
 * Collection of useful constants that are known at compile time.
 */

namespace smash {

/**
 * GeV <-> fm conversion factor.
 */
constexpr double hbarc = 0.197327053;

/// mb <-> fm^2 conversion factor.
constexpr double fm2_mb = 0.1;

/// GeV^-2 <-> mb conversion factor.
constexpr double gev2_mb = hbarc * hbarc / fm2_mb;

/// Numerical error tolerance.
constexpr double really_small = 1.0e-6;

/**
 * \f$ 2\pi \f$.
 */
constexpr double twopi = 2. * M_PI;

/// Ground state density of symmetric nuclear matter [fm^-3].
constexpr double nuclear_density = 0.168;

/// Physical error tolerance.
constexpr double small_number = 1.0e-4;

/**
 * Nucleon mass in GeV.
 *
 * Note that this should be the same as in particles.txt.
 */
constexpr double nucleon_mass = 0.938;

/**
 * Pion mass in GeV.
 *
 * Note that this should be the same as in particles.txt.
 */
constexpr double pion_mass = 0.138;

/**
 * Kaon mass in GeV.
 *
 * Note that this should be the same as in particles.txt.
 */
constexpr double kaon_mass = 0.494;

/**
 * omega mass in GeV.
 *
 * Note that this should be the same as in particles.txt.
 */
constexpr double omega_mass = 0.783;

/**
 * a1 mass in GeV.
 *
 * Note that this should be the same as in particles.txt.
 */
constexpr double a1_mass = 1.26;
/**
 * Delta mass in GeV.
 *
 * Note that this should be the same as in particles.txt.
 */
constexpr double delta_mass = 1.232;

/// Fine-struture constant, approximately 1/137.
constexpr double fine_structure = 7.2973525698e-3;

/**
 * The maximal cross section (in mb) for which it is guaranteed that all
 * collisions with this cross section will be found.
 *
 * This means that all particle pairs, where the transverse distance is smaller
 * or equal to \f$ \sqrt{200mb/\pi} \f$, will be checked for collions.
 *
 * This maximum occurs in the Delta peak of the pi+p cross section.
 * The only exception of physical cross sections going above 200 mb are the
 * elastic NN and KN cross sections, which diverge at threshold.
 */
constexpr double maximum_cross_section = 200.;  // mb

/**
 * Process ID for any photon process.
 *
 * It is chosen such that it will not conflict with any other process.
 */
constexpr std::uint32_t ID_PROCESS_PHOTON =
    std::numeric_limits<std::uint32_t>::max();

}  // namespace smash

#endif  // SRC_INCLUDE_CONSTANTS_H_
