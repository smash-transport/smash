/*
 *    Copyright (c) 2013-2021
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_SMASH_CONSTANTS_H_
#define SRC_INCLUDE_SMASH_CONSTANTS_H_

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

/// MeV to GeV conversion factor.
constexpr double mev_to_gev = 1.e-3;

/// Numerical error tolerance.
constexpr double really_small = 1.0e-6;

/// A very small double, used to avoid division by zero
constexpr double very_small_double = 1.0e-15;

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
/**
 * Deuteron mass in GeV.
 *
 * Note that this should be the same as in particles.txt.
 */
constexpr double deuteron_mass = 1.8756;

/// Fine-struture constant, approximately 1/137.
constexpr double fine_structure = 7.2973525698e-3;

/// Elementary electric charge in natural units, approximately 0.3
const double elementary_charge = std::sqrt(fine_structure * 4 * M_PI);

/**
 * The maximum value of the random seed used in PYTHIA.
 */
constexpr int maximum_rndm_seed_in_pythia = 900000000;

/**
 * Energy in GeV, below which hard reactions via pythia are impossible.
 * This constraint is technical and comes from the pythia model itself.
 * At the same time, physics-wise, hard cross-sections at the low
 * energies are so small, that this constrant is well justified.
 */
constexpr double minimum_sqrts_pythia_can_handle = 10.0;  // GeV

/**
 * Process ID for any photon process.
 *
 * It is chosen such that it will not conflict with any other process.
 */
constexpr std::uint32_t ID_PROCESS_PHOTON =
    std::numeric_limits<std::uint32_t>::max();

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_CONSTANTS_H_
