/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_ACTION_GLOBALS_H_
#define SRC_INCLUDE_ACTION_GLOBALS_H_

#include "lattice.h"
#include "potentials.h"

namespace smash {

/* This should silence the warnings, but it does not work due to a bug in GCC.
 * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=53431 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

/// Pointer to the skyrme potential on the lattice
extern RectangularLattice<double> *UB_lat_pointer;

/// Pointer to the symmmetry potential on the lattice
extern RectangularLattice<double> *UI3_lat_pointer;

/// Pointer to a Potential class
extern Potentials *pot_pointer;
#pragma GCC diagnostic pop

}  // namespace smash

#endif  // SRC_INCLUDE_ACTION_GLOBALS_H_
