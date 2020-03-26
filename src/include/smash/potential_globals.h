/*
 *
 *    Copyright (c) 2014-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_POTENTIAL_GLOBALS_H_
#define SRC_INCLUDE_POTENTIAL_GLOBALS_H_

#include "lattice.h"
#include "potentials.h"

namespace smash {

/// Pointer to the skyrme potential on the lattice
extern RectangularLattice<FourVector> *UB_lat_pointer;

/// Pointer to the symmmetry potential on the lattice
extern RectangularLattice<FourVector> *UI3_lat_pointer;

/// Pointer to a Potential class
extern Potentials *pot_pointer;

}  // namespace smash

#endif  // SRC_INCLUDE_POTENTIAL_GLOBALS_H_
