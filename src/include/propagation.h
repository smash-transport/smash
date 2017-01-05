/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_PROPAGATION_H_
#define SRC_INCLUDE_PROPAGATION_H_

#include "particles.h"
#include "potentials.h"
#include "lattice.h"

namespace Smash {

/** Propagates the positions of all particles on a straight line
  * through the current time step.
  *
  * For each particle, the position is shifted:
  * \f[\vec x^\prime = \vec x + \vec v \cdot \Delta t\f]
  * where \f$\vec x\f$ is the current position, \f$\vec v\f$ its
  * velocity and \f$\Delta t\f$ the duration of this timestep.
  *
  * \param[in,out] particles The particle list in the event
  * \param[in] t1 final time
  * \param[out] dt time interval of propagation
  */
double propagate_straight_line(Particles *particles, double t1);

/**
 * Updates the momenta of all particles at the current
 * time step according to the equations of motion:
 *
 * \f[ \frac{dp}{dt} = -dU(r)/dr \f]
 *
 * \param[in,out] particles The particle list in the event
 * \param[in] dt timestep
 * \param pot The potentials in the system
 * \param UB_grad_lat Lattice for Skyrme potential gradient
 * \param UI3_grad_lat Lattice for symmetry potential gradient
 */
void update_momenta(Particles *particles, double dt,
                    const Potentials &pot,
                    RectangularLattice<ThreeVector>* UB_grad_lat,
                    RectangularLattice<ThreeVector>* UI3_grad_lat);

}  // namespace Smash
#endif  // SRC_INCLUDE_PROPAGATION_H_
