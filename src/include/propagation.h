/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_PROPAGATION_H_
#define SRC_INCLUDE_PROPAGATION_H_

#include "experimentparameters.h"
#include "particles.h"
#include "potentials.h"
#include "lattice.h"

namespace Smash {

double calc_hubble (double time);

/** Propagates the positions of all particles on a straight line
  * through the current time step.
  *
  * For each particle, the position is shifted:
  * \f[\vec x^\prime = \vec x + \vec v \cdot \Delta t\f]
  * where \f$\vec x\f$ is the current position, \f$\vec v\f$ its
  * velocity and \f$\Delta t\f$ the duration of this timestep.
  *
  * \param[in,out] particles The particle list in the event
  * \param[in] parameters parameters for the experiment
  */
void propagate_straight_line(Particles *particles,
                             const ExperimentParameters &parameters);

/**
 * Propagates the positions and momenta of all particles through the current
 * time step, according to the equations of motion.
 *
 * For each particle, the position is shifted:
 * \f[\vec x^\prime = \vec x + \vec v \cdot \Delta t\f]
 * where \f$\vec x\f$ is the current position, \f$\vec v\f$ its
 * velocity and \f$\Delta t\f$ the duration of this timestep.
 *
 * The following equations of motion are solved:
 * \f[ \frac{dr}{dt} = p/E \f]
 * \f[ \frac{dp}{dt} = -dU(r)/dr \f]
 *
 * \param[in,out] particles The particle list in the event
 * \param parameters Parameters for the experiment
 * \param pot The potentials in the system
 * \param UB_grad_lat Lattice for Skyrme potential gradient
 * \param UI3_grad_lat Lattice for symmetry potential gradient
 */
void propagate(Particles *particles, const ExperimentParameters &parameters,
               const Potentials &pot,
               RectangularLattice<ThreeVector>* UB_grad_lat,
               RectangularLattice<ThreeVector>* UI3_grad_lat);

}  // namespace Smash
#endif  // SRC_INCLUDE_PROPAGATION_H_
