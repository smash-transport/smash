/*
 *    Copyright (c) 2015-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_PROPAGATION_H_
#define SRC_INCLUDE_PROPAGATION_H_

#include <vector>

#include "lattice.h"
#include "particles.h"
#include "potentials.h"

namespace Smash {

struct ExpansionProperties {
  // Defines the metric to be used
  ExpansionMode mode_;
  // Defines the expansion parameter (faster expansion for larger values)
  double b_;

  ExpansionProperties(ExpansionMode mode, double b) : mode_(mode), b_(b) {}
};

double calc_hubble(double time, const ExpansionProperties &metric);

/** Propagates the positions of all particles on a straight line
  * through the current time step.
  *
  * For each particle, the position is shifted:
  * \f[\vec x^\prime = \vec x + \vec v \cdot \Delta t\f]
  * where \f$\vec x\f$ is the current position, \f$\vec v\f$ its
  * velocity and \f$\Delta t\f$ the duration of this timestep.
  *
  * \param[in,out] particles The particle list in the event
  * \param[in] to_time final time
  * \param[in] beam_momentum This vector of 4-momenta should have
  *            non-zero size only if "frozen Fermi motion" is on.
  *            The the Fermi momenta are only used for collisions,
  *            but not for propagation. In this case beam_momentum
  *            is used for propagation.
  * \return dt time interval of propagation
  */
double propagate_straight_line(Particles *particles, double to_time,
                               const std::vector<FourVector> &beam_momentum);

/** Modifies positions and momentum of all particles to account for
  * space-time deformation.
  */
void expand_space_time(Particles *particles,
                       const ExperimentParameters &parameters,
                       const ExpansionProperties &metric);

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
void update_momenta(Particles *particles, double dt, const Potentials &pot,
                    RectangularLattice<ThreeVector> *UB_grad_lat,
                    RectangularLattice<ThreeVector> *UI3_grad_lat);

}  // namespace Smash
#endif  // SRC_INCLUDE_PROPAGATION_H_
