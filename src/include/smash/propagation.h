/*
 *    Copyright (c) 2015-2018
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

namespace smash {

/**
 * Struct containing the type of the metric and the expansion parameter of
 * the metric. These elements shall be used in the Expansion Mode, which is
 * implemented in SMASH to compare with the analytical solution
 * \iref{Bazow:2015dha} to the Boltzmann equation with a Hubble expansion.
 */
struct ExpansionProperties {
  /// Type of metric used
  ExpansionMode mode_;
  /// Expansion parameter in the metric (faster expansion for larger values)
  double b_;
  /**
   * Constructor of ExpansionProperties
   *
   * \param[in] mode Type of metric used in the Expansion Mode.
   * \param[in] b Expansion parameter in the metric
   */
  ExpansionProperties(ExpansionMode mode, double b) : mode_(mode), b_(b) {}
};

/**
 * Calculate the Hubble parameter \f$H(t)\f$, which describes how large
 * the expansion flow is. The flow \f$\vec v=H(t) \vec x\f$
 * \iref{Tindall:2016try}
 *
 * \param[in] time time in the computational frame. [fm]
 * \param[in] metric Struct containing the parameters needed to
 *            calculate the metric.
 * \return Hubble parameter [fm^\f${-1}\f$]
 */
double calc_hubble(double time, const ExpansionProperties &metric);

/**
 * Propagates the positions of all particles on a straight line
 * to a given moment.
 *
 * For each particle, the position is shifted:
 * \f[ \vec x^\prime = \vec x + \vec v \Delta t \f]
 * where \f$\vec x\f$ is the current position, \f$\vec v\f$ its
 * velocity and \f$\Delta t\f$ the duration of this timestep.
 *
 * \param[out] particles The particle list in the event
 * \param[in] to_time final time [fm]
 * \param[in] beam_momentum This vector of 4-momenta should have
 *            non-zero size only if "frozen Fermi motion" is on.
 *            The the Fermi momenta are only used for collisions,
 *            but not for propagation. In this case beam_momentum
 *            is used for propagating the initial nucleons. [GeV]
 * \return dt time interval of propagation which is equal to the
 *            difference between the final time and the initial
 *            time read from the 4-position of the particle.
 */
double propagate_straight_line(Particles *particles, double to_time,
                               const std::vector<FourVector> &beam_momentum);

/**
 * Modifies positions and momentum of all particles to account for
 * space-time deformation.
 *
 * \param[out] particles All the particles in the event
 * \param[in] parameters A struct containing the parameters from which
 *            we extract the time in the computational frame.
 * \param[in] metric A struct containing the parameters need to calculate
 *            the metric
 */
void expand_space_time(Particles *particles,
                       const ExperimentParameters &parameters,
                       const ExpansionProperties &metric);

/**
 * Updates the momenta of all particles at the current
 * time step according to the equations of motion:
 *
 * \f[ \frac{dp}{dt} = \vec E + \vec v \times \vec B \f]
 *
 * \param[out] particles The particle list in the event
 * \param[in] dt timestep
 * \param[in] pot The potentials in the system
 * \param[in] FB_lat Lattice for the electric and magnetic
 *            components of the Skyrme force
 * \param[in] FI3_lat Lattice for the electric and magnetic
 *            components of the symmetry force
 */
void update_momenta(
    Particles *particles, double dt, const Potentials &pot,
    RectangularLattice<std::pair<ThreeVector, ThreeVector>> *FB_lat,
    RectangularLattice<std::pair<ThreeVector, ThreeVector>> *FI3_lat);

}  // namespace smash
#endif  // SRC_INCLUDE_PROPAGATION_H_
