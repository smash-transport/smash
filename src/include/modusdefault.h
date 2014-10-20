/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#ifndef SRC_INCLUDE_MODUSDEFAULT_H_
#define SRC_INCLUDE_MODUSDEFAULT_H_

#include "forwarddeclarations.h"
#include "threevector.h"

#include <stdexcept>

namespace Smash {

/**
 * \ingroup modus
 * Baseclass for Modus classes that provides default function implementations.
 *
 * This is only a base class for actual Modus classes. Meaning there will never
 * be objects, references, or pointers to ModusDefault. Therefore, it does not
 * have - and will never need any virtual functions.
 *
 * The rules for adding functions to this class are as follows:
 * - This class is empty per default.
 * - You can add a function if you have a function that is different in at least
 *   two subclasses and equal in at least two subclasses.
 * - Code that is common to all goes into ExperimentImplementation.
 */
class ModusDefault {
 public:
  // never needs a virtual destructor

  // Missing functions for concrete Modus implementations:
  // float initial_conditions(Particles *particles);

  /** Enforces sensible positions for the particles.
   *
   * Currently, this is only needed for BoxModus; the other Modi do
   * nothing.
   *
   * \see BoxModus::sanity_check
   */
  int sanity_check(Particles * /*p*/) { return 0; }

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
  void propagate(Particles *particles, const ExperimentParameters &parameters,
                                       const OutputsList &);
  /** Calculates baryon 4-current in the computational frame.
   *  \f[j^{\mu} = (\sqrt{2\pi} \sigma )^{-3} \sum_{i=1}^N B_i u^{\mu}_i
   *  exp \left(- \frac{(\vec r -\vec r_i + \frac{\gamma_i^2}{1 + \gamma_i}
   *  \vec \beta_i (\vec \beta_i, \vec r - \vec r_i))^2}{2\sigma^2} \right)\f]
   *  \f[ \rho_B^{Eckart} = \sqrt{j^{\mu} j_{\mu}} \f]
   *
   * \param[in] r Arbitrary space point where baryon 4-current is calculated
   * \param[in] plist List of all particles to be used in \f$j^{\mu}\f$
   *            calculation. If the distance between particle and calculation
   *            point r, \f$ |r-r_i| > 6 \sigma \f$ then particle input
   *            to density will be ignored.
   * \param[in] gs_sigma Width of the gaussian (\f$ \sigma \f$),
   *  which represents particle Wigner density.
   */
  FourVector baryon_jmu(ThreeVector r, const ParticleList &plist,
                                                  double gs_sigma);

  /** Evaluates potential at point r. Potential is always taken in the local
   * Eckart rest frame, but point r is in the computational frame.
   * 
   * \param[in] r Arbitrary space point where potential is calculated
   * \param[in] plist List of all particles to be used in \f$j^{\mu}\f$
   *            calculation. If the distance between particle and calculation
   *            point r, \f$ |r-r_i| > 6 \sigma \f$ then particle input
   *            to density will be ignored.
   * \param[in] gs_sigma Width of the gaussian (\f$ \sigma \f$),
   *  which represents particle Wigner density.
   **/
  double potential(ThreeVector r, const ParticleList &plist, double gs_sigma);
 
  /** Evaluates potential gradient at point r. Potential is always taken in
   * the local Eckart rest frame, but point r is in the computational frame.
   * 
   * \param[in] r Arbitrary space point where potential gradient is calculated
   * \param[in] plist List of all particles to be used in \f$j^{\mu}\f$
   *            calculation. If the distance between particle and calculation
   *            point r, \f$ |r-r_i| > 6 \sigma \f$ then particle input
   *            to density will be ignored.
   * \param[in] gs_sigma Width of the gaussian (\f$ \sigma \f$),
   *  which represents particle Wigner density.
   **/
   ThreeVector potential_gradient(ThreeVector r, const ParticleList &plist,
                                                          double gs_sigma);

  /** \ingroup exception
   *  BadInput is an error to throw if the configuration options are invalid.
   *
   **/
  struct BadInput : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };
  /// \ingroup exception
  /// Thrown when the requested energy is smaller than the masses
  /// of two particles.
  struct InvalidEnergy : public BadInput {
    using BadInput::BadInput;
  };
};

}  // namespace Smash

#endif  // SRC_INCLUDE_MODUSDEFAULT_H_
