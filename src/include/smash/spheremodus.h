/*
 *    Copyright (c) 2013-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_SMASH_SPHEREMODUS_H_
#define SRC_INCLUDE_SMASH_SPHEREMODUS_H_

#include <stdint.h>

#include <cmath>
#include <list>
#include <map>

#include "forwarddeclarations.h"
#include "modusdefault.h"

namespace smash {

/**
 * \ingroup modus
 * SphereModus: Provides a modus for expanding matter calculations
 *
 * Matter is put in a sphere of radius R with uniform density;
 * isotropic thermal momenta are typically used for initialization,
 * although other initial momentum states are also included,
 * see \iref{Bazow:2016oky} and \iref{Tindall:2016try}
 *
 * To use this modus, choose
 * \code
 * General:
 *      Modus: Sphere
 * \endcode
 * in the configuration file.
 *
 * Options for SphereModus go in the "Modi"→"Sphere" section of the
 * configuration:
 *
 * \code
 * Modi:
 *      Sphere:
 *              # definitions here
 * \endcode
 *
 * The following configuration options are understood: \ref
 * doxypage_input_conf_modi_sphere
 */
class SphereModus : public ModusDefault {
 public:
  /**
   * Constructor
   *
   * Takes all there is to take from the (truncated!) configuration
   * object (only contains configuration for this modus).
   *
   * \param[in] modus_config The configuration object that sets all
   *                         initial conditions of the experiment.
   * \param[in] parameters Unused, but necessary because of templated
   *                       initialization
   */
  explicit SphereModus(Configuration modus_config,
                       const ExperimentParameters &parameters);

  /**
   * Generates initial state of the particles in the system according to
   * specified parameters: number of particles of each species, momentum
   * and coordinate space distributions. Susbsequently makes the total
   * 3-momentum 0.
   *
   * \param[out] particles An empty list that gets filled up by this function
   * \param[in] parameters The initialization parameters of the box
   * \return The starting time of the simulation
   */
  double initial_conditions(Particles *particles,
                            const ExperimentParameters &parameters);

  /// \return If the modus is sphere modus, which is always true
  bool is_sphere() const { return true; }
  /// \return radius
  double radius() const { return radius_; }

 private:
  /// Sphere radius (in fm)
  double radius_;
  /// Temperature for momentum distribution (in GeV)
  double sphere_temperature_;
  /// Starting time for the Sphere
  const double start_time_ = 0.;
  /**
   *  Whether to use a thermal initialization for all particles
   *  instead of specific numbers
   */
  const bool use_thermal_ = false;
  /**
   *  Baryon chemical potential for thermal initialization;
   *  only used if \key use_thermal_ is true
   */
  const double mub_;
  /**
   * Strange chemical potential for thermal initialization;
   * only used if \key use_thermal_ is true
   */
  const double mus_;
  /**
   *  Charge chemical potential for thermal initialization;
   *  only used if \key use_thermal_ is true
   */
  const double muq_;
  /**
   * Multiplicative factor for thermal multiplicity of heavy flavored hadrons;
   * only used if \key use_thermal_ is true
   */
  const double hf_multiplier_;
  /**
   * In case of thermal initialization:
   * - true -- account for resonance spectral functions, while computing
   * multiplicities and sampling masses,
   * - false -- simply use pole masses.
   */
  const bool account_for_resonance_widths_;
  /**
   * Particle multiplicities at initialization;
   * required if use_thermal_ is false
   */
  const std::map<PdgCode, int> init_multipl_;
  /**
   * Average multiplicities in case of thermal initialization.
   * Saved to avoid recalculating at every event
   */
  std::map<PdgCode, double> average_multipl_;
  /**
   * Initialization scheme for momenta in the sphere;
   * used for expanding metric setup
   */
  const SphereInitialCondition init_distr_;
  /**
   * Parameter \f$ u_0\f$ in the initial flow velocity profile of particles in
   * the sphere, which has the form \f$ u = u_0 (r / R)^n\f$.
   */
  const double radial_velocity_;
  /**
   * Parameter \f$ n\f$ in the initial flow velocity profile of particles in
   * the sphere, which has the form \f$ u = u_0 (r / R)^n\f$.
   */
  const double radial_velocity_exponent_;
  /**
   * Optional PDG code of the particle to use as a jet, i.e. a single high
   * energy particle at the center (0,0,0) of the expanding sphere. This
   * particle will be placed at t=0 and will initially be moving along the x
   * axis, in the positive or negative direction depending on the sign of
   * its momentum.
   */
  const std::optional<PdgCode> jet_pdg_;
  /**
   * Initial momentum of the jet particle; only used if jet_pdg_ is not nullopt
   */
  const double jet_mom_;
  /**
   * Initial position of the jet particle; only used if jet_pdg_ is not nullopt
   */
  const ThreeVector jet_pos_;
  /**
   * Create the back to back jet with the corresponding antiparticle; only used
   * if jet_pdg_ is not nullopt
   */
  const bool jet_back_;
  /**
   * Initial separation between the back to back jets; can only be set by the
   * user if jet_back_ is true
   */
  const double jet_back_separation_ = smash_NaN<double>;
  /// Spin interaction type
  const SpinInteractionType spin_interaction_type_;
  /**\ingroup logging
   * Writes the initial state for the Sphere to the output stream.
   *
   * \param[in] out The ostream into which to output
   * \param[in] m The SphereModus object to write into out
   */
  friend std::ostream &operator<<(std::ostream &out, const SphereModus &m);
};
}  // namespace smash
#endif  // SRC_INCLUDE_SMASH_SPHEREMODUS_H_
