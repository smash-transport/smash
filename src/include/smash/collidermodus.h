/*
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_COLLIDERMODUS_H_
#define SRC_INCLUDE_COLLIDERMODUS_H_

#include <cstring>
#include <memory>
#include <utility>

#include "deformednucleus.h"
#include "forwarddeclarations.h"
#include "fourvector.h"
#include "interpolation.h"
#include "modusdefault.h"
#include "nucleus.h"
#include "pdgcode.h"

namespace smash {

struct ExperimentParameters;

/**
 * \ingroup modus
 * ColliderModus: Provides a modus for colliding nuclei.
 *
 * To use this modus, choose
 * \code
 * General:
 *      Modus: Collider
 * \endcode
 * in the configuration file.
 *
 * Options for ColliderModus go in the "Modi"→"Collider" section of the
 * configuration.
 *
 * The following configuration options are understood: \ref input_modi_collider_
 */
class ColliderModus : public ModusDefault {
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
   * \throw ColliderEmpty if projectile or nucleus are empty (i.e. do
   *                      not contain particles)
   * \throw InvalidEnergy if sqrts from config is not large enough to support
   *                      the colliding masses of the nuclei, or if E_kin or
   *                      P_lab are negative
   * \throw domain_error if more or less than exactly one of the
   *                     input energy options is specified, or if custom
   *                     impact parameter Values and Yields are improperly
   *                     supplied
   * \todo JB:remove the second parameter?
   **/
  explicit ColliderModus(Configuration modus_config,
                         const ExperimentParameters &parameters);

  /**
   * Generates initial state of the particles in the system.
   * In particular, it initializes the momenta and positions of nucleons
   * withing the colliding nuclei.
   *
   * \param[out] particles An empty list that gets filled up by this function
   * \param[in] parameters The initialization parameters of the system
   * \return The starting time of the simulation (negative, so that nuclei
   *         collide exactly at t=0)
   * \throw domain_error if the velocities of each nucleus are >= 1, or if
   *                     input for Fermi motion is invalid
   */
  double initial_conditions(Particles *particles,
                            const ExperimentParameters &parameters);
  /// \return The total number of test particles in the initial nuclei
  int total_N_number() const { return target_->size() + projectile_->size(); }
  /// \return The number of test particles in the projectile nucleus
  int proj_N_number() const { return projectile_->size(); }
  /**
   * \return the beam velocity of the projectile, which will be used to
   *         calculate the beam momenta in experiment.cc if Fermi motion is
   *         frozen.
   */
  double velocity_projectile() const { return velocity_projectile_; }
  /**
   * \return the beam velocity of the target, which will be used to calculate
   *         the beam momenta in experiment.cc if Fermi motion is frozen.
   */
  double velocity_target() const { return velocity_target_; }
  /**
   * \return A flag: whether to allow first collisions within the same nucleus.
   */
  bool cll_in_nucleus() { return cll_in_nucleus_; }
  /// \return The Fermi motion type
  FermiMotion fermi_motion() { return fermi_motion_; }
  /// \return whether the modus is collider (which is, yes, trivially true)
  bool is_collider() const { return true; }
  /// \return impact parameter of the collision
  double impact_parameter() const { return impact_; }
  /**
   * \ingroup exception
   *  Thrown when either \a projectile_ or \a target_ nuclei are empty.
   */
  struct ColliderEmpty : public ModusDefault::BadInput {
    using ModusDefault::BadInput::BadInput;
  };

 private:
  /**
   * Projectile.
   *
   * The object that goes from negative z-values to positive z-values
   * with positive velocity.
   **/
  std::unique_ptr<Nucleus> projectile_;
  /**
   * Target.
   *
   * The object that goes from positive z-values to negative z-values
   * with negative velocity. In fixed target experiments, the target is
   * at rest.
   **/
  std::unique_ptr<Nucleus> target_;
  /**
   * Center-of-mass energy squared of the nucleus-nucleus collision.
   *
   * Needs to be double to allow for calculations at LHC energies
   * **/
  double total_s_;
  /**
   * Center-of-mass energy of a nucleon-nucleon collision.
   *
   * Needs to be double to allow for calculations at LHC energies
   * **/
  double sqrt_s_NN_;
  /**
   * Impact parameter.
   *
   * The nuclei projectile_ and target_ will be shifted along the x-axis
   * so that their centers move on antiparallel lines that are this
   * distance apart from each other.
   **/
  double impact_ = 0.;
  /// Method used for sampling of impact parameter.
  Sampling sampling_ = Sampling::Quadratic;
  /// Minimum value of impact parameter.
  double imp_min_ = 0.0;
  /// Maximum value of impact parameter.
  double imp_max_ = 0.0;
  /// Maximum value of yield. Needed for custom impact parameter sampling.
  double yield_max_ = 0.0;
  /// Pointer to the impact parameter interpolation.
  std::unique_ptr<InterpolateDataLinear<double>> impact_interpolation_ =
      nullptr;
  /** Sample impact parameter.
   *
   * Samples the impact parameter from values between imp_min_ and imp_max_, if
   * linear or quadratic sampling is used. By specifying impact parameters and
   * corresponding yields, custom sampling can be used.
   * This depends on the value of sampling_.
   *
   * Note that imp_max_ less than imp_min_ also works fine.
   *
   **/
  void sample_impact();
  /** Initial z-displacement of nuclei.
   *
   * Projectile is shifted on -(this value) in z-direction
   * and target on +(this value)*v_target/v_projectile. In this way
   * projectile and target touch at t=0 in z=0.
   **/
  double initial_z_displacement_ = 2.0;
  /**
   * Reference frame for the system, as specified from config
   */
  /// \todo enum classes do not appear in doxygen, since they are
  /// defined in fowarddeclarations which is intentionally neglected
  /// from doxygen (s. fowarddeclerations.h)
  CalculationFrame frame_ = CalculationFrame::CenterOfVelocity;
  /**
   * An option to include Fermi motion ("off", "on", "frozen")
   */
  FermiMotion fermi_motion_ = FermiMotion::Off;
  /**
   * An option to accept first collisions within the same nucleus
   */
  bool cll_in_nucleus_ = false;
  /**
   * Beam velocity of the projectile
   */
  double velocity_projectile_ = 0.0;
  /**
   * Beam velocity of the target
   */
  double velocity_target_ = 0.0;
  /**
   * Get the frame dependent velocity for each nucleus, using
   * the current reference frame. \see frame_
   *
   * \param[in] mandelstam_s The total center-of-mass energy of the system.
   * \param[in] m_a The (positive) mass of the projectile.
   * \param[in] m_b The (positive) mass of the target.
   * \return A pair < v_a, v_b > containing the velocities of the nuclei.
   * \throw domain_error if the reference frame is not properly specified
   */
  std::pair<double, double> get_velocities(double mandelstam_s, double m_a,
                                           double m_b);

  /**\ingroup logging
   * Writes the initial state for the ColliderModus to the output stream.
   *
   * \param[in] out The ostream into which to output
   * \param[in] m The ColliderModus object to write into out
   */
  friend std::ostream &operator<<(std::ostream &, const ColliderModus &);
};

}  // namespace smash

#endif  // SRC_INCLUDE_COLLIDERMODUS_H_
