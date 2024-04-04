/*
 *    Copyright (c) 2012-2020,2022-2023
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_SMASH_COLLIDERMODUS_H_
#define SRC_INCLUDE_SMASH_COLLIDERMODUS_H_

#include <cstring>
#include <memory>
#include <string>
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
 * At the moment there are two ways to specify input for initial conditions in
 * the configuration, one of which is deprecated and will be removed in a next
 * release. For the moment, these variables are of type `std::optional<double>`
 * to *allow* for the key duplication consistently. When more types of IC are
 * implemented in the future, this will allow setting only the appropriate
 * parameters.
 */
struct InitialConditionParameters {
  /// Hypersurface proper time in IC
  std::optional<double> proper_time = std::nullopt;
  /// Lower bound for proper time in IC
  std::optional<double> lower_bound = std::nullopt;
  /// Rapidity cut on hypersurface in IC
  std::optional<double> rapidity_cut = std::nullopt;
  /// Transverse momentum cut on hypersurface IC
  std::optional<double> pT_cut = std::nullopt;
};

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
 * The following configuration options are understood: \ref
 * doxypage_input_conf_modi_collider
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
   **/
  explicit ColliderModus(Configuration modus_config,
                         const ExperimentParameters &parameters);
  /**
   * Creates full path string consisting of file_directory and file_name
   * Needed to initialize a customnucleus.
   * \param[in] file_directory is the path to the external file
   * \param[in] file_name is the name of the external file
   */
  std::string custom_file_path(const std::string &file_directory,
                               const std::string &file_name);
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

  /// Time until nuclei have passed through each other
  double nuclei_passing_time() const {
    const double passing_distance =
        projectile_->get_nuclear_radius() + target_->get_nuclear_radius();
    const double passing_time =
        passing_distance /
        std::sqrt(sqrt_s_NN_ * sqrt_s_NN_ /
                      ((2 * nucleon_mass) * (2 * nucleon_mass)) -
                  1);
    return passing_time;
  }
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
  /// \return The Fermi motion type
  FermiMotion fermi_motion() { return fermi_motion_; }
  /// \return whether the modus is collider (which is, yes, trivially true)
  bool is_collider() const { return true; }
  /// \return center of mass energy per nucleon pair
  double sqrt_s_NN() const { return sqrt_s_NN_; }
  /// \return impact parameter of the collision
  double impact_parameter() const { return impact_; }
  /// \return Whether the calculation frame is the fixed target frame
  bool calculation_frame_is_fixed_target() const {
    return frame_ == CalculationFrame::FixedTarget ? true : false;
  }
  /// \return Proper time of the hypersurface for IC
  std::optional<double> proper_time() const {
    return IC_parameters_.proper_time;
  }
  /// \return Lower bound on proper time of the hypersurface for IC
  std::optional<double> lower_bound() const {
    return IC_parameters_.lower_bound;
  }
  /// \return Maximum rapidity for IC
  std::optional<double> rapidity_cut() const {
    return IC_parameters_.rapidity_cut;
  }
  /// \return Maximum transverse momentum for IC
  std::optional<double> pT_cut() const { return IC_parameters_.pT_cut; }

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
  /** Configure Deformed Nucleus
   *
   * Sets up a deformed nucleus object based on the input parameters in the
   * configuration file.
   * \param[in] nucleus_cfg Subset of configuration, projectile or target
   * section.
   * \param[in] ntest Number of test particles
   * \param[in] nucleus_type String 'projectile' or 'target'. To display an
   * appropriate error message.
   * \return Pointer to the created deformed nucleus object.
   * **/
  static std::unique_ptr<DeformedNucleus> create_deformed_nucleus(
      Configuration &nucleus_cfg, const int ntest,
      const std::string &nucleus_type);
  /**
   * Checks if target and projectile are read from the same external file if
   * they are both initialized as a customnucleus. Function is only called if,
   * projectile is customnucleus.
   * /param[in] proj_config Configuration of projectile nucleus
   * /param[in] targ_config Configuration of target nucleus
   **/
  bool same_inputfile(Configuration &proj_config, Configuration &targ_config);
  /**
   * Impact parameter.
   *
   * The nuclei projectile_ and target_ will be shifted along the x-axis
   * so that their centers move on antiparallel lines that are this
   * distance apart from each other.
   */
  double impact_ = 0.;
  /// Whether the reaction plane should be randomized
  bool random_reaction_plane_;
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

  /**
   * Rotate the reaction plane about the angle phi
   *
   * \param[in] phi Angle about which to rotate
   * \param[in] particles Particles, whose position is rotated
   */
  void rotate_reaction_plane(double phi, Particles *particles);

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
  CalculationFrame frame_ = CalculationFrame::CenterOfVelocity;
  /**
   * An option to include Fermi motion ("off", "on", "frozen")
   */
  FermiMotion fermi_motion_ = FermiMotion::Off;
  /**
   * Beam velocity of the projectile
   */
  double velocity_projectile_ = 0.0;
  /**
   * Beam velocity of the target
   */
  double velocity_target_ = 0.0;
  /// Plain Old Data type to hold parameters for initial conditions
  InitialConditionParameters IC_parameters_;

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

#endif  // SRC_INCLUDE_SMASH_COLLIDERMODUS_H_
