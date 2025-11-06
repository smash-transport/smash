/*
 *
 *    Copyright (c) 2019-2022,2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_BREMSSTRAHLUNGACTION_H_
#define SRC_INCLUDE_SMASH_BREMSSTRAHLUNGACTION_H_

#include <utility>

#include "scatteraction.h"

namespace smash {
/**
 * \ingroup action
 * BremsAction is a special action which takes two incoming particles
 * and performs a perturbative scattering where a Bremsstrahlung photon is
 * produced.
 * The final state particles are not further propagated, only written
 * to the photon output.
 */
class BremsstrahlungAction : public ScatterAction {
 public:
  /**
   * Construct a ScatterActionBrems object.
   *
   * \param[in] in ParticleList of incoming particles.
   * \param[in] time Time relative to underlying hadronic action.
   * \param[in] n_frac_photons Number of photons to produce for each hadronic
   *                            scattering.
   * \param[in] hadronic_cross_section_input Cross-section of
   *                                          underlying hadronic cross-section.
   * \param[in] spin_interaction_type Which type of spin interaction to use.
   * \return The constructed object.
   */

  BremsstrahlungAction(const ParticleList &in, const double time,
                       const int n_frac_photons,
                       const double hadronic_cross_section_input,
                       const SpinInteractionType spin_interaction_type =
                           SpinInteractionType::Off);
  /**
   * Create the final state and write to output.
   *
   * \param[in] outputs List of all outputs. Does not have to be a specific
   *                      photon output, the function will take care of this.
   */
  void perform_bremsstrahlung(const OutputsList &outputs);

  /**
   * Generate the final-state for the Bremsstrahlung process. Generates only
   * 3-body final state.
   */
  void generate_final_state() override;

  /**
   * Sample the final state anisotropically, considering the differential
   * cross sections with respect to theta and k.
   */
  void sample_3body_phasespace();

  /**
   * Return the weight of the last created photon.
   *
   * \return The total weight.
   */
  double get_total_weight() const override { return weight_; }

  /**
   * Return the total cross section of the underlying hadronic scattering
   * It is necessary for the weighting procedure.
   *
   * \return total cross-section [mb]
   */
  double hadronic_cross_section() const { return hadronic_cross_section_; }

  /**
   * Adds one hadronic process with a given cross-section.
   *
   * The intended use is to add the hadronic cross-section from the already
   * performed hadronic
   * action without recomputing it.
   *
   * \param[in] reaction_cross_section Total cross-section of underlying
   *                                    hadronic process [mb]
   */
  void add_dummy_hadronic_process(double reaction_cross_section);

  /**
   * Add the photonic process. Also compute the total cross section as a side
   * effect.
   */
  void add_single_process() {
    add_processes<CollisionBranch>(brems_cross_sections(),
                                   collision_processes_bremsstrahlung_,
                                   cross_section_bremsstrahlung_);
  }

  /**
   * Enum for encoding the photon process. It is uniquely determined by the
   * incoming particles. The naming scheme is :
   * Incoming_1_Incoming_2_.
   */
  enum class ReactionType {
    no_reaction,
    pi_z_pi_m,
    pi_z_pi_p,
    pi_p_pi_m,
    pi_m_pi_m,
    pi_p_pi_p,
    pi_z_pi_z
  };

  /**
   * Determine photon process from incoming particles.
   *
   * If incoming particles are not part of any implemented photonic process,
   * return no_reaction.
   *
   * \param[in] in ParticleList of incoming particles.
   * \return ReactionType enum-member
   */
  static ReactionType bremsstrahlung_reaction_type(const ParticleList &in);

  /**
   * Check if particles can undergo an implemented photon process.
   *
   * This function does not check the involved kinematics.
   *
   * \param[in] in ParticleList of incoming particles.
   * \return bool if photon reaction implemented.
   */
  static bool is_bremsstrahlung_reaction(const ParticleList &in) {
    return bremsstrahlung_reaction_type(in) != ReactionType::no_reaction;
  }

 private:
  /**
   * Holds the bremsstrahlung branch. As of now, this will always
   * hold only one branch.
   */
  CollisionBranchList collision_processes_bremsstrahlung_;

  /// Reaction process as determined from incoming particles.
  const ReactionType reac_;

  /**
   * Number of photons created for each hadronic scattering, needed for correct
   * weighting. Note that in generate_final_state() only one photon + two pions
   * are created.
   */
  const int number_of_fractional_photons_;

  /// Weight of the produced photon.
  double weight_ = 0.0;

  /// Total cross section of bremsstrahlung process.
  double cross_section_bremsstrahlung_ = 0.0;

  /// Total hadronic cross section
  const double hadronic_cross_section_;

  /// Sampled value of k (photon momentum)
  double k_;

  /// Sampled value of theta (angle of the photon)
  double theta_;

  /// Type of spin interaction to use
  const SpinInteractionType spin_interaction_type_;

  /**
   * Create interpolation objects for tabularized cross sections:
   * total cross section, differential dSigma/dk, differential dSigma/dtheta
   */
  void create_interpolations();

  /**
   * Computes the total cross section of the bremsstrahlung process.
   *
   * \returns List of photon reaction branches.
   */
  CollisionBranchList brems_cross_sections();

  /**
   * Computes the differential cross sections dSigma/dk and dSigma/dtheta of the
   * bremsstrahlung process.
   *
   * \returns Pair containing dSigma/dk as a first argument and dSigma/dtheta
   *          as a second argument
   */
  std::pair<double, double> brems_diff_cross_sections();
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_BREMSSTRAHLUNGACTION_H_
