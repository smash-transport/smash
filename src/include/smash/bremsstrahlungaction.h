/*
 *
 *    Copyright (c) 2016-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_BREMSSTRAHLUNG_H_
#define SRC_INCLUDE_BREMSSTRAHLUNG_H_

#include "scatteraction.h"

#include "smash/crosssectionsbrems.h"

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
   * \return The constructed object.
   */

  BremsstrahlungAction(const ParticleList &in, const double time,
                       const int n_frac_photons,
                       const double hadronic_cross_section_input);
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
    pi_p_pi_p
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
   * weighting. Note that in generate_final_state() only one photon + hadron is
   * created.
   */
  const int number_of_fractional_photons_;

  /// Weight of the produced photon.
  double weight_ = 0.0;

  /// Total cross section of bremsstrahlung process.
  double cross_section_bremsstrahlung_ = 0.0;

  /// Total hadronic cross section
  const double hadronic_cross_section_;

  /**
   * Create interpolation objects for tabularized cross sections. There are only
   * two channels; one involving two charged pions (marked pipi) and one
   * involving a charged and a neutral pion (marked pi0pi). All other cross
   * sections can be determined from crossing symmetries.
   */
  void create_interpolations() {
    // Read in tabularized cross section data (sqrt(s) and sigma)
    std::vector<double> x_pipi_opp_charge = BREMS_PIPI_OPP_C_SQRTS;
    std::vector<double> y_pipi_opp_charge = BREMS_PIPI_OPP_C_SIG;

    std::vector<double> x_pipi_same_charge = BREMS_PIPI_SAME_C_SQRTS;
    std::vector<double> y_pipi_same_charge = BREMS_PIPI_SAME_C_SIG;

    std::vector<double> x_pi0pi = BREMS_PI0PI_SQRTS;
    std::vector<double> y_pi0pi = BREMS_PI0PI_SIG;

    // Create interpolation object containing linear interpolations
    pipi_opp_charge_interpolation = make_unique<InterpolateDataLinear<double>>(
        x_pipi_opp_charge, y_pipi_opp_charge);
    pipi_same_charge_interpolation = make_unique<InterpolateDataLinear<double>>(
        x_pipi_same_charge, y_pipi_same_charge);
    pi0pi_interpolation =
        make_unique<InterpolateDataLinear<double>>(x_pi0pi, y_pi0pi);
  }

  /**
   * Computes the total cross section of the bremsstrahlung process.
   *
   * \returns List of photon reaction branches.
   */
  CollisionBranchList brems_cross_sections();
};

}  // namespace smash

#endif  // SRC_INCLUDE_BREMSSTRAHLUNG_H_
