/*
 *
 *    Copyright (c) 2016-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SCATTERACTIONPHOTON_H_
#define SRC_INCLUDE_SCATTERACTIONPHOTON_H_

#include <utility>

#include "scatteraction.h"

namespace smash {

/**
 * \ingroup action
 * ScatterActionPhoton is a special action which takes two incoming particles
 * and performs a perturbative electromagnetic scattering.
 * The final state particles are not further propagated, only written
 * to the output.
 */

class ScatterActionPhoton : public ScatterAction {
 public:
  /**
   * Construct a ScatterActionPhoton object.
   *
   * \param[in] in ParticleList of incoming particles.
   * \param[in] time Time relative to underlying hadronic action.
   * \param[in] n_frac_photons Number of photons to produce for each hadronic
   *                            scattering. 
   * \param[in] hadronic_cross_section_input Cross-section of
   *                                          underlying hadronic cross-section.
   * \return The constructed object.
   */

  ScatterActionPhoton(const ParticleList &in, const double time,
                      const int n_frac_photons,
                      const double hadronic_cross_section_input);

  /**
   * Create the photon final state and write to output.
   *
   * \param[in] outputs List of all outputs. Does not have to be a specific
   *                     photon output, the function will take care of this.
   */
  void perform_photons(const OutputsList &outputs);

  /**
   * Generate the final-state for the photon scatter process. Generates only one
   * photon / hadron pair
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
   *
   * \return total cross-section [mb]
   */
  double hadronic_cross_section() const { return hadronic_cross_section_; }

  /**
   * Sample the mass of the outgoing hadron. Returns the pole mass if
   * particle is stable.
   *
   * \param[in] out_type TypePtr of the outgoing hadron.
   *
   * \return Mass of outgoing hadron [GeV]
   *
   */
  double sample_out_hadron_mass(const ParticleTypePtr out_type);

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
    add_processes<CollisionBranch>(photon_cross_sections(),
                                   collision_processes_photons_,
                                   cross_section_photons_);
  }

  /**
   * Enum for encoding the photon process. It is uniquely determined by the
   * incoming particles. The naming scheme is :
   * Incoming_1__Incoming_2__Outgoing_hadron. The photon is omitted in the
   * naming.
   */
  enum class ReactionType {
    no_reaction,
    pi_z_pi_p_rho_p,
    pi_z_pi_m_rho_m,
    pi_p_rho_z_pi_p,
    pi_m_rho_z_pi_m,
    pi_m_rho_p_pi_z,
    pi_p_rho_m_pi_z,
    pi_z_rho_p_pi_p,
    pi_z_rho_m_pi_m,
    pi_p_pi_m_rho_z,
    pi_z_rho_z_pi_z
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
  static ReactionType photon_reaction_type(const ParticleList &in);

  /**
   * Check if particles can undergo an implemented photon process.
   * 
   * This function does not check the involved kinematics.
   *
   * \param[in] in ParticleList of incoming particles.
   * \return bool if photon reaction implemented.
   */
  static bool is_photon_reaction(const ParticleList &in) {
    return photon_reaction_type(in) != ReactionType::no_reaction;
  }

  /**
   * Return ParticleTypePtr of hadron in the out channel, given the incoming
   * particles.
   *
   * This function is overloaded since we need the hadron type in different
   * places.
   *
   * \param [in] in ParticleList of incoming particles.
   * \return ParticeTypePtr to hadron in outgoing channel.
   */
  static ParticleTypePtr outgoing_hadron_type(const ParticleList &in);

  /**
   * Return ParticleTypePtr of hadron in the out channel,
   * given the ReactionType. 
   *
   * This function is overloaded since we need the hadron type in different
   * places.
   *
   * \param [in] reaction ReactionType, determined from incoming particles.
   * \returns ParticeTypePtr to hadron in outgoing channel.
   */
  static ParticleTypePtr outgoing_hadron_type(const ReactionType reaction);

  /**
   *  Check if CM-energy is sufficient to produce hadron in final state.
   *
   *  \param [in] s_sqrt CM-energy [GeV]
   *  \param [in] in ParticleList of incoming hadrons
   *  \returns true if particles can be produced.
   */
  static bool is_kinematically_possible(const double s_sqrt,
                                        const ParticleList &in);

 private:
  /**
   * Holds the photon branch. As of now, this will always
   * hold only one branch.
   */
  CollisionBranchList collision_processes_photons_;

  /// Photonic process as determined from incoming particles.
  const ReactionType reac_;

  /**
   * Number of photons created for each hadronic scattering, needed for correct
   * weighting. Note that in generate_final_state() only one photon + hadron is
   * created.
   */
  const int number_of_fractional_photons_;

  /// ParticleTypePtr to the type of the outgoing hadron.
  const ParticleTypePtr hadron_out_t_;

  /// Mass of outgoing hadron
  const double hadron_out_mass_;

  /**
   * Compile-time switch for setting the handling of processes which can happen
   * via different mediating particles. Relevant only for the processes
   * pi0 + rho => pi + y and pi + rho => pi0 + gamma, which both can happen
   * via exchange of (rho, a1, pi) or omega.
   * If MediatorType::SUM is set, the cross section for both processes is added.
   * If MediatorType::PION/ OMEGA is set, only the respective processes are
   * computed.
   */
  enum class MediatorType { SUM, PION, OMEGA };
  /// Value used for default exchange particle. See MediatorType.
  static constexpr MediatorType default_mediator_ = MediatorType::SUM;

  /// Weight of the produced photon.
  double weight_ = 0.0;

  /// Total cross section of photonic process.
  double cross_section_photons_ = 0.0;

  /// Total hadronic cross section
  const double hadronic_cross_section_;

  /**
   * Calculate the differential cross section of  photon process.
   * Formfactors are not included
   *
   * \param[in] t Mandelstam-t [GeV^2].
   * \param[in] m_rho Mass of the incoming or outgoing rho-particle [GeV]
   * \param[in] mediator Switch for determing which mediating particle to use
   *
   * \return Differential cross section. [mb/\f$GeV^2\f$]
   */
  double diff_cross_section(const double t, const double m_rho,
                            MediatorType mediator = default_mediator_) const;

  /**
   * Find the mass of the participating rho-particle. 
   *
   * In case of a rho in the incoming channel it is the mass of the incoming 
   * rho, in case of an rho in the outgoing channel it is the mass sampled in
   * the constructor. When an rho acts in addition as a mediator, its mass is 
   * the same as the incoming / outgoing rho. This function returns the alrady
   * sampled mass or the mass of the incoming rho, depending on the process.
   *
   * \returns mass of participating rho [GeV]
   */
  double rho_mass() const;

  /**
   * Computes the total cross section of the photon process.
   *
   * \param[in] mediator Switch for determing which mediating particle to use.
   * \returns List of photon reaction branches.
   */
  CollisionBranchList photon_cross_sections(
      MediatorType mediator = default_mediator_);

  /**
   * For processes which can happen via (pi, a1, rho) and omega exchange,
   * return the differential cross section for the (pi, a1, rho) process in
   * the first argument, for the omega process in the second. If only
   * one process exists, both values are the same.
   *
   * \param[in] t Mandelstam-t [GeV^2]
   * \param[in] m_rho Mass of the incoming or outgoing rho-particle [GeV]
   *
   * \returns diff. cross section for (pi,a1,rho) in the first argument,
   *           for omega in the second.
   */
  std::pair<double, double> diff_cross_section_single(const double t,
                                                      const double m_rho);

  /**
   * For processes which can happen via (pi, a1, rho) and omega exchange,
   * return the form factor for the (pi, a1, rho) process in
   * the first argument, for the omega process in the second. If only
   * one process exists, both values are the same.
   *
   * \param[in] E_photon Energy of the photon [GeV]
   *
   * \return Form factor for (pi,a1,rho) in the first argument,
   * for omega in the second.
   */
  std::pair<double, double> form_factor_single(const double E_photon);

  /**
   * Compute the form factor for a process with a pion as the lightest exchange
   * particle.
   *
   * See wiki for details how form factors are handled.
   *
   * \param[in] E_photon Energy of photon [GeV]
   * \returns form factor
   */
  double form_factor_pion(const double E_photon) const;

  /**
   * Compute the form factor for a process with a omega as the lightest exchange
   * particle.
   *
   * See wiki for details how form factors are handled.
   *
   * \param[in] E_photon Energy of photon [GeV]
   * \returns form factor
   */
  double form_factor_omega(const double E_photon) const;

  /**
   * Compute the differential cross section with form factors included.
   *
   * Takes care of correct handling of reactions with multiple processes by
   * reading the default_mediator_ member variable.
   *
   * \param[in] t Mandelstam-t [GeV^2]
   * \param[in] m_rho Mass of the incoming or outgoing rho-particle [GeV]
   * \param[in] E_photon of outgoing photon [GeV]
   *
   * \returns diff. cross section [mb / GeV \f$^2\f$]
   */
  double diff_cross_section_w_ff(const double t, const double m_rho,
                                 const double E_photon);
};

}  // namespace smash

#endif  // SRC_INCLUDE_SCATTERACTIONPHOTON_H_
