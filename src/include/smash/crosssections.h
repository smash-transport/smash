/*
 *
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_CROSSSECTIONS_H_
#define SRC_INCLUDE_CROSSSECTIONS_H_

#include "forwarddeclarations.h"
#include "isoparticletype.h"
#include "particles.h"
#include "processstring.h"

namespace smash {

/**
 * The cross section class assembels everything that is needed to
 * calculate the cross section and returns a list of all possible reactions
 * for the incoming particles at the given energy with the calculated cross
 * sections.
 */
class CrossSections {
 public:
  /**
   * Construct CrossSections instance.
   *
   * \param[in] incoming_particles Particles that are reacting.
   * \param[in] sqrt_s Center-of-mass energy of the reaction.
   */
  CrossSections(const ParticleList& incoming_particles, const double sqrt_s);

  /**
   * Generate a list of all possible collisions between the incoming particles
   * with the given c.m. energy and the calculated cross sections.
   * The string processes are not added at this step if it's not triggerd
   * according to the probability. It will then be added in
   * add_all_scatterings in scatteraction.cc
   *
   * \param[in] elastic_parameter Value of the constant global elastic cross
   *            section, if it is non-zero.
   *            The parametrized elastic cross section is used otherwise.
   * \param[in] two_to_one_switch 2->1 reactions enabled?
   * \param[in] included_2to2 Which 2->2 ractions are enabled?
   * \param[in] low_snn_cut Elastic collisions with CME below are forbidden.
   * \param[in] strings_switch Are string processes enabled?
   * \param[in] use_AQM Is the Additive Quark Model enabled?
   * \param[in] strings_with_probability Are string processes triggered
   *            according to a probability?
   * \param[in] nnbar_treatment NNbar treatment through resonance, strings or
   *                                                        none
   * \param[in] string_process a pointer to the StringProcess object,
   *            which is used for string excitation and fragmentation.
   * \return List of all possible collisions.
   */
  CollisionBranchList generate_collision_list(
      double elastic_parameter, bool two_to_one_switch,
      ReactionsBitSet included_2to2, double low_snn_cut, bool strings_switch,
      bool use_AQM, bool strings_with_probability,
      NNbarTreatment nnbar_treatment, StringProcess* string_process);

  /**
   * Determine the elastic cross section for this collision. If elastic_par is
   * given (and positive), we just use a constant cross section of that size,
   * otherwise a parametrization of the elastic cross section is used
   * (if available).
   *
   * \param[in] elast_par Elastic cross section parameter from the input file.
   * \param[in] use_AQM Whether to extend elastic cross-sections with AQM.
   *
   * \return A ProcessBranch object containing the cross section and
   * final-state IDs.
   */
  CollisionBranchPtr elastic(double elast_par, bool use_AQM);

  /**
   * Find all resonances that can be produced in a 2->1 collision of the two
   * input particles and the production cross sections of these resonances.
   *
   * Given the data and type information of two colliding particles,
   * create a list of possible resonance production processes
   * and their cross sections.
   *
   * \return A list of processes with resonance in the final state.
   * Each element in the list contains the type of the final-state particle
   * and the cross section for that particular process.
   */
  CollisionBranchList two_to_one();

  /**
   * Return the 2-to-1 resonance production cross section for a given resonance.
   *
   * \param[in] type_resonance Type information for the resonance to be
   * produced.
   * \param[in] cm_momentum_sqr Square of the center-of-mass momentum of the
   * two initial particles.
   *
   * \return The cross section for the process
   * [initial particle a] + [initial particle b] -> resonance.
   */
  double formation(const ParticleType& type_resonance, double cm_momentum_sqr);

  /**
   * Find all inelastic 2->2 processes for the given scattering.
   *
   * This function calls the different, more specific functions for
   * the different scatterings.
   *
   * \param[in] included_2to2 Which 2->2 reactions are enabled?
   * \return List of all possibe inelastic 2->2 processes.
   */
  CollisionBranchList two_to_two(ReactionsBitSet included_2to2);

  /**
   * Determine the cross section for string excitations, which is given by the
   * difference between the parametrized total cross section and all the
   * explicitly implemented channels at low energy (elastic, resonance
   * excitation, etc).
   *
   * \param[in] total_string_xs Total cross section for the string process [mb].
   * \param[in] string_process a pointer to the StringProcess object,
   *            which is used for string excitation and fragmentation.
   * \param[in] use_AQM whether to extend string cross-sections with AQM
   * \return List of subprocesses (single-diffractive,
   *        double-diffractive and non-diffractive) with their cross sections.
   *
   * \throw std::runtime_error
   *        if string_process is a null pointer.
   *
   * This method has to be called after all other processes
   * have been determined.
   * \todo Same assumption made by NNbar_annihilation. Resolve.
   */
  CollisionBranchList string_excitation(double total_string_xs,
                                        StringProcess* string_process,
                                        bool use_AQM);

  /**
   * Determine the cross section for NNbar annihilation, which is given by the
   * difference between the parametrized total cross section and all the
   * explicitly implemented channels at low energy (in this case only elastic).
   * \param[in] current_xs Sum of all cross sections of already determined
   *                                                     processes
   * \return Collision Branch with NNbar annihilation process and its cross
   *   section
   *
   * This method has to be called after all other processes
   * have been determined.
   * \todo Same assumption made by string_excitation. Resolve.
   */
  CollisionBranchPtr NNbar_annihilation(const double current_xs);

  /**
   * Determine the cross section for NNbar creation, which is given by
   * detailed balance from the reverse reaction. See
   * NNbar_annihilation_cross_section
   * \return Collision Branch with NNbar creation process and its cross
   * section
   */
  CollisionBranchList NNbar_creation();

  /**
   * Determine the parametrized total cross section at high energies
   * for the given collision, which is non-zero for Baryon-Baryon and
   * Nucleon-Pion scatterings currently.
   *
   * This is rescaled by AQM factors.
   */
  double high_energy() const;

  /**
   * Return, if the scattering between the incoming particles are scattering
   * via string fragmentation or not.
   *
   * If use_transition_probability is true:
   * The string fragmentation is implemented in the same way in GiBUU (Physics
   * Reports 512(2012), 1-124, pg. 33). If the center of mass energy is low, two
   * particles scatter through the resonance channels. If high, the outgoing
   * particles are generated by string fragmentation. If in between, the out-
   * going particles are generated either through the resonance channels or
   * string fragmentation by chance. In detail, the low energy region is from
   * the threshold to (mix_scatter_type_energy - mix_scatter_type_window_width),
   * while the high energy region is from (mix_scatter_type_energy +
   * mix_scatter_type_window_width) to infinity. In between, the probability for
   * string fragmentation increases smoothly from 0 to 1 as the c.m. energy.
   *
   * If use_transition_probability is false:
   * The string fragmentation is implemented similarly to what is in UrQMD
   * (\iref{Bass:1998ca}). If sqrts is lower than some cutoff value, there are
   * no strings. If higher, strings are allowed, with the cross-section being
   * the difference between some parametrized total cross-section and the sum
   * of all other channels, if this parametrization is larger than the sum of
   * the channels. If not, strings are not allowed (this cross-section check
   * is performed directly after the function is called, for technical reasons).
   *
   * Both of these methods are initially implemented for NN and Npi cross-
   * sections, and extended using the AQM to all BB, BM and MM interactions.
   *
   * Baryon-antibaryon annihilation also uses this function to decide
   * whether to produce strings or not.
   * Since there are no other contributions for this process,
   * there are no cutoffs or gradual increase in the probability of this process
   * happening or not, it just requires the proper combination of incoming
   * particles and config parameters.
   *
   * \param[in] strings_switch Is string fragmentation enabled?
   * \param[in] use_transition_probability which algorithm to use for string
   *                         treatment (see Switch_on_String_with_Probability)
   * \param[in] use_AQM whether AQM is activated
   * \param[in] treat_nnbar_with_strings use strings for nnbar treatment?
   *
   * \return Is the scattering between the incoming particles done via string·
   * fragmentation or not?
   */
  bool decide_string(bool strings_switch, bool use_transition_probability,
                     bool use_AQM, bool treat_nnbar_with_strings) const;

  /**
   * \return if the species of the two incoming particles are allowed to
   * interact via string fragmentation. Currently, only nucleon-nucleon
   * and nucleon-pion can interact via string.
   */
  bool included_in_string() const;

 private:
  /**
   * Choose the appropriate parametrizations for given incoming particles and
   * return the (parametrized) elastic cross section.
   *
   * \param[in] use_AQM whether AQM is activated
   * \return Elastic cross section
   */
  double elastic_parametrization(bool use_AQM);

  /**
   * Determine the (parametrized) elastic cross section for a
   * nucleon-nucleon (NN) collision.
   * \return Elastic cross section for NN
   *
   * \throw std::runtime_error
   *        if positive cross section cannot be specified.
   */
  double nn_el();

  /**
   * Determine the elastic cross section for a nucleon-pion (Npi) collision.
   * It is given by a parametrization of experimental data.
   * \return Elastic cross section for Npi
   *
   * \throw std::runtime_error
   *        if incoming particles are not nucleon+pion.
   * \throw std::runtime_error
   *        if positive cross section cannot be specified.
   */
  double npi_el();

  /**
   * Determine the elastic cross section for a nucleon-kaon (NK) collision.
   * It is given by a parametrization of experimental data.
   * \return Elastic cross section for NK
   *
   * \throw std::runtime_error
   *        if incoming particles are not nucleon+kaon.
   * \throw std::runtime_error
   *        if positive cross section cannot be specified.
   */
  double nk_el();

  /**
   * Find all inelastic 2->2 processes for Baryon-Baryon (BB) Scattering
   * except the more specific Nucleon-Nucleon Scattering.
   * \param[in] included_2to2 Which 2->2 reactions are enabled?
   * \return List of all possible BB reactions with their cross sections
   */
  CollisionBranchList bb_xx_except_nn(ReactionsBitSet included_2to2);

  /**
   * Find all inelastic 2->2 processes for Nucelon-Nucelon Scattering.
   * Calculate cross sections for resonance production from
   * nucleon-nucleon collisions (i.e. N N -> N R, N N -> Delta R).
   *
   * Checks are processed in the following order:
   * 1. Charge conservation
   * 2. Isospin factors (Clebsch-Gordan)
   * 3. Enough energy for all decay channels to be available for the resonance
   *
   * \param[in] included_2to2 Which 2->2 reactions are enabled?
   *
   * \return List of resonance production processes possible in the collision
   * of the two nucleons. Each element in the list contains the type(s) of the
   * final state particle(s) and the cross section for that particular process.
   */
  CollisionBranchList nn_xx(ReactionsBitSet included_2to2);

  /**
   * Find all inelastic 2->2 background processes for Nucleon-Kaon (NK)
   * Scattering.
   * \param[in] included_2to2 Which 2->2 reactions are enabled?
   * \return List of all possible NK reactions with their cross sections
   */
  CollisionBranchList nk_xx(ReactionsBitSet included_2to2);

  /**
   * Find all inelastic 2->2 processes for Delta-Kaon (DeltaK) Scattering.
   * \param[in] included_2to2 Which 2->2 reactions are enabled?
   * \return List of all possible DeltaK reactions with their cross sections
   * */
  CollisionBranchList deltak_xx(ReactionsBitSet included_2to2);

  /**
   * Find all inelastic 2->2 processes for Hyperon-Pion (Ypi) Scattering.
   * \param[in] included_2to2 Which 2->2 reactions are enabled?
   * \return List of all possible Ypi reactions with their cross sections
   */
  CollisionBranchList ypi_xx(ReactionsBitSet included_2to2);

  /**
   * Find all inelastic 2->2 processes involving Pion and (anti-) Deuteron
   * (dpi), specifically dπ→ NN, d̅π→ N̅N̅; πd→ πd' (mockup for πd→ πnp), πd̅→ πd̅'
   * and reverse.
   * \param[in] included_2to2 Which 2->2 reactions are enabled?
   * \return List of all possible dpi reactions with their cross sections
   */
  CollisionBranchList dpi_xx(ReactionsBitSet included_2to2);

  /**
   * Find all inelastic 2->2 processes involving Nucleon and (anti-) Deuteron
   * (dN), specifically Nd → Nd', N̅d →  N̅d', N̅d̅→ N̅d̅', Nd̅→ Nd̅' and reverse (e.g.
   * Nd'→ Nd).
   * \param[in] included_2to2 Which 2->2 reactions are enabled?
   * \return List of all possible dN reactions with their cross sections
   */
  CollisionBranchList dn_xx(ReactionsBitSet included_2to2);

  /**
   * Determine the (parametrized) hard non-diffractive string cross section
   * for this collision.
   *
   * \return Parametrized cross section without AQM scaling.
   */
  double string_hard_cross_section() const;

  /**
   * Calculate cross sections for resonance absorption
   * (i.e. NR->NN and ΔR->NN).
   *
   * \param[in] is_anti_particles Whether the colliding particles are
   * antiparticles
   *
   * \return List of possible resonance absorption processes. Each element of
   * the list contains the types of the final-state particles and the cross
   * section for that particular process.
   */
  CollisionBranchList bar_bar_to_nuc_nuc(const bool is_anti_particles);

  /**
   * Scattering matrix amplitude squared (divided by 16π) for resonance
   * production processes like NN → NR and NN → ΔR, where R is a baryon
   * resonance (Δ, N*, Δ*). Includes no spin or isospin factors.
   *
   * \param[in] sqrts sqrt(Mandelstam-s), i.e. collision CMS energy.
   * \param[in] type_a Type information for the first final-state particle.
   * \param[in] type_b Type information for the second final-state particle.
   * \param[in] twoI Twice the total isospin of the involved state.
   *
   * \return Matrix amplitude squared \f$ |\mathcal{M}(\sqrt{s})|^2/16\pi \f$.
   */
  static double nn_to_resonance_matrix_element(double sqrts,
                                               const ParticleType& type_a,
                                               const ParticleType& type_b,
                                               const int twoI);

  /**
   * Utility function to avoid code replication in nn_xx().
   * \param[in] type_res_1 List of possible first final resonance types
   * \param[in] type_res_2 List of possible second final resonance types
   * \param[in] integrator Used to integrate over the kinematically allowed
   * mass range of the Breit-Wigner distribution
   * \return List of all possible NN reactions with their cross sections
   * with different final states
   */
  template <class IntegrationMethod>
  CollisionBranchList find_nn_xsection_from_type(
      const ParticleTypePtrList& type_res_1,
      const ParticleTypePtrList& type_res_2,
      const IntegrationMethod integrator);

  /**
   * Determine the momenta of the incoming particles in the
   * center-of-mass system.
   * \return Center-of-mass momentum
   */
  double cm_momentum() const {
    const double m1 = incoming_particles_[0].effective_mass();
    const double m2 = incoming_particles_[1].effective_mass();
    return pCM(sqrt_s_, m1, m2);
  }

  /// List with data of scattering particles.
  ParticleList incoming_particles_;

  /// Total energy in the center-of-mass frame.
  double sqrt_s_;

  /// Whether incoming particles are a baryon-antibaryon pair
  const bool is_BBbar_pair_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_CROSSSECTIONS_H_
