/*
 *
 *    Copyright (c) 2013-2014,2018-2024,2026
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_CROSSSECTIONS_H_
#define SRC_INCLUDE_SMASH_CROSSSECTIONS_H_

#include <memory>
#include <optional>
#include <utility>

#include "smash/forwarddeclarations.h"
#include "smash/isoparticletype.h"
#include "smash/particles.h"
#include "smash/potential_globals.h"
#include "smash/scatteractionsfinderparameters.h"
#include "smash/stringprocess.h"

namespace smash {

/**
 * The CrossSections class assembles everything that is needed to calculate
 * cross sections and gathers a list of all possible reactions for the incoming
 * particles at the given energy with the calculated cross sections.
 */
class CrossSections {
 public:
  /**
   * Construct CrossSections instance.
   *
   * \param[in] incoming_particles Particles that are interacting.
   * \param[in] sqrt_s Center-of-mass energy of the reaction.
   * \param[in] potentials Potentials at the interacting point. They are
   *            used to calculate the corrections on the thresholds.
   */
  CrossSections(const ParticleList& incoming_particles, double sqrt_s,
                const std::pair<FourVector, FourVector> potentials);

  /**
   * Generate a list of all possible collisions between the incoming particles
   * with the given c.m. energy and the calculated cross sections.
   * The string processes are not added at this step if it is not triggerd
   * according to the probability. It will then be added in
   * add_all_scatterings in scatteraction.cc
   *
   * \param[in] finder_parameters parameters for collision finding.
   * \param[in] string_process a pointer to the StringProcess object,
   * which is used for string excitation and fragmentation.
   *
   * \return List of all possible collisions.
   */
  CollisionBranchList generate_collision_list(
      const ScatterActionsFinderParameters& finder_parameters,
      StringProcess* string_process) const;

  /**
   * Select the parametrization for the total cross section, given the types of
   * incoming particles.
   *
   * \param[in] finder_parameters Parameters for collision finding, containing
   * cut for low energy NN interactions.
   *
   * \return The appropriate total cross section value.
   */
  double parametrized_total(
      const ScatterActionsFinderParameters& finder_parameters) const;

  /**
   * Helper function:
   * Sum all cross sections of the given process list.
   */
  static double sum_xs_of(const CollisionBranchList& list) {
    double xs_sum = 0.0;
    for (auto& proc : list) {
      xs_sum += proc->weight();
    }
    return xs_sum;
  }

  /**
   * Determine the elastic cross section for this collision. If elastic_par is
   * given (and positive), we just use a constant cross section of that size,
   * otherwise a parametrization of the elastic cross section is used
   * (if available). Optional a constant additional elastic cross section is
   * added
   *
   * \param[in] finder_parameters parameters for collision finding, including
   * cross section modifications from config file.
   *
   * \note The additional constant elastic cross section contribution is added
   * after the scaling of the cross section.
   *
   * \return A ProcessBranch object containing the cross section and
   * final-state IDs.
   */
  CollisionBranchPtr elastic(
      const ScatterActionsFinderParameters& finder_parameters) const;

  /**
   * Find all resonances that can be produced in a 2->1 collision of the two
   * input particles and the production cross sections of these resonances.
   *
   * Given the data and type information of two colliding particles,
   * create a list of possible resonance production processes
   * and their cross sections.
   *
   * If <tt>\ref key_CT_charm_rescattering_ "Charm_Rescattering_Method"</tt> is
   * not set to `"Resonances"` then interactions of charmed hadrons will not be
   * considered in two to one processes.
   *
   * \param[in] charm_rescattering Type of charm rescattering
   *
   * \return A list of processes with resonance in the final state.
   * Each element in the list contains the type of the final-state particle
   * and the cross section for that particular process.
   */
  CollisionBranchList two_to_one(CharmRescattering charm_rescattering) const;

  /**
   * Calculates the 2-to-1 resonance production cross section for a given
   * resonance using the Breit-Wigner distribution as probability amplitude.
   * See eq. (176) in \iref{Buss:2011mx}.
   *
   * \param[in] type_resonance Type information for the resonance to be
   * produced.
   * \param[in] cm_momentum_sqr Square of the center-of-mass momentum of the
   * two initial particles.
   *
   * \return The cross section for the process
   * [initial particle a] + [initial particle b] -> resonance.
   */
  double formation(const ParticleType& type_resonance,
                   double cm_momentum_sqr) const;

  /**
   * Find all 2->2 processes which are suppressed at high energies when
   * strings are turned on with probabilites, but important for the
   * production of rare species such as strange particles.
   *
   * This function should call the different, more specific functions for
   * the different scatterings. But so far, only Nucleon-Pion to Hyperon-
   * Kaon scattering is implemented.
   *
   * \return List of all possible rare 2->2 processes.
   */
  CollisionBranchList rare_two_to_two() const;

  /**
   * Find all inelastic 2->2 processes for the given scattering.
   *
   * This function calls the different, more specific functions for
   * the different scatterings.
   *
   * \param[in] included_2to2 Which 2->2 reactions are enabled?
   * \param[in] KN_offset Offset to the minimum energy for string production in
   * KN scatterings
   * \param[in] charm_rescattering Type of charm rescattering
   *
   * \return List of all possible inelastic 2->2 processes.
   */
  CollisionBranchList two_to_two(const ReactionsBitSet& included_2to2,
                                 double KN_offset,
                                 CharmRescattering charm_rescattering) const;

  /**
   * Find all 2->3 processes for the given scattering.
   *
   * This function calls the different, more specific functions for
   * the different scatterings.
   *
   * \return List of all possible 2->3 processes.
   */
  CollisionBranchList two_to_three() const;

  /**
   * Find all 2->4 processes for the given scattering.
   *
   * This function calls the different, more specific functions for
   * the different scatterings.
   *
   * \return List of all possible 2->4 processes.
   */
  CollisionBranchList two_to_four() const;

  /**
   * Determine the cross section for string excitations, which is given by the
   * difference between the parametrized total cross section and all the
   * explicitly implemented channels at low energy (elastic, resonance
   * excitation, etc).
   *
   * \param[in] total_string_xs Total cross section for the string process [mb]
   * \param[in] string_process a pointer to the StringProcess object,
   *            which is used for string excitation and fragmentation
   * \param[in] finder_parameters parameters for collision finding and cross
   * sections
   * \return List of subprocesses (single-diffractive, double-diffractive and
   * non-diffractive) with their cross sections
   *
   * \throw std::runtime_error if string_process is a null pointer
   *
   * This method has to be called after all other processes
   * have been determined.
   * \todo Same assumption made by NNbar_annihilation. Resolve.
   */
  CollisionBranchList string_excitation(
      double total_string_xs, StringProcess* string_process,
      const ScatterActionsFinderParameters& finder_parameters) const;

  /**
   * Determine the cross section for NNbar annihilation, which is given by the
   * difference between the parametrized total cross section and all the
   * explicitly implemented channels at low energy (in this case only elastic).
   * \param[in] current_xs Sum of all cross sections of already determined
   *                                                     processes
   * \param[in] scale_xs Factor by which all (partial) cross sections are scaled
   * \return Collision Branch with NNbar annihilation process and its cross
   *   section
   *
   * This method has to be called after all other processes
   * have been determined.
   * \todo Same assumption made by string_excitation. Resolve.
   */
  CollisionBranchPtr NNbar_annihilation(double current_xs,
                                        double scale_xs) const;

  /**
   * Determine the cross section for NNbar creation, which is given by
   * detailed balance from the reverse reaction. See
   * NNbar_annihilation.
   * \return Collision Branch with NNbar creation process and its cross
   * section
   */
  CollisionBranchList NNbar_creation() const;

  /**
   * Create collision branch for NNbar annihilation going directly into 5 pions.
   * The cross section is given by the parametrized ppbar cross section, which
   * is also used for the reverse 5-to-2 process.
   *
   * \param[in] scale_xs Factor by which all (partial) cross sections are scaled
   * \return Collision Branch with NNbar annihilation process
   *
   */
  CollisionBranchPtr NNbar_to_5pi(double scale_xs) const;

  /**
   * Determine 2->3 cross section for the scattering of the given particle
   * types.
   *
   * That the function only depends on the types of particles (plus sqrt(s)) and
   * not on the specific particles, is an assumption needed in order to treat
   * the 3->2 back-reaction with the stochastic criterion, where this function
   * also needs to be called for 3-to-2 collision probability with only types
   * and sqrt(s) known at this point. Therefore the function is also made
   * static.
   *
   * \param[in] type_in1 first scatterning particle type
   * \param[in] type_in2 second scatterning particle type
   * \param[in] sqrts center-of-mass energy of scattering
   * \return cross section for 2->3 process
   */
  static double two_to_three_xs(const ParticleType& type_in1,
                                const ParticleType& type_in2, double sqrts);

  /**
   * Determine 2->4 cross section for the scattering of the given particle
   * types.
   *
   * Same assumptions as for 2->3 cross section, see respective documentation.
   *
   * \param[in] type_in1 first scatterning particle type
   * \param[in] type_in2 second scatterning particle type
   * \param[in] sqrts center-of-mass energy of scattering
   * \return cross section for 2->4 process
   */
  static double two_to_four_xs(const ParticleType& type_in1,
                               const ParticleType& type_in2, double sqrts);

  /**
   * Determine the parametrized total cross section at high energies
   * for the given collision, which is non-zero for Baryon-Baryon and
   * Nucleon-Pion scatterings currently.
   *
   * This is rescaled by AQM factors.
   *
   * \param[in] finder_parameters parameters for collision finding and cross
   *                              sections.
   */
  double high_energy(
      const ScatterActionsFinderParameters& finder_parameters) const;

  /**
   * \return the probability whether the scattering between the incoming
   * particles is via string fragmentation or not.
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
   * Baryon-antibaryon annihilation also uses this function to decide whether to
   * produce strings or not.
   * Since there are no other contributions for this process, there are no
   * cutoffs or gradual increase in the probability of this process happening or
   * not, it just requires the proper combination of incoming particles and
   * config parameters.
   *
   * \param[in] finder_parameters parameters for collision finding and cross
   * sections.
   */
  double string_probability(
      const ScatterActionsFinderParameters& finder_parameters) const;
  /**
   * Computes a smooth transition probability as a function of sqrt(s).
   *
   * The probability is:
   *   - 0 for sqrt(s) < region_lower
   *   - 1 for sqrt(s) > region_upper
   *   - smoothly varying between 0 and 1 inside the transition region
   *     according to a sinusoidal profile
   *
   * This probability can be used to determine the relative contribution of
   * two processes or regimes (e.g. soft ↔ hard string excitation).
   *
   * \param[in] region_lower Lower bound of the transition region in sqrt(s)
   *                         [GeV]
   * \param[in] region_upper Upper bound of the transition region in sqrt(s)
   *                         [GeV]
   *
   * \return Transition probability in the interval [0,1]
   */
  double transition_probability_at_sqrts(double region_lower,
                                         double region_upper) const;

 private:
  /**
   * Choose the appropriate parametrizations for given incoming particles and
   * return the (parametrized) elastic cross section.
   *
   * \param[in] finder_parameters parameters for collision finding and cross
   * sections, containing whether to extend string cross-sections with AQM and
   * the offset to the minimum energy for string production in \f$\pi\pi \f$
   * scatterings
   *
   * \return Elastic cross section
   *
   * \throw std::runtime_error if elastic cross section is negative.
   */
  double elastic_parametrization(
      const ScatterActionsFinderParameters& finder_parameters) const;

  /**
   * Determine the (parametrized) elastic cross section for a
   * nucleon-nucleon (NN) collision.
   * \return Elastic cross section for NN
   *
   * \throw std::runtime_error if positive cross section cannot be specified.
   */
  double nn_el() const;

  /**
   * Determine the elastic cross section for a nucleon-pion (Npi) collision.
   * It is given by a parametrization of experimental data.
   * \return Elastic cross section for Npi
   *
   * \throw std::runtime_error if incoming particles are not nucleon+pion.
   * \throw std::runtime_error if positive cross section cannot be specified.
   */
  double npi_el() const;

  /**
   * Determine the elastic cross section for a nucleon-kaon (NK) collision.
   * It is given by a parametrization of experimental data.
   * \return Elastic cross section for NK
   *
   * \throw std::runtime_error if incoming particles are not nucleon+kaon.
   * \throw std::runtime_error if positive cross section cannot be specified.
   */
  double nk_el() const;

  /**
   * Determine the elastic cross section for a D meson-pion (Dpi) or a D*-pion
   * (D*pi) collision. If the center-of-mass energy for the collision is below
   * the lower bound of the energy range of the underlying cross section data,
   * the return value is zero. If it is above the upper bound, the return value
   * will be `std::nullopt`.
   *
   * \return Elastic cross section for Dpi or D*pi.
   *
   * \throw std::runtime_error if incoming particles are not Dpi or D*pi.
   * \throw std::runtime_error if cross section is negative.
   */
  std::optional<double> Dpi_and_Dstarpi_elastic() const;

  /**
   * Determine the inelastic cross section for a D meson-pion (Dpi) or a D*-pion
   * (D*pi) collision.
   *
   * \return Inlastic cross section for Dpi or D*pi.
   *
   * \throw std::runtime_error if incoming particles are not Dpi or D*pi.
   * \throw std::runtime_error if cross section is negative.
   */
  double Dpi_and_Dstarpi_inelastic() const;

  /**
   * Determine the elastic cross section for a D meson-eta (Deta) or a D*-eta
   * (D*eta) collision. If the center-of-mass energy for the collision is below
   * the lower bound of the energy range of the underlying cross section data,
   * the return value is zero. If it is above the upper bound, the return value
   * will be `std::nullopt`.
   *
   * \return Elastic cross section for Deta or D*eta.
   *
   * \throw std::runtime_error if incoming particles are not Deta or D*eta.
   * \throw std::runtime_error if cross section is negative.
   */
  std::optional<double> Deta_and_Dstareta_elastic() const;

  /**
   * Determine the elastic cross section for a D meson-kaon (DK) or a D*-kaon
   * (D*K) collision. If the center-of-mass energy for the collision is below
   * the lower bound of the energy range of the underlying cross section data,
   * the return value is zero. If it is above the upper bound, the return value
   * will be `std::nullopt`.
   *
   * \return Elastic cross section for DK or D*K.
   *
   * \throw std::runtime_error if incoming particles are not DK or D*K.
   * \throw std::runtime_error if cross section is negative.
   */
  std::optional<double> DK_and_DstarK_elastic() const;

  /**
   * Determine the inelastic cross section for a D meson-kaon (DK) or a D*-kaon
   * (D*K) collision.
   *
   * \return Inlastic cross section for DK or D*K.
   *
   * \throw std::runtime_error if incoming particles are not DK or D*K.
   * \throw std::runtime_error if cross section is negative.
   */
  double DK_and_DstarK_inelastic() const;

  /**
   * Determine the elastic cross section for a D meson-nucleon (DN) collision,
   * If the center-of-mass energy for the collision is below the lower bound of
   * the energy range of the underlying cross section data, the return value is
   * zero. If it is above the upper bound, the return value will be
   * `std::nullopt`.
   *
   * \return Elastic cross section for DN.
   *
   * \throw std::runtime_error if incoming particles are not DN.
   * \throw std::runtime_error if cross section is negative.
   */
  std::optional<double> DN_elastic() const;

  /**
   * Find all processes for Nucleon-Pion to Hyperon-Kaon Scattering.
   * These scatterings are suppressed at high energies when strings are
   * turned on with probabilities, so they need to be added back manually.
   *
   * \return List of all possible Npi -> YK reactions with their cross sections
   */
  CollisionBranchList npi_yk() const;

  /**
   * Find all inelastic 2->2 processes for Baryon-Baryon (BB) Scattering
   * except the more specific Nucleon-Nucleon Scattering.
   * \param[in] included_2to2 Which 2->2 reactions are enabled?
   * \return List of all possible BB reactions with their cross sections
   */
  CollisionBranchList bb_xx_except_nn(
      const ReactionsBitSet& included_2to2) const;

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
  CollisionBranchList nn_xx(const ReactionsBitSet& included_2to2) const;

  /**
   * Find all inelastic 2->2 background processes for Nucleon-Kaon (NK)
   * Scattering.
   *
   * \param[in] included_2to2 Which 2->2 reactions are enabled?
   * \param[in] KN_offset Offset to the minimum energy for string production in
   * KN scatterings
   *
   * \return List of all possible NK reactions with their cross sections
   */
  CollisionBranchList nk_xx(const ReactionsBitSet& included_2to2,
                            double KN_offset) const;

  /**
   * Find all inelastic 2->2 processes for D meson-pion (Dpi) and D*-pion (D*pi)
   * scattering.
   *
   * \param[in] included_2to2 Which 2->2 reactions are enabled?
   * \param[in] charm_rescattering Type of charm rescattering
   *
   * \return List of all possible Dpi or D*pi reactions with their cross
   * sections
   */
  CollisionBranchList Dpi_and_Dstarpi_xx(
      const ReactionsBitSet& included_2to2,
      CharmRescattering charm_rescattering) const;

  /**
   * Find all inelastic 2->2 processes for D meson-kaon (DK) and D*-kaon (D*K)
   * scattering.
   *
   * \param[in] included_2to2 Which 2->2 reactions are enabled?
   * \param[in] charm_rescattering Type of charm rescattering
   *
   * \return List of all possible DK or D*K reactions with their cross sections
   */
  CollisionBranchList DK_and_DstarK_xx(
      const ReactionsBitSet& included_2to2,
      CharmRescattering charm_rescattering) const;

  /**
   * Find all inelastic 2->2 processes for Delta-Kaon (DeltaK) Scattering.
   * \param[in] included_2to2 Which 2->2 reactions are enabled?
   * \return List of all possible DeltaK reactions with their cross sections
   */
  CollisionBranchList deltak_xx(const ReactionsBitSet& included_2to2) const;

  /**
   * Find all inelastic 2->2 processes for Hyperon-Pion (Ypi) Scattering.
   * \param[in] included_2to2 Which 2->2 reactions are enabled?
   * \return List of all possible Ypi reactions with their cross sections
   */
  CollisionBranchList ypi_xx(const ReactionsBitSet& included_2to2) const;

  /**
   * Find all inelastic 2->2 processes involving Pion and (anti-) Deuteron
   * (dpi), specifically dπ→ NN, d̅π→ N̅N̅; πd→ πd' (mockup for πd→ πnp), πd̅→ πd̅'
   * and reverse.
   * \param[in] included_2to2 Which 2->2 reactions are enabled?
   * \return List of all possible dpi reactions with their cross sections
   */
  CollisionBranchList dpi_xx(const ReactionsBitSet& included_2to2) const;

  /**
   * Find all inelastic 2->2 processes involving Nucleon and (anti-) Deuteron
   * (dN), specifically Nd → Nd', N̅d →  N̅d', N̅d̅→ N̅d̅', Nd̅→ Nd̅' and reverse (e.g.
   * Nd'→ Nd).
   * \param[in] included_2to2 Which 2->2 reactions are enabled?
   * \return List of all possible dN reactions with their cross sections
   */
  CollisionBranchList dn_xx(const ReactionsBitSet& included_2to2) const;

  /**
   * Parametrized cross section for πd→ πd' (mockup for πd→ πnp), πd̅→ πd̅' and
   * reverse, see \iref{Oliinychenko:2018ugs} for details.
   * \param[in] sqrts square-root of mandelstam s
   * \param[in] cm_mom center of mass momentum of incoming particles
   * \param[in] produced_nucleus type of outgoing deuteron or d-prime
   * \param[in] type_pi type of scattering pion
   * \return cross section for given scattering
   */
  static double xs_dpi_dprimepi(double sqrts, double cm_mom,
                                ParticleTypePtr produced_nucleus,
                                const ParticleType& type_pi);

  /**
   * Parametrized cross section for Nd → Nd', N̅d →  N̅d', N̅d̅→ N̅d̅', Nd̅→ Nd̅' and
   * reverse (e.g. Nd'→ Nd), see \iref{Oliinychenko:2018ugs} for details.
   * \param[in] sqrts square-root of mandelstam s
   * \param[in] cm_mom center of mass momentum of incoming particles
   * \param[in] produced_nucleus type of outgoing deuteron or d-prime
   * \param[in] type_nucleus type of scattering (incoming) deuteron or d-prime
   * \param[in] type_N type of scattering nucleon
   * \return cross section for given scattering
   */
  static double xs_dn_dprimen(double sqrts, double cm_mom,
                              ParticleTypePtr produced_nucleus,
                              const ParticleType& type_nucleus,
                              const ParticleType& type_N);

  /**
   * Determine the (parametrized) hard non-diffractive string cross section
   * for this collision.
   *
   * \return Parametrized cross section (without AQM scaling).
   */
  double string_hard_cross_section() const;

  /**
   * Calculate cross sections for 2 → 2 resonance absorption (i.e. NR → NN and
   * ΔR → NN). See eqs. (B.6), (B.9) and (181) in \iref{Buss:2011mx}.
   *
   * \param[in] is_anti_particles Whether the colliding particles are
   * antiparticles
   *
   * \return List of possible resonance absorption processes. Each element of
   * the list contains the types of the final-state particles and the cross
   * section for that particular process.
   */
  CollisionBranchList bar_bar_to_nuc_nuc(bool is_anti_particles) const;

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
                                               int twoI);

  /**
   * Utility function to avoid code replication in nn_xx().
   *
   * \param[in] type_res_1 List of possible first final resonance types
   * \param[in] type_res_2 List of possible second final resonance types
   * \param[in] integrator Used to integrate over the kinematically allowed
   *                       mass range of the Breit-Wigner distribution
   *
   * \return List of all possible NN reactions with their cross sections
   * with different final states
   */
  template <class IntegrationMethod>
  CollisionBranchList find_nn_xsection_from_type(
      const ParticleTypePtrList& type_res_1,
      const ParticleTypePtrList& type_res_2,
      const IntegrationMethod integrator) const;

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
  const ParticleList incoming_particles_;

  /// Total energy in the center-of-mass frame.
  const double sqrt_s_;

  /**
   * Potentials at the interacting point.
   * They are used to calculate the corrections on the threshold energies.
   */
  const std::pair<FourVector, FourVector> potentials_;

  /**
   * Whether incoming particles are a pair of a baryon and an antibaryon
   * (could be different baryon types)
   */
  const bool is_BBbar_pair_;

  /// Whether incoming particles are a nulecon-antinucleon pair (same isospin)
  const bool is_NNbar_pair_;

  /**
   * Helper function:
   * Add a 2-to-2 channel to a collision branch list given a cross section.
   *
   * The cross section is only calculated if there is enough energy
   * for the process. If the cross section is small, the branch is not added.
   */
  template <typename F>
  void add_channel(CollisionBranchList& process_list, F&& get_xsection,
                   double sqrts, const ParticleType& type_a,
                   const ParticleType& type_b) const {
    const double sqrt_s_min =
        type_a.min_mass_spectral() + type_b.min_mass_spectral();
    /* Determine wether the process is below the threshold. */
    double scale_B = 0.0;
    double scale_I3 = 0.0;
    bool is_below_threshold;
    FourVector incoming_momentum = FourVector();
    if (pot_pointer != nullptr) {
      for (const auto& p : incoming_particles_) {
        incoming_momentum += p.momentum();
        scale_B += pot_pointer->force_scale(p.type()).first;
        scale_I3 +=
            pot_pointer->force_scale(p.type()).second * p.type().isospin3_rel();
      }
      scale_B -= pot_pointer->force_scale(type_a).first;
      scale_I3 -=
          pot_pointer->force_scale(type_a).second * type_a.isospin3_rel();
      scale_B -= pot_pointer->force_scale(type_b).first;
      scale_I3 -=
          pot_pointer->force_scale(type_b).second * type_b.isospin3_rel();
      is_below_threshold = (incoming_momentum + potentials_.first * scale_B +
                            potentials_.second * scale_I3)
                               .abs() <= sqrt_s_min;
    } else {
      is_below_threshold = (sqrts <= sqrt_s_min);
    }
    if (is_below_threshold) {
      return;
    }
    const auto xsection = get_xsection();
    if (xsection > really_small) {
      process_list.push_back(std::make_unique<CollisionBranch>(
          type_a, type_b, xsection, ProcessType::TwoToTwo));
    }
  }
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_CROSSSECTIONS_H_
