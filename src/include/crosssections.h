/*
 *
 *    Copyright (c) 2015-2017
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

// TODO(staudenmaier): Class documentation
class cross_sections {
 public:
  cross_sections(const ParticleList& scat_particles, const double sqrt_s);

  /**
   * Generate a list of all possible collisions between the incoming particles
   * with the given c.m. energy and the calculated cross sections.
   * \param[in] elastic_parameter If non-zero, given global elastic cross
   *                                                            section.
   * \param[in] two_to_one_switch 2->1 reactions enabled?
   * \param[in] included_2to2 Which 2->2 ractions are enabled?
   * \param[in] low_snn_cut Elastic collisions with CME below are forbidden.
   * \param[in] strings_switch Are string processes enabled?
   * \param[in] nnbar_treatment NNbar treatment through resonance, strings or
   *                                                        none
   * \param[in] string_process String process used for string fragmentaion.
   * \return List of all possible collisions.
   */
  CollisionBranchList generate_collision_list(
      double elastic_parameter, bool two_to_one_switch,
      ReactionsBitSet included_2to2, double low_snn_cut, bool strings_switch,
      NNbarTreatment nnbar_treatment, StringProcess* string_process);

  /**
   * Determine the elastic cross section for this collision. If elastic_par is
   * given (and positive), we just use a constant cross section of that size,
   * otherwise a parametrization of the elastic cross section is used
   * (if available).
   *
   * \param[in] elast_par Elastic cross section parameter from the input file.
   *
   * \return A ProcessBranch object containing the cross section and
   * final-state IDs.
   */
  CollisionBranchPtr elastic(double elast_par);

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

  /** Find all inelastic 2->2 processes for the given scattering.
   * This function calls the different, more specific functions for
   * the different scatterings.
   * \param[in] included_2to2 Which 2->2 ractions are enabled?
   * \return List of all possibe inelastic 2->2 processes.
   */
  CollisionBranchList two_to_two(ReactionsBitSet included_2to2);

  /**
   * Determine the cross section for string excitations, which is given by the
   * difference between the parametrized total cross section and all the
   * explicitly implemented channels at low energy (elastic, resonance
   * excitation, etc).
   * \param[in] string_process String process used for string fragmentaion.
   *
   * \return List of subprocesses (single-diffractive,
   * double-diffractive and non-diffractive) with their cross sections
   *
   * This method has to be called after all other processes
   * have been determined.
   * \todo Same assumption made by NNbar_annihilation. Resolve.
   */
  CollisionBranchList string_excitation(StringProcess* string_process);

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
   * Determine the cross section for NNbar annihilation, which is given by
   * detailed balance from the reverse reaction. See
   * NNbar_annihilation_cross_section
   */
  CollisionBranchList NNbar_creation();

 private:
  /// Choose between parametrization for elastic cross sections.
  double elastic_parametrization();

  /** TODO(staudenmaier): Revise documentation from here on. Also .cc file.
   * Determine the (parametrized) elastic cross section for a
   * nucleon-nucleon collision.
   */
  double nn_el();

  /**
   * Determine the elastic cross section for a nucleon-pion collision.
   * It is given by a parametrization of experimental data.
   */
  double npi_el();

  /**
   * Determine the elastic cross section for a nucleon-kaon collision.
   * It is given by a parametrization of experimental data.
   */
  double nk_el();

  /**
   * Find all inelastic 2->2 processes for Baryon-Baryon Scattering
   * except the more specific Nucleon-Nucleon Scattering.
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
   * \return List of resonance production processes possible in the collision
   * of the two nucleons. Each element in the list contains the type(s) of the
   * final state particle(s) and the cross section for that particular process.
   */
  CollisionBranchList nn_xx(ReactionsBitSet included_2to2);

  /** Find all inelastic 2->2 processes for Nucelon-Kaon Scattering. */
  CollisionBranchList nk_xx(ReactionsBitSet included_2to2);

  /** Find all inelastic 2->2 processes for Delta-Kaon Scattering. */
  CollisionBranchList deltak_xx(ReactionsBitSet included_2to2);

  /** Find all inelastic 2->2 processes for Hyperon-Pion Scattering. */
  CollisionBranchList ypi_xx(ReactionsBitSet included_2to2);

  /** Determine the parametrized total cross section at high energies
   * for the given collision, which is non-zero for Baryon-Baryon and
   * Nucleon-Pion scatterings currently.
   */
  double high_energy() const;

  /**
   * Determine the (parametrized) hard non-diffractive string cross section
   * for this collision.
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
   * \param[in] sqrts	sqrt(Mandelstam-s), i.e. collision CMS energy.
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
   */
  template <class IntegrationMethod>
  CollisionBranchList find_nn_xsection_from_type(
      const ParticleTypePtrList& type_res_1,
      const ParticleTypePtrList& type_res_2,
      const IntegrationMethod integrator);

  /** Return, if the scattering between the incoming particles are scattering
   * via string fragmentaion or not.
   * The string fragmentation is implemented in the same way in GiBUU (Physics
   * Reports 512(2012), 1-124, pg. 33). If the center of mass energy is low, two
   * particles scatter through the resonance channels. If high, the out going
   * particles are generated by string fragmentation. If in between, the out
   * going particles are generated either through the resonance channels or
   * string fragmentation by chance. In detail, the low energy regoin is from
   * the threshold to (mix_scatter_type_energy - mix_scatter_type_window_width),
   * while the high energy region is from (mix_scatter_type_energy +
   * mix_scatter_type_window_width) to infinity. In between, the probability for
   * string fragmentation increases linearly from 0 to 1 as the c.m. energy.
   */
  bool decide_string(bool strings_switch, const bool both_are_nucleons) const;

  /** Determine the momenta of the incoming particles in the
   * center-of-mass system.
   */
  double cm_momentum() const {
    const double m1 = incoming_particles_[0].effective_mass();
    const double m2 = incoming_particles_[1].effective_mass();
    return pCM(sqrt_s_, m1, m2);
  }

  /** List with data of scattering particles.  */
  ParticleList incoming_particles_;

  /** total energy in the center-of-mass frame. */
  double sqrt_s_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_CROSSSECTIONS_H_
