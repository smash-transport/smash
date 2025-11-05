/*
 *
 *    Copyright (c) 2014-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_SCATTERACTION_H_
#define SRC_INCLUDE_SMASH_SCATTERACTION_H_

#include <memory>
#include <set>
#include <string>
#include <utility>

#include "action.h"
#include "isoparticletype.h"
#include "scatteractionsfinderparameters.h"
#include "stringprocess.h"

namespace smash {

/**
 * \ingroup action
 * ScatterAction is a special action which takes two incoming particles
 * and performs a scattering, producing one or more final-state particles.
 */
class ScatterAction : public Action {
 public:
  /**
   * Construct a ScatterAction object.
   *
   * \param[in] in_part1 first scattering partner
   * \param[in] in_part2 second scattering partner
   * \param[in] time Time at which the action is supposed to take place
   * \param[in] isotropic if true, do the collision isotropically
   * \param[in] string_formation_time Time string fragments take to form
   * \param[in] box_length Passing box length to determine
   *            coordinate of the collision, in case it happened through
   *            the wall in a box. If negative, then there is no wrapping.
   * \param[in] is_total_parametrized Whether the total cross section used for
   * collision finding is parametrized
   * \param[in] spin_interaction_type Which type of spin interaction to use
   */
  ScatterAction(const ParticleData& in_part1, const ParticleData& in_part2,
                double time, bool isotropic = false,
                double string_formation_time = 1.0, double box_length = -1.0,
                bool is_total_parametrized = false,
                const SpinInteractionType spin_interaction_type =
                    SpinInteractionType::Off);
  /**
   * Add a new collision channel.
   *
   * \param[in] p Channel to be added.
   */
  void add_collision(CollisionBranchPtr p);

  /**
   * Add several new collision channels at once.
   *
   * \param[in] pv list of channels to be added.
   */
  void add_collisions(CollisionBranchList pv);

  /**
   * Calculate the transverse distance of the two incoming particles in their
   * local rest frame.
   *
   * According to UrQMD criterion, \iref{Bass:1998ca} eq. (3.27):
   *  - position of particle a: \f$\mathbf{x}_a\f$
   *  - position of particle b: \f$\mathbf{x}_b\f$
   *  - momentum of particle a: \f$\mathbf{p}_a\f$
   *  - momentum of particle b: \f$\mathbf{p}_b\f$
   *
   * \f[
   * d^2_\mathrm{coll}
   * = (\mathbf{x}_a - \mathbf{x}_b)^2 -
   *   \frac{\bigl[(\mathbf{x}_a - \mathbf{x}_b) \cdot
   *               (\mathbf{p}_a - \mathbf{p}_b)\bigr]^2 }
   *        {(\mathbf{p}_a - \mathbf{p}_b)^2}
   * \f]
   *
   * \return  squared distance \f$d^2_\mathrm{coll}\f$.
   */
  double transverse_distance_sqr() const;

  /**
   * Calculate the transverse distance of the two incoming particles in their
   * local rest frame written in a covariant form. Equivalent to the UrQMD
   * transverse distance. See \iref{Hirano:2012yy} (5.6)-(5.11).
   *
   * \return squared distance  \f$d^2_\mathrm{coll}\f$.
   */
  double cov_transverse_distance_sqr() const;
  /**
   * Determine the Mandelstam s variable,
   *
   * \f[s = (p_a + p_b)^2\f]
   * Equal to the square of CMS energy.
   *
   * \return Mandelstam s
   */
  double mandelstam_s() const;

  /**
   * Get the relative velocity of the two incoming
   * particles. For a defintion see e.g. \iref{Seifert:2017oyb}, eq. (5)
   *
   * \return relative velocity.
   */
  double relative_velocity() const;

  /**
   * Generate the final-state of the scattering process.
   * Performs either elastic or inelastic scattering.
   *
   * \throw InvalidScatterAction

   */
  void generate_final_state() override;

  /**
   * Get the total cross section of scattering particles.
   *
   * \return total cross section.
   */
  double get_total_weight() const override;

  /**
   * Get the partial cross section of the chosen channel.
   *
   * \return partial cross section.
   */
  double get_partial_weight() const override;

  /**
   * Sample final-state angles in a 2->2 collision (possibly anisotropic).
   */
  void sample_angles(std::pair<double, double> masses,
                     double kinetic_energy_cm) override;

  /**
   * Add all possible scattering subprocesses for this action object. This can
   * only be called once per ScatterAction instance.
   *
   * \param[in] finder_parameters parameters for collision finding.
   */
  void add_all_scatterings(
      const ScatterActionsFinderParameters& finder_parameters);

  /**
   * Given the incoming particles, assigns the correct parametrization of the
   * total cross section.
   *
   * \param[in] finder_parameters Parameters for collision finding.
   */
  void set_parametrized_total_cross_section(
      const ScatterActionsFinderParameters& finder_parameters);

  /**
   * Get list of possible collision channels.
   *
   * \return list of possible collision channels.
   */
  const CollisionBranchList& collision_channels() {
    return collision_channels_;
  }

  /**
   * \ingroup exception
   * Thrown when ScatterAction is called to perform with unknown ProcessType.
   */
  class InvalidScatterAction : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };

  /**
   * Set the StringProcess object to be used.
   *
   * The StringProcess object is used to handle string excitation and to
   * generate final state particles.
   *
   * \param[in] str_proc String process object to be used.
   */
  void set_string_interface(StringProcess* str_proc) {
    string_process_ = str_proc;
  }

  /**
   * Get the total cross section of the scattering particles, either from a
   * parametrization, or from the sum of partials.
   *
   * \return total cross section.
   */
  virtual double cross_section() const {
    if (is_total_parametrized_) {
      return *parametrized_total_cross_section_;
    }
    return sum_of_partial_cross_sections_;
  }

 protected:
  /**
   * Get the momentum of the center of mass of the incoming particles
   * in the calculation frame.
   *
   * \return center of mass momentum.
   */
  double cm_momentum() const;
  /**
   * Get the squared momentum of the center of mass of the incoming
   * particles in the calculation frame.
   *
   * \return center of mass momentum squared.
   */
  double cm_momentum_squared() const;

  /**
   * Get the velocity of the center of mass of the scattering/incoming particles
   * in the calculation frame.
   *
   * Note: Do not use this function to boost the outgoing
   * particles. Use total_momentum_of_outgoing_particles(), which corrects for
   * the effect of potentials on intial and final state.
   *
   * \return boost velocity between center of mass and calculation frame.
   */
  ThreeVector beta_cm() const;
  /**
   * Get the gamma factor corresponding to a boost to the center of mass frame
   * of the colliding particles.
   *
   * \return gamma factor.
   */
  double gamma_cm() const;

  /// Perform an elastic two-body scattering, i.e. just exchange momentum.
  void elastic_scattering();

  /// Perform an inelastic two-body scattering, i.e. new particles are formed
  void inelastic_scattering();

  /// Perform an inelastic two-to-many-body scattering (more than 2)
  void two_to_many_scattering();

  /// Creates the final states for string-processes after they are performed
  void create_string_final_state();
  /**
   * Todo(ryu): document better - it is not really UrQMD-based, isn't it?
   * Perform the UrQMD-based string excitation and decay
   */
  void string_excitation();

  /**
   * Perform spin interaction in binary interactions. At the moment, we include
   * a spin-flip in the y component of elastic scatterings if enabled.
   */
  void spin_interaction();

  /**
   * \ingroup logging
   * Writes information about this scatter action to the \p out stream.
   */
  void format_debug_output(std::ostream& out) const override;

  /// List of possible collisions
  CollisionBranchList collision_channels_;

  /// Current sum of partial hadronic cross sections
  double sum_of_partial_cross_sections_;

  /// Partial cross-section to the chosen outgoing channel
  double partial_cross_section_;

  /// Do this collision isotropically?
  bool isotropic_ = false;

  /// Time fragments take to be fully formed in hard string excitation.
  double string_formation_time_ = 1.0;

 private:
  /**
   * Check if the scattering is elastic.
   *
   * \return whether the scattering is elastic.
   */
  bool is_elastic() const;

  /**
   * Perform a 2->1 resonance-formation process.
   * \throw InvalidResonanceFormation
   */
  void resonance_formation();

  /**
   * Loop over the possible branches and rescales their weight according to the
   * desired total cross section. In case the current sum of partials is close
   * to 0, a warning is issued as this would not happen in an usual run, and an
   * elastic process is added to match the total.
   */
  void rescale_outgoing_branches();

  /**
   * Try to find a pseudo-resonance that can be created from the incoming
   * particles using a given method.
   *
   * \param[in] method used to select the pseudo-resonance among possible
   * candidates. \see_key{key_CT_pseudoresonance_}
   * \param[in] transition parameters for the string transition region, which
   * are also used to determine when a pseudo-resonance can be created.
   * \return the appropriate pseudo-resonance, if there is any, or an invalid
   * pointer otherwise.
   */
  ParticleTypePtr try_find_pseudoresonance(
      const PseudoResonance method,
      const StringTransitionParameters& transition) const;

  /// Pointer to interface class for strings
  StringProcess* string_process_ = nullptr;

  /// Whether the total cross section is parametrized
  bool is_total_parametrized_ = false;

  /// If cross section is parametrized, store the value
  std::optional<double> parametrized_total_cross_section_ = std::nullopt;

  /// What kind of spin interaction to use
  SpinInteractionType spin_interaction_type_ = SpinInteractionType::Off;

  /// Lock for calling add_all_scatterings only once
  bool were_processes_added_ = false;

  /// Warn about zero cross section only once per particle type pair
  static inline std::set<std::set<ParticleTypePtr>>
      warned_no_rescaling_available{};
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_SCATTERACTION_H_
