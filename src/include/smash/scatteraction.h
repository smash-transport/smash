/*
 *
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SCATTERACTION_H_
#define SRC_INCLUDE_SCATTERACTION_H_

#include <memory>
#include <set>
#include <string>
#include <utility>

#include "action.h"
#include "cxx14compat.h"
#include "isoparticletype.h"
#include "kinematics.h"
#include "processstring.h"

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
   */
  ScatterAction(const ParticleData& in_part1, const ParticleData& in_part2,
                double time, bool isotropic = false,
                double string_formation_time = 1.0);

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
   * According to UrQMD criterion
   * position of particle a: x_a \n
   * position of particle b: x_b \n
   * momentum of particle a: p_a \n
   * momentum of particle b: p_b \n
   * \f[d^2_\mathrm{coll} = (\vec{x_a} - \vec{x_b})^2 - \frac{((\vec{x_a} -
   * \vec{x_b}) \cdot (\vec{p_a} - \vec{p_b}))^2 } {(\vec{p_a} -
   * \vec{p_b})^2}\f]
   *
   * \return  squared distance \f$d^2_\mathrm{coll}\f$.
   */
  double transverse_distance_sqr() const;

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
  void sample_angles(std::pair<double, double> masses) override;

  /**
   * Add all possible scattering subprocesses for this action object.
   *
   * \param[in] elastic_parameter If non-zero, given global
   *            elastic cross section.
   * \param[in] two_to_one  2->1 reactions enabled?
   * \param[in] included_2to2 Which 2->2 reactions are enabled?
   * \param[in] low_snn_cut Elastic collisions with CME below are forbidden.
   * \param[in] strings_switch Are string processes enabled?
   * \param[in] use_AQM use elastic cross sections via AQM?
   * \param[in] strings_with_probability Are string processes triggered
   *            according to a probability?
   * \param[in] nnbar_treatment NNbar treatment through resonance, strings or
   *                                                        none
   */
  void add_all_scatterings(double elastic_parameter, bool two_to_one,
                           ReactionsBitSet included_2to2, double low_snn_cut,
                           bool strings_switch, bool use_AQM,
                           bool strings_with_probability,
                           NNbarTreatment nnbar_treatment);

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
   * Get the total cross section of the scattering particles.
   *
   * \return total cross section.
   */
  virtual double cross_section() const { return total_cross_section_; }

 protected:
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
   * Get the velocity of the center of mass of the scattering particles
   * in the calculation frame.
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

  /**
   * Todo(ryu): document better - it is not really UrQMD-based, isn't it?
   * Perform the UrQMD-based string excitation and decay
   */
  void string_excitation();

  /**
   * \ingroup logging
   * Writes information about this scatter action to the \p out stream.
   */
  void format_debug_output(std::ostream& out) const override;

  /// List of possible collisions
  CollisionBranchList collision_channels_;

  /// Total hadronic cross section
  double total_cross_section_;

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

  /// Pointer to interface class for strings
  StringProcess* string_process_ = nullptr;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SCATTERACTION_H_
