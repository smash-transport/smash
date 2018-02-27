/*
 *
 *    Copyright (c) 2015-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SCATTERACTION_H_
#define SRC_INCLUDE_SCATTERACTION_H_

#include <memory>
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
   * \param[in] string_formation_time the time a string takes to form
   */
  ScatterAction(const ParticleData& in_part1, const ParticleData& in_part2,
                double time, bool isotropic = false,
                double string_formation_time = 1.0);

  /** Add a new collision channel. */
  void add_collision(CollisionBranchPtr p);
  /** Add several new collision channels at once. */
  void add_collisions(CollisionBranchList pv);

  /**
   * Calculate the transverse distance of the two incoming particles in their
   * local rest frame.
   *
   * Returns the squared distance.
   */
  double transverse_distance_sqr() const;

  /**
   * Generate the final-state of the scattering process.
   * Performs either elastic or inelastic scattering.
   *
   * \throws InvalidResonanceFormation
   */
  void generate_final_state() override;

  double raw_weight_value() const override;

  double partial_weight() const override;

  /**
   * Sample final-state angles in a 2->2 collision (possibly anisotropic).
   *
   * \throws InvalidResonanceFormation
   */
  void sample_angles(std::pair<double, double> masses) override;

  /** Add all possible scattering subprocesses for this action object. */
  void add_all_scatterings(double elastic_parameter, bool two_to_one,
                                   bool two_to_two, double low_snn_cut,
                                   bool strings_switch,
                                   NNbarTreatment nnbar_treatment);

  /// Returns list of possible collision channels
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

  void set_string_interface(StringProcess* str_proc) {
    string_process_ = str_proc;
  }

  virtual double cross_section() const { return total_cross_section_; }

 protected:
  /** Determine the Mandelstam s variable,
   * s = (p_a + p_b)^2 = square of CMS energy.
   */
  double mandelstam_s() const;
  /** Determine the momenta of the incoming particles in the
   * center-of-mass system.
   */
  double cm_momentum() const;
  /** Determine the squared momenta of the incoming particles in the
   * center-of-mass system.
   */
  double cm_momentum_squared() const;
  /// determine the velocity of the center-of-mass frame in the lab
  ThreeVector beta_cm() const;
  /// determine the corresponding gamma factor
  double gamma_cm() const;

  /** Perform an elastic two-body scattering, i.e. just exchange momentum. */
  void elastic_scattering();

  /** Perform an inelastic two-body scattering, i.e. new particles are formed*/
  void inelastic_scattering();

  /** Todo(ryu): document better - it is not really UrQMD-based, isn't it?
   *  Perform the UrQMD-based string excitation and decay */
  void string_excitation_soft();

  /** Perform the string excitation and decay via Pythia. */
  void string_excitation_pythia();

  /**
   * \ingroup logging
   * Writes information about this scatter action to the \p out stream.
   */
  void format_debug_output(std::ostream& out) const override;

  /** List of possible collisions  */
  CollisionBranchList collision_channels_;

  /** Total cross section */
  double total_cross_section_;

  /** Partial cross-section to the chosen outgoing channel */
  double partial_cross_section_;

  /** Do this collision isotropically. */
  bool isotropic_ = false;

  /** Formation time parameter for string fragmentation*/
  double string_formation_time_ = 1.0;

 private:
  /** Check if the scattering is elastic. */
  bool is_elastic() const;

  /** Perform a 2->1 resonance-formation process. */
  void resonance_formation();

  /** Pointer to interface class for strings */
  StringProcess* string_process_ = nullptr;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SCATTERACTION_H_
