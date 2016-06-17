/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SCATTERACTION_H_
#define SRC_INCLUDE_SCATTERACTION_H_

#include "action.h"


namespace Smash {


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
   */
  ScatterAction(const ParticleData &in_part1, const ParticleData &in_part2,
                float time, bool isotropic = false,
                float formation_time = 1.0f);

  /** Add a new collision channel. */
  void add_collision(CollisionBranchPtr p);
  /** Add several new collision channels at once. */
  void add_collisions(CollisionBranchList pv);

  /**
   * Calculate the transverse distance of the two incoming particles in their
   * local rest frame.
   *
   * Returns the squared distance.
   *
   * \fpPrecision Why \c double?
   */
  double transverse_distance_sqr() const;

  /**
   * Generate the final-state of the scattering process.
   * Performs either elastic or inelastic scattering.
   *
   * \throws InvalidResonanceFormation
   */
  void generate_final_state() override;

  float raw_weight_value() const override;

  /** Add all possible subprocesses for this action object. */
  void add_all_processes(float elastic_parameter,
                         bool two_to_one, bool two_to_two, bool strings_switch);

  /**
   * Determine the (parametrized) total cross section for this collision. This
   * is currently only used for calculating the string excitation cross section.
   */
  virtual float total_cross_section() const {
    return 0.;
  }

  /**
   * Determine the (parametrized) elastic cross section for this collision.
   * It is zero by default, but can be overridden in the child classes.
   */
  virtual float elastic_parametrization() {
    return 0.;
  }

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
  CollisionBranchPtr elastic_cross_section(float elast_par);

  /**
   * Determine the cross section for string excitations, which is given by the
   * difference between the parametrized total cross section and all the
   * explicitly implemented channels at low energy (elastic, resonance
   * excitation, etc). This method has to be called after all other processes
   * have been added to the Action object.
   */
  virtual CollisionBranchPtr string_excitation_cross_section();

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
  virtual CollisionBranchList resonance_cross_sections();

  /**
   * Return the 2-to-1 resonance production cross section for a given resonance.
   *
   * \param[in] type_resonance Type information for the resonance to be produced.
   * \param[in] srts Total energy in the center-of-mass frame.
   * \param[in] cm_momentum_sqr Square of the center-of-mass momentum of the
   * two initial particles.
   *
   * \return The cross section for the process
   * [initial particle a] + [initial particle b] -> resonance.
   *
   * \fpPrecision Why \c double?
   */
  double two_to_one_formation(const ParticleType &type_resonance,
                              double srts, double cm_momentum_sqr);

  /** Find all inelastic 2->2 processes for this reaction. */
  virtual CollisionBranchList two_to_two_cross_sections() {
    return CollisionBranchList();
  }

  /** Determine the total energy in the center-of-mass frame,
   * i.e. sqrt of Mandelstam s.  */
  double sqrt_s() const override;
  /**
   * \ingroup exception
   * Thrown when ScatterAction is called to perform with unknown ProcessType.
   */
  class InvalidScatterAction : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };

  float cross_section() const {
    return total_cross_section_;
  }

 protected:
  /** Determine the Mandelstam s variable,
   * s = (p_a + p_b)^2 = square of CMS energy.
   *
   * \fpPrecision Why \c double?
   */
  double mandelstam_s() const;
  /** Determine the momenta of the incoming particles in the
   * center-of-mass system.
   * \fpPrecision Why \c double?
   */
  double cm_momentum() const;
  /** Determine the squared momenta of the incoming particles in the
   * center-of-mass system.
   * \fpPrecision Why \c double?
   */
  double cm_momentum_squared() const;
  /// determine the velocity of the center-of-mass frame in the lab
  ThreeVector beta_cm() const;
  /// determine the corresponding gamma factor
  double gamma_cm() const;

  /** Perform an elastic two-body scattering, i.e. just exchange momentum. */
  void elastic_scattering();

  /** Perform the string excitation and decay via Pythia
   */
  void string_excitation();
  
  /**
   * \ingroup logging
   * Writes information about this scatter action to the \p out stream.
   */
  void format_debug_output(std::ostream &out) const override;

  /** List of possible collisions  */
  CollisionBranchList collision_channels_;

  /** Total cross section */
  float total_cross_section_;

  /** Do this collision isotropically. */
  bool isotropic_ = false;

  /** Formation time parameter for string fragmentation*/
  float formation_time_ = 1.0f;

 private:
  /** Check if the scattering is elastic. */
  bool is_elastic() const;

  /** Perform a 2->1 resonance-formation process. */
  void resonance_formation();
};


}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTION_H_
