/*
 *    Copyright (c) 2013-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_SMASH_EXPERIMENTPARAMETERS_H_
#define SRC_INCLUDE_SMASH_EXPERIMENTPARAMETERS_H_

#include <memory>
#include <set>

#include "clock.h"

namespace smash {

/**
 * Helper structure for Experiment.
 *
 * Experiment has one member of this struct. In essence the members of this
 * struct are members of Experiment, but combined in one structure for easier
 * function argument passing.
 */
struct ExperimentParameters {
  /// System clock (for simulation time keeping in the computational frame)
  std::unique_ptr<Clock> labclock;

  /// Output clock to keep track of the next output time
  std::unique_ptr<Clock> outputclock;

  /// Number of parallel ensembles
  int n_ensembles;

  /// Number of test-particles
  int testparticles;

  // mode of calculating gradients
  const DerivativesMode derivatives_mode;

  /// Width of gaussian Wigner density of particles
  double gaussian_sigma;

  /// Distance at which gaussian is cut, i.e. set to zero, IN SIGMA (not fm)
  double gauss_cutoff_in_sigma;

  /// Employed collision criterion
  const CollisionCriterion coll_crit;

  /// This indicates whether two to one reactions are switched on.
  bool two_to_one;

  /// This indicates which two to two reactions are switched off.
  const ReactionsBitSet included_2to2;

  /// This indicates which multi-particle reactions are switched on.
  const MultiParticleReactionsBitSet included_multi;

  /// This indicates whether string fragmentation is switched on.
  bool strings_switch;

  /// Whether to use the AQM or not
  bool use_AQM;

  /**
   * Multiplicative factor to be applied to resonance lifetimes; in the case of
   * thermal multiplicities this should also be applied to initial
   * multiplicities of resonances, so that one does not artificially introduce
   * a non-zero pion chemical potential.
   */
  double res_lifetime_factor;

  /**
   * This indicates whether the string fragmentation is swiched on with
   * a probability smoothly increasing with energy. If it's set equal to
   * false, the cross section of the string fragmentation is counted by
   * taking the difference between the parametrized total cross section
   * and the sum of the non-string cross sections.
   */
  bool strings_with_probability;
  /**
   * This indicates how NN̅ annihilation should be treated; options are to
   * neglect it, make it conserve detailed balance using NN̅ → h₁(1170)ρ
   * (which goes to 5 pions on average) or use strings.
   */
  NNbarTreatment nnbar_treatment;

  /**
   * Elastic collisions between the nucleons with the square root s
   * below low_snn_cut are excluded.
   */
  double low_snn_cut;

  /**
   * This indicates whether the mean field potentials affect the scattering
   * or decaying processes by shifting the threshold energies.
   */
  bool potential_affect_threshold;

  /**
   * Length of the box in fm in case of box modus, otherwise -1
   */
  double box_length;

  /**
   * The maximal cross section (in mb) for which it is guaranteed that all
   * collisions with this cross section will be found.
   *
   * This means that all particle pairs, where the transverse distance
   * is smaller or equal to \f$ \sqrt{\sigma_{max}/\pi} \f$,
   * will be checked for collions.
   *
   * The maximal cross section is also scaled with the cross section
   * scaling factor.
   */
  double maximum_cross_section;  // mb

  /// Allow or forbid the first collisions within the same nucleus
  bool allow_collisions_within_nucleus;
  /**
   * Global factor which all cross sections are scaled with.
   *
   * Using it will break agreement with experimental data for cross sections
   * that are constrained with data.
   */
  double scale_xs;

  /**
   * Additional constant contribution (in mb) to the elastic cross sections.
   *
   * Using it will break agreement with experimental data for elastic cross
   * sections that are constrained with data.
   */
  double additional_el_xs;  // mb
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_EXPERIMENTPARAMETERS_H_
