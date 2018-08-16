/*
 *    Copyright (c) 2013-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_EXPERIMENTPARAMETERS_H_
#define SRC_INCLUDE_EXPERIMENTPARAMETERS_H_

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
  Clock labclock;

  /// Output clock to keep track of the next output time
  Clock outputclock;

  /// Number of test particle
  int testparticles;

  /// Width of gaussian Wigner density of particles
  double gaussian_sigma;

  /// Distance at which gaussian is cut, i.e. set to zero, IN SIGMA (not fm)
  double gauss_cutoff_in_sigma;

  /// This indicates whether two to one reactions are switched on.
  bool two_to_one;

  /// This indicates which two to two reactions are switched off.
  const ReactionsBitSet included_2to2;

  /// This indicates whether string fragmentation is switched on.
  bool strings_switch;

  /// Whether to use the AQM or not
  bool use_AQM;

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

  /// This indicates whether photons are switched on.
  bool photons_switch;

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
};

}  // namespace smash

#endif  // SRC_INCLUDE_EXPERIMENTPARAMETERS_H_
