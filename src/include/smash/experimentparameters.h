/*
 *    Copyright (c) 2014-2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_SMASH_EXPERIMENTPARAMETERS_H_
#define SRC_INCLUDE_SMASH_EXPERIMENTPARAMETERS_H_

#include <memory>
#include <optional>
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

  /// mode of calculating gradients for density calculation
  DerivativesMode derivatives_mode;
  /// mode of calculating rest frame density gradients (on or off)
  RestFrameDensityDerivativesMode rho_derivatives_mode;
  /// mode of calculating field derivatives
  FieldDerivativesMode field_derivatives_mode;

  /// mode of smearing for density calculation
  const SmearingMode smearing_mode;

  /// Width of gaussian Wigner density of particles
  double gaussian_sigma;

  /// Distance at which gaussian is cut, i.e. set to zero, IN SIGMA (not fm)
  double gauss_cutoff_in_sigma;

  /// Weight applied to the center cell in the discrete smearing
  double discrete_weight;

  /// Range of lattice nodes, in units of lattice spacing, that the
  /// triangular smearing uses
  double triangular_range;

  /// Employed collision criterion
  const CollisionCriterion coll_crit;

  /// This indicates whether string fragmentation is switched on.
  bool two_to_one;

  /// This indicates which two to two reactions are switched off.
  const ReactionsBitSet included_2to2;

  /// This indicates which multi-particle reactions are switched on.
  const MultiParticleReactionsBitSet included_multi;

  /// This indicates whether string fragmentation is switched on.
  bool strings_switch;

  /**
   * Multiplicative factor to be applied to resonance lifetimes; in the case of
   * thermal multiplicities this should also be applied to initial
   * multiplicities of resonances, so that one does not artificially introduce
   * a non-zero pion chemical potential.
   */
  double res_lifetime_factor;

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

  /**
   * Fixed minimal grid cell length (in fm). Only used and useful in case of the
   * stochastic criterion, where the grid cell size is calculation parameter.
   * Only particles within the cells are checked for collisions. The cell
   * length is not scaled by the number of test-particles.
   *
   * Default of 2.5 fm produces (without test-particles) equivalent-sized grid
   * as the maximum cross section default of 200 fm.
   *
   * In case of the geometric criteria, the cell size is determined by the
   * maximum cross section.
   */
  double fixed_min_cell_length;  // fm

  /**
   * Global factor which all cross sections are scaled with.
   *
   * Using it will break agreement with experimental data for cross sections
   * that are constrained with data.
   */
  double scale_xs;

  /** In thermodynamics outputs, it decides whether to use only participants
   * (true) or also spectators (false, default value).
   */
  bool only_participants;

  /// Do non-strong decays at the end
  bool do_non_strong_decays;

  /**
   * Whether to decay initial state particles.
   */
  bool decay_initial_particles;

  /**
   * Whether to include spin interactions.
   */
  SpinInteractionType spin_interaction_type;

  /** Bool for the default usage of the monash tune in the collider modus.
   * The used type is std::optional since its value might not be known at
   * creation time. E.g. in Experiment this flag is set after the instance is
   * created.
   */
  std::optional<bool> use_monash_tune_default;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_EXPERIMENTPARAMETERS_H_
