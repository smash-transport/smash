/*
 *    Copyright (c) 2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_SMASH_SCATTERACTIONSFINDERPARAMETERS_H_
#define SRC_INCLUDE_SMASH_SCATTERACTIONSFINDERPARAMETERS_H_

#include <memory>
#include <set>

namespace smash {

/**
 * Helper structure for ScatterActionsFinder.
 *
 * ScatterActionsFinder has one member of this struct, which just collects
 * general parameters, for easier function argument passing.
 */
struct ScatterActionsFinderParameters {
  /// Elastic cross section parameter (in mb).
  const double elastic_parameter;
  /**
   * Elastic collsions between two nucleons with sqrt_s below low_snn_cut_ are
   * excluded.
   */
  const double low_snn_cut;
  /// Factor by which all (partial) cross sections are scaled
  const double scale_xs;
  /**
   * Additional constant contribution (in mb) to the elastic cross sections.
   *
   * Using it will break agreement with experimental data for elastic cross
   * sections that are constrained with data.
   */
  const double additional_el_xs;  // mb
  /// \see input_collision_term_
  const double maximum_cross_section;
  /// Specifies which collision criterion is used
  const CollisionCriterion coll_crit;
  /// Switch for NNbar reactions
  const NNbarTreatment nnbar_treatment;
  /// List of included 2<->2 reactions
  const ReactionsBitSet included_2to2;
  /// List of included multi-particle reactions
  const MultiParticleReactionsBitSet included_multi;
  /// Number of test particles.
  const int testparticles;
  /// Enables resonance production
  const bool two_to_one;
  /// If particles within the same nucleus are allowed to collide for their
  /// first time
  const bool allow_collisions_within_nucleus;
  /// Indicates whether string fragmentation is switched on
  const bool strings_switch;
  /// Switch to control whether to use AQM or not
  const bool use_AQM;
  /**
   * This indicates whether the string fragmentation is swiched on with
   * a probability smoothly increasing with energy. If it's set equal to
   * false, the cross section of the string fragmentation is counted by
   * taking the difference between the parametrized total cross section
   * and the sum of the non-string cross sections.
   */
  const bool strings_with_probability;
  /**
   * Switch to turn off throwing an exception for collision probabilities larger
   * than 1. In larger production runs it is ok, if the probability rarely slips
   * over 1.
   */
  const bool only_warn_for_high_prob;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_SCATTERACTIONSFINDERPARAMETERS_H_
