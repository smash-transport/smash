/*
 *    Copyright (c) 2024-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_SMASH_ICPARAMETERS_H_
#define SRC_INCLUDE_SMASH_ICPARAMETERS_H_

#include <optional>

namespace smash {
/**
 * The variables in this POD struct are of type `std::optional<double>` so that
 * only the relevant parameters are set for the different types of initial
 * conditions. Only \c type is the exception, as it is a required key for the IC
 * runs.
 */
struct InitialConditionParameters {
  /// Type of initialization
  FluidizationType type;
  /// Which processes can have outgoing particles transformed into fluid in
  /// dynamic IC
  std::optional<FluidizableProcessesBitSet> fluidizable_processes =
      std::nullopt;
  /// Hypersurface proper time in IC
  std::optional<double> proper_time = std::nullopt;
  /// Lower bound for proper time in IC
  std::optional<double> lower_bound = std::nullopt;
  /// Scaling factor for proper time in IC
  std::optional<double> proper_time_scaling = std::nullopt;
  /// Rapidity cut on hypersurface IC
  std::optional<double> rapidity_cut = std::nullopt;
  /// Transverse momentum cut on hypersurface IC
  std::optional<double> pT_cut = std::nullopt;
  /// Minimum energy density for dynamic IC
  std::optional<double> energy_density_threshold = std::nullopt;
  /// Minimum time (in lab frame) in fm for dynamic IC
  std::optional<double> min_time = std::nullopt;
  /// Maximum time (in lab frame) in fm for dynamic IC
  std::optional<double> max_time = std::nullopt;
  /// Number of interpolating cells in each direction for dynamic IC
  std::optional<int> num_fluid_cells = std::nullopt;
  /**
   * Fraction of formation time to pass before particles can fluidize in
   * dynamic IC
   */
  std::optional<double> formation_time_fraction = std::nullopt;
  /// Smearing kernel at 0 for dynamic IC
  std::optional<double> smearing_kernel_at_0 = std::nullopt;
  /// Whether the first elastic interaction of an initial nucleon is fluidizable
  std::optional<bool> delay_initial_elastic = std::nullopt;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_ICPARAMETERS_H_
