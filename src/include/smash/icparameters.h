/*
 *    Copyright (c) 2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_SMASH_ICPARAMETERS_H_
#define SRC_INCLUDE_SMASH_ICPARAMETERS_H_

#include <optional>

namespace smash {

/*
 * This struct is separated from the rest of ColliderModus because its
 * implementation needs to be known by Experiment, so including this file won't
 * include the full collidermodus.h
 */
/**
 * At the moment there are two ways to specify input for initial conditions in
 * the configuration, one of which is deprecated and will be removed in a next
 * release. For the moment, these variables are of type
 * `std::optional<double>` to *allow* for the key duplication consistently.
 * When more types of IC are implemented in the future, this will allow
 * setting only the appropriate parameters.
 */
struct InitialConditionParameters {
  /// Hypersurface proper time in IC
  std::optional<double> proper_time = std::nullopt;
  /// Lower bound for proper time in IC
  std::optional<double> lower_bound = std::nullopt;
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
  /// Number of interpolating cells for dynamic IC
  std::optional<int> num_fluid_cells = std::nullopt;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_ICPARAMETERS_H_
