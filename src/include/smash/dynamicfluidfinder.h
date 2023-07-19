/*
 *
 *    Copyright (c) 2023
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_DYNAMICFLUIDFINDER_H_
#define SRC_INCLUDE_SMASH_DYNAMICFLUIDFINDER_H_ 

#include <map>
#include <vector>

#include "actionfinderfactory.h"
#include "input_keys.h"

namespace smash {

/**
 * \ingroup action
 * Finder for dynamic fluidizations.
 * Loops through all particles and checks if they reach the energy density
 * threshold. This happens at the end of every time step for all hadrons that
 * originate in a decay or string fragmentation. For the latter process,
 * fluidization itself happens only after the formation time of the particle,
 * but it will happen if it passes the threshold at any point in the
 * propagation.
 */
class DynamicFluidizationFinder : public ActionFinderInterface {
 public:
  /**
   * Construct finder for fluidization action.
   * \param[in] e_den_lat pointer to the lattice used for the energy density interpolation
   * \param[in] e_den_background pointer to the background mapping between particle indices and the corresponding background
   * \param[in] energy_threshold minimum energy density required for a hadron to fluidize
   * \param[in] min_time minimum time to start fluidization,
   * \param[in] max_time
   * \param[in] fluid_cells
   */
  DynamicFluidizationFinder(RectangularLattice<EnergyMomentumTensor> *e_den_lat,
                            std::map<int32_t, double> *e_den_background,
                            double energy_threshold,
                            double min_time,
                            double max_time,
                            int fluid_cells)
      : e_den_lat_{*e_den_lat},
        background_{*e_den_background},
        energy_density_threshold_{energy_threshold},
        min_time_{min_time},
        max_time_{max_time},
        fluid_cells_{fluid_cells} {};

  /**
   * Find particles to fluidize, depending on the energy density around them.
   * \param[in] search_list List of candidate particles to fluidize
   * \param[in] dt Time step duration \unit{in fm}
   * \param[in] gcell_vol volume of the grid cell \unit{in fm³}
   * \param[in] beam_momentum unused
   */
  ActionList find_actions_in_cell(
      const ParticleList &search_list, double dt, const double gcell_vol,
      const std::vector<FourVector> &beam_momentum) const override;

  /// Ignore the neighbor search for fluidization
  ActionList find_actions_with_neighbors(
      const ParticleList &, const ParticleList &, const double,
      const std::vector<FourVector> &) const override {
    return {};
  }

  /// Ignore the surrounding searches for fluidization
  ActionList find_actions_with_surrounding_particles(
      const ParticleList &, const Particles &, double,
      const std::vector<FourVector> &) const override {
    return {};
  }

  /// No final actions after fluidizing
  ActionList find_final_actions(const Particles &, bool) const override {
    return {};
  }

 private:
  /**
   * Lattice where energy momentum tensor is computed
   * \note It is a reference so that it can be updated outside the class.
   */
  RectangularLattice<EnergyMomentumTensor> &e_den_lat_;
  /// Background energy density at positions of particles, using the id as key
  std::map<int32_t, double> &background_;
  /// Minimum energy density surrounding the particle to fluidize it
  const double energy_density_threshold_ = InputKeys::output_initialConditions_eDenThreshold.default_value();
  /// Minimum time (in lab frame) in fm to allow fluidization
  const double min_time_ = InputKeys::output_initialConditions_minTime.default_value();
  /// Maximum time (in lab frame) in fm to allow fluidization
  const double max_time_ = InputKeys::output_initialConditions_maxTime.default_value();
  /// Number of cells to interpolate the energy density
  const int fluid_cells_ = InputKeys::output_initialConditions_fluidCells.default_value();
  /**
   * Queue for future fluidizations, which will take place after the formation time of particles. Keys are particle indices and values are absolute formation time in the lab frame.
   * \note It must be mutable so that finder_actions_in_cell, a const method, can modify it.
   */
  mutable std::map<int32_t, double> queue_{};

  /**
   * Determine fluidization
   * \param[in] pdata particle to be checked for fluidization
   * \return whether energy density around pdata is high enough
   */
  bool above_threshold(const ParticleData &pdata) const;
};

/// Build energy momentum tensor
void build_fluidization_lattice(
    RectangularLattice<EnergyMomentumTensor> *e_den_lat, double t,
    const std::vector<Particles> &ensembles, const DensityParameters &dens_par);

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_DYNAMICFLUIDFINDER_H_
