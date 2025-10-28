/*
 *
 *    Copyright (c) 2022-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_DYNAMICFLUIDFINDER_H_
#define SRC_INCLUDE_SMASH_DYNAMICFLUIDFINDER_H_

#include <limits>
#include <map>
#include <vector>

#include "actionfinderfactory.h"
#include "icparameters.h"
#include "input_keys.h"

namespace smash {

/**
 * \ingroup action
 * Finder for dynamic fluidizations.
 * Loops through all particles and checks if they reach the energy density
 * threshold. This happens at the end of every time step for all hadrons that
 * participated in a fluidizable process. For string fragmentation products,
 * fluidization happens only after the fraction of formation time of the
 * particle, controlled by \key Formation_Time_Fraction.
 */
class DynamicFluidizationFinder : public ActionFinderInterface {
 public:
  /**
   * Construct finder for fluidization action.
   * \param[in] lattice Lattice used for the energy density interpolation
   * \param[in] background Background map between particle indices and the
   * corresponding background energy density
   * \param[in] ic_params Parameters for dynamic fluidization
   *
   * \note \c energy_density_lattice and \c energy_density_background are both
   * "in" parameters because the class stores references, but the values are non
   * constant.
   */
  DynamicFluidizationFinder(
      const RectangularLattice<EnergyMomentumTensor> &lattice,
      const std::map<int32_t, double> &background,
      const InitialConditionParameters &ic_params)
      : energy_density_lattice_{lattice},
        background_{background},
        energy_density_threshold_{ic_params.energy_density_threshold.value()},
        min_time_{ic_params.min_time.value()},
        max_time_{ic_params.max_time.value()},
        formation_time_fraction_{ic_params.formation_time_fraction.value()},
        smearing_kernel_at_0_{ic_params.smearing_kernel_at_0.value()},
        fluid_cells_{ic_params.num_fluid_cells.value()},
        fluidizable_processes_{ic_params.fluidizable_processes.value()},
        delay_initial_elastic_{ic_params.delay_initial_elastic.value()} {};

  /**
   * Find particles to fluidize, depending on the energy density around them.
   * \param[in] search_list List of candidate particles to fluidize.
   * \param[in] dt Time step duration \unit{in fm}.
   * \param[in] gcell_vol Volume of the grid cell \unit{in fm³},
   * \param[in] beam_momentum Unused.
   * \return List of fluidization actions that will happen at the corresponding
   * formation time.
   */
  ActionList find_actions_in_cell(
      const ParticleList &search_list, double dt, double gcell_vol,
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

  /**
   * Prepare corona particles left in the IC for the afterburner.
   *
   * - Remove core particles.
   * - Move corona particles back to their last interaction point so they start
   *   at the correct spacetime position for afterburner rescattering with
   *   hadrons sampled from the fluid.
   * - Do not modify spectators.
   *
   * \param[in] search_list All particles at the end of the simulation.
   * \return ActionList with removals and backpropagation steps.
   *
   * \note The backpropagation is encoded as two FreeForAll actions:
   *       remove the corona particle, then add it at the target spacetime
   *       point. This will appear in the Collisions output.
   */
  ActionList find_final_actions(const Particles &search_list) const override;

  /**
   * Determine if fluidization condition is satisfied.
   *
   * \param[in] pdata Particle to be checked for fluidization.
   * \return Whether energy density around pdata is high enough.
   */
  bool above_threshold(const ParticleData &pdata) const;

  /**
   * Checks if a given process type is in \ref fluidizable_processes_. In
   * particular, initially sampled hadrons are not fluidizable and have
   * ProcessType::None, which falls in the default case for the switch.
   *
   * \return whether the process is fluidizable
   */
  bool is_process_fluidizable(const HistoryData &history) const;

 private:
  /**
   * Lattice where energy momentum tensor is computed
   *
   * \note It must be a reference so that it can be updated outside the class,
   * without creating a new Finder object.
   */
  const RectangularLattice<EnergyMomentumTensor> &energy_density_lattice_;
  /**
   * Background energy density at positions of particles, using the id as key
   *
   * \note It is a reference so that it can be updated outside the class, e.g.
   * by an external manager using SMASH as a library.
   */
  const std::map<int32_t, double> &background_;
  /// Minimum energy density surrounding the particle to fluidize it
  const double energy_density_threshold_ = smash_NaN<double>;
  /// Minimum time (in lab frame) in fm to allow fluidization
  const double min_time_ = smash_NaN<double>;
  /// Maximum time (in lab frame) in fm to allow fluidization
  const double max_time_ = smash_NaN<double>;
  /// Fraction of formation time after which a particles can fluidize
  const double formation_time_fraction_ = smash_NaN<double>;
  /// Smearing kernel at the position of the particle of interest
  const double smearing_kernel_at_0_ = smash_NaN<double>;
  /// Number of cells to interpolate the energy density
  const int fluid_cells_ = smash_NaN<int>;
  /// Processes that create a fluidizable particle
  const FluidizableProcessesBitSet fluidizable_processes_;
  /// Whether the first elastic interaction of an initial nucleon is fluidizable
  const bool delay_initial_elastic_ = true;

  /// Accumulated number of core particles
  mutable int particles_in_core_ = 0;
  /// Accumulated energy of core particles
  mutable double energy_in_core_ = 0.;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_DYNAMICFLUIDFINDER_H_
