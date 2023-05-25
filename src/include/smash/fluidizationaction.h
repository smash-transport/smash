/*
 *
 *    Copyright (c) 2023
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_HYPERSURFACECROSSINGACTION_H_
#define SRC_INCLUDE_SMASH_HYPERSURFACECROSSINGACTION_H_

#include <vector>

#include "action.h"
#include "actionfinderfactory.h"
#include "input_keys.h"

namespace smash {

/**
 * \ingroup action
 */
class FluidizationAction : public Action {
 public:
  /**
   * Construct action of fluidization.
   * \param[in] in_part Data of incoming particle.
   * \param[in] out_part Data of particles which surpass the energy density threshold
   */
  FluidizationAction(const ParticleData &in_part,
                             const ParticleData &out_part)
      : Action(in_part, out_part, time_until,
               ProcessType::HyperSurfaceCrossing) {}
  double get_total_weight() const override { return 0.0; };
  double get_partial_weight() const override { return 0.0; };
  void format_debug_output(std::ostream &out) const override {
    out << "Fluidization of " << incoming_particles_;
  }

  /**
   * Generate the final state of particles to be fluidized and removes them from the
   * evolution.
   */
  void generate_final_state() override;

  /** 
   * This function overrides Action::check_conservation that returns
   * the amount of energy density violation due to Pythia processes,
   * which is 0. here.
   *
   * \param[in] id_process process id only used for debugging output.
   * \return 0.
   */
  double check_conservation(uint32_t id_process) const override;
};

/**
 * \ingroup action
 */
class FluidizationActionsFinder : public ActionFinderInterface {
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
  explicit FluidizationFinder(const double min_energy=0.5, const double min_time=0, const double max_time=100) : energy_density_threshold_{min_energy}, min_time_{min_time}, max_time_{max_time} {
    std::array<double, 3> l{20., 20., 20.};
    std::array<int, 3> n{10, 10, 10};
    std::array<double, 3> origin{-10., -10., -10.};
    bool periodic = false;
    lat_ = make_unique<RectangularLattice<EnergyMomentumTensor>>(l, n, origin, periodic, LatticeUpdate::EveryTimestep);
  };
  FluidizationActionsFinder(RectangularLattice<EnergyMomentumTensor> *e_den_lat,
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

  /// Ignore the cell search for fluidization
  ActionList find_actions_in_cell(
      const ParticleList &, double, const double,
      const std::vector<FourVector> &) const override {
	  return {};
  }

  /**
   * Find particles to fluidize positions due to neighboring particles.
   * \param[in] search_list List of candidate particles to fluidize
   * \param[in] neighbors_list List of neighboring particles that contribute 
   * \param[in] dt duration of the current time step [fm]
   * \param[in] beam_momentum Unimportant?
   */
  ActionList find_actions_with_neighbors(
      const ParticleList &search_list, const ParticleList &neighbors_list, const double dt,
      const std::vector<FourVector> &beam_momentum) const override;

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
  /// Parameters for smearing particles in the lattice
  const DensityParameters dens_par_;
  /// Lattice where the energy density will be computed
  std::unique_ptr<RectangularLattice<EnergyMomentumTensor>> lat_;
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
  bool above_threshold(ParticleData &pdata) const;
  /// Build energy momentum tensor
  inline void build_lattice(const ParticleList &neighbors_list){
    update_lattice(lat_.get(), LatticeUpdate::EveryTimestep, DensityType::None, dens_par_, neighbors_list, false);
  }
};

/// Build energy momentum tensor
void build_fluidization_lattice(
    RectangularLattice<EnergyMomentumTensor> *e_den_lat, double t,
    const std::vector<Particles> &ensembles, const DensityParameters &dens_par);

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_HYPERSURFACECROSSINGACTION_H_
