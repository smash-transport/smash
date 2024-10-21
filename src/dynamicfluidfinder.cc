/*
 *
 *    Copyright (c) 2023
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/dynamicfluidfinder.h"

#include "smash/fluidizationaction.h"
#include "smash/logging.h"

namespace smash {
static constexpr int LFluidization = LogArea::HyperSurfaceCrossing::id;

ActionList DynamicFluidizationFinder::find_actions_in_cell(
    const ParticleList &search_list, double dt,
    [[maybe_unused]] const double gcell_vol,
    [[maybe_unused]] const std::vector<FourVector> &beam_momentum) const {
  ActionList actions;

  for (const ParticleData &p : search_list) {
    const double t0 = p.position().x0();
    const double t_end = t0 + dt;
    // Particles should not be removed before the nuclei collide, and after some
    // time max_time_ there won't be any fluidization, so this saves resources
    if (t0 < min_time_ || t_end > max_time_)
      break;

    const int32_t id = p.id();
    if (queue_.count(id)) {
      if (queue_[id] < t_end) {
        actions.emplace_back(
            std::make_unique<FluidizationAction>(p, p, queue_[id] - t0));
        queue_.erase(id);
      }
    } else {
      const auto process_type = p.get_history().process_type;
      if (process_type == ProcessType::Decay ||
          is_string_soft_process(process_type) ||
          process_type == ProcessType::StringHard) {
        if (above_threshold(p)) {
          double formation = p.formation_time();
          if (formation >= t_end) {
            queue_.emplace(id, formation);
          } else {
            double time_until = (1 - p.xsec_scaling_factor() <= really_small)
                                    ? 0
                                    : formation - t0;
            // todo(hirayama): add option to not wait for formation (also
            // leading hadrons) or scale it
            actions.emplace_back(
                std::make_unique<FluidizationAction>(p, p, time_until));
          }
        }
      }
    }
  }  // search_list loop
  return actions;
}

bool DynamicFluidizationFinder::above_threshold(
    const ParticleData &pdata) const {
  EnergyMomentumTensor Tmunu;
  // value_at returns false if pdata is out of bounds, this is desirable here
  bool inside =
      energy_density_lattice_.value_at(pdata.position().threevec(), Tmunu);
  if (inside) {
    // If the particle is not in the map, the background evaluates to 0.
    const double background = background_[pdata.id()];
    const double e_den_particles =
        Tmunu.boosted(Tmunu.landau_frame_4velocity())[0];
    if (e_den_particles + background >= energy_density_threshold_) {
      logg[LFluidization].debug()
          << "Fluidize " << pdata.id() << " with " << e_den_particles
          << " and background " << background << " GeV/fm^3 at "
          << pdata.formation_time();
      return true;
    }
  }
  return false;
}

void build_fluidization_lattice(
    RectangularLattice<EnergyMomentumTensor> *energy_density_lattice,
    const double t, const std::vector<Particles> &ensembles,
    const DensityParameters &dens_par) {
  if (energy_density_lattice == nullptr) {
    return;
  }
  // In most scenarios where dynamic fluidization is applicable, t > 20 fm is
  // dilute enough to not need a very fine lattice.
  if (t > 20) {
    std::array<double, 3> new_l{2 * t, 2 * t, 2 * t};
    std::array<double, 3> new_orig{-t, -t, -t};
    energy_density_lattice->reset_and_resize(new_l, new_orig);
    logg[LFluidization].debug() << "Lattice resizing at " << t;
  }

  update_lattice(energy_density_lattice, LatticeUpdate::EveryTimestep,
                 DensityType::Hadron, dens_par, ensembles, false);
}

}  // namespace smash
