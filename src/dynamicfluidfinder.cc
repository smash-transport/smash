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
  std::vector<ActionPtr> actions;

  for (const ParticleData &p : search_list) {
    double t0 = p.position().x0();
    double t_end = t0 + dt;  // Time at the end of timestep
    // Particles should not be removed before the nuclei collide, and after some
    // time max_time_ there won't be any fluidization, so this saves resources
    if (t0 < min_time_ || t_end > max_time_)
      break;

    const int32_t id = p.id();
    if (queue_.count(id)) {
      if (queue_.at(id) < t_end) {
        ActionPtr action =
            std::make_unique<FluidizationAction>(p, p, queue_.at(id) - t0);
        actions.emplace_back(std::move(action));
        queue_.erase(id);
      }
    } else {
      const auto proc = p.get_history().process_type;
      if (proc == ProcessType::Decay || is_string_soft_process(proc) ||
          proc == ProcessType::StringHard) {
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
            ActionPtr action =
                std::make_unique<FluidizationAction>(p, p, time_until);
            actions.emplace_back(std::move(action));
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
  bool fluidize = false;

  double background = background_[pdata.id()];
  // value_at returns false if pdata is out of bounds, this is desirable here
  bool inside = e_den_lat_.value_at(pdata.position().threevec(), Tmunu);
  if (inside) {
    double e_den_particles = Tmunu.boosted(Tmunu.landau_frame_4velocity())[0];
    if (e_den_particles + background >= energy_density_threshold_) {
      fluidize = true;
      logg[LFluidization].debug()
          << "Fluidize " << pdata.id() << " with " << e_den_particles
          << " and background " << background << " GeV/fm^3 at "
          << pdata.formation_time();
    }
  }
  return fluidize;
}

void build_fluidization_lattice(
    RectangularLattice<EnergyMomentumTensor> *e_den_lat, const double t,
    const std::vector<Particles> &ensembles, const DensityParameters &dens_par) {
  if (t > 20) {
    std::array<double, 3> new_l{2 * t, 2 * t, 2 * t};
    std::array<double, 3> new_orig{-t, -t, -t};
    e_den_lat->reset_and_resize(new_l, new_orig);
  }

  update_lattice(e_den_lat, LatticeUpdate::EveryTimestep, DensityType::Hadron,
                 dens_par, ensembles, false);
}

} // namespace smash
