/*
 *
 *    Copyright (c) 2024
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
    if (t_end < min_time_ || t0 > max_time_) {
      break;
    }

    const int32_t id = p.id();
    if (queue_.count(id)) {
      if (queue_[id] < t_end) {
        actions.emplace_back(
            std::make_unique<FluidizationAction>(p, p, queue_[id] - t0));
        queue_.erase(id);
      }
    } else {
      if (is_process_fluidizable(p.get_history().process_type)) {
        if (above_threshold(p)) {
          double fluidization_time =
              t0 + formation_time_fraction_ * (p.formation_time() - t0);
          if (fluidization_time >= t_end) {
            queue_.emplace(id, fluidization_time);
            continue;
          } else {
            double time_until = (1 - p.xsec_scaling_factor() <= really_small)
                                    ? 0
                                    : fluidization_time - t0;
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
  // value_at returns false if pdata is out of bounds
  const bool inside =
      energy_density_lattice_.value_at(pdata.position().threevec(), Tmunu);
  if (inside) {
    // If the particle is not in the map, the background evaluates to 0
    const double background =
        background_.count(pdata.id()) ? background_.at(pdata.id()) : 0;
    const double e_den_particles =
        Tmunu.boosted(Tmunu.landau_frame_4velocity())[0];
    if (e_den_particles + background >= energy_density_threshold_) {
      logg[LFluidization].debug()
          << "Fluidize " << pdata.id() << " with " << e_den_particles << "+"
          << background << " GeV/fm^3 at " << pdata.position().x0()
          << " fm, formed at " << pdata.formation_time() << " fm";
      return true;
    }
  }
  return false;
}

bool DynamicFluidizationFinder::is_process_fluidizable(
    const ProcessType &type) const {
  if (is_string_soft_process(type)) {
    return fluidizable_processes_
        [IncludedFluidizableProcesses::From_SoftString];
  }
  switch (type) {
    case ProcessType::Elastic:
      return fluidizable_processes_[IncludedFluidizableProcesses::From_Elastic];
      break;
    case ProcessType::Decay:
      return fluidizable_processes_[IncludedFluidizableProcesses::From_Decay];
      break;
    case ProcessType::TwoToOne:
    case ProcessType::TwoToTwo:
    case ProcessType::TwoToThree:
    case ProcessType::TwoToFour:
    case ProcessType::TwoToFive:
    case ProcessType::MultiParticleThreeMesonsToOne:
    case ProcessType::MultiParticleThreeToTwo:
    case ProcessType::MultiParticleFourToTwo:
    case ProcessType::MultiParticleFiveToTwo:
      return fluidizable_processes_
          [IncludedFluidizableProcesses::From_Inelastic];
      break;
    case ProcessType::StringHard:
      return fluidizable_processes_
          [IncludedFluidizableProcesses::From_HardString];
      break;
    default:
      return false;
  }
  return false;
}

}  // namespace smash
