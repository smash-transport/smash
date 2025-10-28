/*
 *
 *    Copyright (c) 2023-2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/dynamicfluidfinder.h"

#include "smash/fluidizationaction.h"
#include "smash/freeforallaction.h"
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
    const double t_creation = p.begin_formation_time();
    const double t_end = t0 + dt;
    if (t_end < min_time_ || t0 > max_time_) {
      break;
    }
    if (p.is_core()) {
      continue;
    }
    const double fluidization_time =
        t_creation +
        formation_time_fraction_ * (p.formation_time() - t_creation);
    if (fluidization_time >= t_end) {
      continue;
    }
    if (!is_process_fluidizable(p.get_history())) {
      continue;
    }
    if (above_threshold(p)) {
      double time_until = (1 - p.xsec_scaling_factor() <= really_small)
                              ? 0
                              : std::max(fluidization_time - t0, 0.);
      actions.emplace_back(
          std::make_unique<FluidizationAction>(p, p, time_until));
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
    if (e_den_particles + background >=
        energy_density_threshold_ + pdata.pole_mass() * smearing_kernel_at_0_) {
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
    const HistoryData &history) const {
  const ProcessType &type = history.process_type;
  if (is_string_soft_process(type)) {
    return fluidizable_processes_
        [IncludedFluidizableProcesses::From_SoftString];
  }
  switch (type) {
    case ProcessType::Elastic:
      if (history.collisions_per_particle == 1 && delay_initial_elastic_) {
        return false;
      }
      return fluidizable_processes_[IncludedFluidizableProcesses::From_Elastic];
    case ProcessType::Decay:
      return fluidizable_processes_[IncludedFluidizableProcesses::From_Decay];
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
    case ProcessType::StringHard:
      return fluidizable_processes_
          [IncludedFluidizableProcesses::From_HardString];
    default:
      return false;
  }
}

ActionList DynamicFluidizationFinder::find_final_actions(
    const Particles &search_list) const {
  ActionList actions;
  const bool are_there_core_particles =
      std::any_of(search_list.begin(), search_list.end(),
                  [](const ParticleData &p) { return p.is_core(); });
  if (are_there_core_particles) {
    for (auto &p : search_list) {
      if (p.is_core()) {
        actions.emplace_back(std::make_unique<FreeforallAction>(
            ParticleList{p}, ParticleList{}, p.position().x0()));
        particles_in_core_++;
        energy_in_core_ += p.momentum()[0];
      }
    }
  } else {
    for (auto &original : search_list) {
      if (original.get_history().collisions_per_particle == 0) {
        // Spectators are not propagated back.
        continue;
      }
      const double t = original.position().x0();
      double corona_time = original.get_history().time_last_collision;
      if (t == corona_time) {
        continue;
      }
      // This prevents particles from decays or strings from ending up at the
      // same position
      corona_time += 0.01;
      const ThreeVector r = original.position().threevec() -
                            (t - corona_time) * original.velocity();
      ParticleData backpropagated{original.type()};
      backpropagated.set_4position(FourVector(corona_time, r));
      backpropagated.set_4momentum(original.momentum());
      // This is done so no further actions can be found for this particle
      backpropagated.set_formation_time(t);
      backpropagated.set_cross_section_scaling_factor(0.0);
      actions.emplace_back(std::make_unique<FreeforallAction>(
          ParticleList{original}, ParticleList{}, t));
      actions.emplace_back(std::make_unique<FreeforallAction>(
          ParticleList{}, ParticleList{backpropagated}, corona_time));
    }
    logg[LFluidization].info()
        << particles_in_core_ << " particles were part of the core with energy "
        << energy_in_core_ << " GeV.";
  }
  return actions;
}

}  // namespace smash
