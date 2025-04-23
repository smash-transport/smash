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
    if (!is_process_fluidizable(p.get_history().process_type)) {
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
        energy_density_threshold_ +
            pdata.pole_mass() * smearing_kernel_at_0_) {
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
    const Particles &search_list, [[maybe_unused]] bool only_res) const {
  ActionList actions;
  static int fluidized_particles = 0;
  static double energy = 0.;
  bool print_warning = true;
  for (auto &p : search_list) {
    if (p.is_core()) {
      actions.emplace_back(std::make_unique<FreeforallAction>(
          ParticleList{p}, ParticleList{}, p.position().x0()));
      fluidized_particles++;
      energy += p.momentum()[0];
      print_warning = false;  // print only at the final call
    }
  }
  if (print_warning) {
    logg[LFluidization].info() << fluidized_particles
                               << " particles were part of the core at the end"
                                  " of the evolution with energy "
                               << energy << " GeV.";
  }
  return actions;
}

}  // namespace smash
