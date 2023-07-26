/*
 *
 *    Copyright (c) 2023
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/fluidizationaction.h"

#include "smash/logging.h"
#include "smash/quantumnumbers.h"

namespace smash {
static constexpr int LFluidization = LogArea::HyperSurfaceCrossing::id;

void HypersurfacecrossingAction::generate_final_state() {
  logg[LFluidization].debug("Process: Fluidization. ");

  // check that there is only 1 incoming particle
  assert(incoming_particles_.size() == 1);

  // Return empty list because we want to remove the incoming particle
  outgoing_particles_ = {};
}

void FluidizationAction::check_conservation(
    const uint32_t id_process) const {
  QuantumNumbers before(incoming_particles_);
  QuantumNumbers after(outgoing_particles_);
  if (before == after) {
    throw std::runtime_error(
        "Conservation laws obeyed during fluidization, which should not happen "
        "as particles are removed. Particle was not properly removed in "
        "process: " +
        std::to_string(id_process));
  }

  if (outgoing_particles_.size() != 0) {
    throw std::runtime_error(
        "Particle was not removed successfully in fluidization action.");
  }
}

ActionList FluidizationActionFinder::find_actions_with_neighbors(
      const ParticleList &search_list, const ParticleList &neighbors_list, const double dt,
      const std::vector<FourVector> &) const {
  std::vector<ActionPtr> actions;

  build_lattice(neighbors_list);

  for (const ParticleData &p : search_list) {
    double t0 = p.position().x0();
    double t_end = t0 + dt;  // Time at the end of timestep

    // Particles should not be removed before the nuclei collide, and after some 
    // time max_time_ there won't be any fluidization
    if (t_end < min_time_ || t0 > max_time_) {
      break;
    }

    bool fluidize = above_threshold(p, neighbors_list);

    /*
       This is from HypersurfacecrossingAction, not sure if will be needed.

       If rapidity or transverse momentum cut is to be employed; check if
       particles are within the relevant region
       Implementation explanation: The default for both cuts is 0.0, as a cut at
       0 implies that not a single particle contributes to the initial
       conditions. If the user specifies a value of 0.0 in the config, SMASH
       crashes with a corresponding error message. The same applies to negtive
       values.
    */
    bool is_within_y_cut = true;
    // Check whether particle is in desired rapidity range
    if (rap_cut_ > 0.0) {
      const double rapidity =
          0.5 * std::log((p.momentum().x0() + p.momentum().x3()) /
                         (p.momentum().x0() - p.momentum().x3()));
      if (std::fabs(rapidity) > rap_cut_) {
        is_within_y_cut = false;
      }
    }
    bool is_within_pT_cut = true;
    // Check whether particle is in desired pT range
    if (pT_cut_ > 0.0) {
      const double transverse_momentum =
          std::sqrt(p.momentum().x1() * p.momentum().x1() +
                    p.momentum().x2() * p.momentum().x2());
      if (transverse_momentum > pT_cut_) {
        is_within_pT_cut = false;
      }
    }

    if (fluidize && is_within_y_cut && is_within_pT_cut) {
      // for now the action should be performed at the start of the time step
      ActionPtr action = std::make_unique<FluidizationAction>(p, 0);
      actions.emplace_back(std::move(action));
    }
  }  // search_list loop
  return actions;
}

bool above_threshold(ParticleData &pdata){
  EnergyMomentumTensor Tmunu;
  bool fluidize = false;
  // value_at returns false if pdata is out of bounds, this is desirable here
  if(lat_.value_at(pdata.position().threevec(), Tmunu)){
      if (Tmunu.boosted(Tmunu.landau_frame_4velocity())[0] >= energy_density_threshold_)
	     fluidize = true; 
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

}  // namespace smash
