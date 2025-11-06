/*
 *
 *    Copyright (c) 2019-2020,2022-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/hypersurfacecrossingfinder.h"

#include "smash/fluidizationaction.h"
#include "smash/logging.h"

namespace smash {
static constexpr int LHyperSurfaceCrossing = LogArea::HyperSurfaceCrossing::id;

ActionList HyperSurfaceCrossActionsFinder::find_actions_in_cell(
    const ParticleList &plist, double dt, const double,
    const std::vector<FourVector> &beam_momentum) const {
  ActionList actions;

  for (const ParticleData &p : plist) {
    ParticleData pdata_before_propagation = p;
    ParticleData pdata_after_propagation = p;  // Will receive updated position
    double t0 = p.position().x0();
    double t_end = t0 + dt;  // Time at the end of timestep

    // We don't want to remove particles before the nuclei have interacted
    // because those would not yet be part of the newly-created medium.
    if (t_end < 0.0) {
      continue;
    }
    if (p.is_core()) {
      continue;
    }

    // For frozen Fermi motion:
    // Fermi momenta are only applied if particles interact. The particle
    // properties p.velocity() and p.momentum() already contain the values
    // corrected by Fermi motion, but those particles that have not yet
    // interacted are propagated along the beam-axis with v = (0, 0, beam_v)
    // (and not with p.velocity()).
    // To identify the corresponding hypersurface crossings the finding for
    // those paricles without prior interactions has to be performed with
    // v = vbeam instead of p.velocity().
    // Note: The beam_momentum vector is empty in case frozen Fermi motion is
    // not applied.
    const bool no_prior_interactions =
        (static_cast<uint64_t>(p.id()) <                  // particle from
         static_cast<uint64_t>(beam_momentum.size())) &&  // initial nucleus
        (p.get_history().collisions_per_particle == 0);
    ThreeVector v;
    if (no_prior_interactions) {
      const FourVector vbeam = beam_momentum[p.id()];
      v = vbeam.velocity();
    } else {
      v = p.velocity();
    }

    // propagate particles to position where they would be at the end of the
    // time step (after dt)
    const FourVector distance = FourVector(0.0, v * dt);
    FourVector position = p.position() + distance;
    position.set_x0(t_end);
    // update coordinates to the position corresponding to t_end
    pdata_after_propagation.set_4position(position);

    bool hypersurface_is_crossed = crosses_hypersurface(
        pdata_before_propagation, pdata_after_propagation, prop_time_);

    /*
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
      const double rapidity = p.rapidity();
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

    if (hypersurface_is_crossed && is_within_y_cut && is_within_pT_cut) {
      // Get exact coordinates where hypersurface is crossed
      FourVector crossing_position = coordinates_on_hypersurface(
          pdata_before_propagation, pdata_after_propagation, prop_time_);

      double time_until_crossing = crossing_position[0] - t0;

      ParticleData outgoing_particle(p);
      outgoing_particle.set_4position(crossing_position);
      ActionPtr action = std::make_unique<FluidizationAction>(
          p, outgoing_particle, time_until_crossing);
      actions.emplace_back(std::move(action));
    }
  }
  return actions;
}

bool HyperSurfaceCrossActionsFinder::crosses_hypersurface(
    ParticleData &pdata_before_propagation,
    ParticleData &pdata_after_propagation, const double tau) const {
  bool hypersurface_is_crossed = false;
  const bool t_greater_z_before_prop =
      (std::fabs(pdata_before_propagation.position().x0()) >
               std::fabs(pdata_before_propagation.position().x3())
           ? true
           : false);
  const bool t_greater_z_after_prop =
      (std::fabs(pdata_after_propagation.position().x0()) >
               std::fabs(pdata_after_propagation.position().x3())
           ? true
           : false);

  if (t_greater_z_before_prop && t_greater_z_after_prop) {
    // proper time before and after propagation
    const double tau_before = pdata_before_propagation.hyperbolic_time();
    const double tau_after = pdata_after_propagation.hyperbolic_time();

    if (tau_before <= tau && tau <= tau_after) {
      hypersurface_is_crossed = true;
    }
  } else if (!t_greater_z_before_prop && t_greater_z_after_prop) {
    // proper time after propagation
    const double tau_after = pdata_after_propagation.hyperbolic_time();
    if (tau_after >= tau) {
      hypersurface_is_crossed = true;
    }
  }

  return hypersurface_is_crossed;
}

FourVector HyperSurfaceCrossActionsFinder::coordinates_on_hypersurface(
    ParticleData &pdata_before_propagation,
    ParticleData &pdata_after_propagation, const double tau) const {
  // find t and z at start of propagation
  const double t1 = pdata_before_propagation.position().x0();
  const double z1 = pdata_before_propagation.position().x3();

  // find t and z after propagation
  const double t2 = pdata_after_propagation.position().x0();
  const double z2 = pdata_after_propagation.position().x3();

  // find slope and intercept of linear function that describes propagation on
  // straight line
  const double m = (z2 - z1) / (t2 - t1);
  const double n = z1 - m * t1;

  // The equation to solve is a quadratic equation which provides two solutions,
  // the latter is usually out of the t-interval we are looking at.
  const double sol1 = n * m / (1 - m * m) +
                      std::sqrt((1 - m * m) * tau * tau + n * n) / (1 - m * m);
  [[maybe_unused]] const double sol2 =  // only used in DEBUG output
      n * m / (1 - m * m) -
      std::sqrt((1 - m * m) * tau * tau + n * n) / (1 - m * m);

  assert((sol1 >= t1 && sol1 <= t2));
  assert(!(sol2 >= t1 && sol2 <= t2));

  // Propagate to point where hypersurface is crossed
  const ThreeVector v = pdata_before_propagation.velocity();
  const FourVector distance = FourVector(0.0, v * (sol1 - t1));
  FourVector crossing_position = pdata_before_propagation.position() + distance;
  crossing_position.set_x0(sol1);

  return crossing_position;
}

void HyperSurfaceCrossActionsFinder::warn_if_some_particles_did_not_cross(
    const size_t number_of_particles, bool impose_kinematic_cut) {
  if (number_of_particles != 0 && !impose_kinematic_cut) {
    logg[LHyperSurfaceCrossing].warn(
        "End time might be too small for initial conditions output. "
        "Hypersurface has not yet been crossed by ",
        number_of_particles, " particle(s).");
  }
}

}  // namespace smash
