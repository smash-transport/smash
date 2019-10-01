/*
 *
 *    Copyright (c) 2019-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/hypersurfacecrossingaction.h"

#include "smash/logging.h"
#include "smash/particledata.h"
#include "smash/particles.h"
#include "smash/quantumnumbers.h"

namespace smash {
inline constexpr int HyperSurfaceCrossing = LogArea::HyperSurfaceCrossing::id;

void HypersurfacecrossingAction::generate_final_state() {
  logg[HyperSurfaceCrossing].debug("Process: Hypersurface Crossing. ");

  ParticleList empty_list;

  // check that there is only 1 incoming particle
  assert(incoming_particles_.size() == 1);

  // Return empty list because we want to remove the incoming particle
  outgoing_particles_ = empty_list;
}

void HypersurfacecrossingAction::check_conservation(
    const uint32_t id_process) const {
  QuantumNumbers before(incoming_particles_);
  QuantumNumbers after(outgoing_particles_);
  if (before == after) {
    // Conservation laws should not be conserved since particles are removed
    // from the evolution
    throw std::runtime_error(
        "Conservation laws conserved in the hypersurface "
        "crossing action. Particle was not properly removed in process: " +
        std::to_string(id_process));
  }

  if (outgoing_particles_.size() != 0) {
    throw std::runtime_error(
        "Particle was not removed successfully in "
        "hypersurface crossing action.");
  }
}

ActionList HyperSurfaceCrossActionsFinder::find_actions_in_cell(
    const ParticleList &plist, double dt, const double,
    const std::vector<FourVector> beam_momentum) const {
  std::vector<ActionPtr> actions;

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

    // For frozen Fermi motion:
    // Fermi momenta are only applied if particles interact. The particle
    // properties p.velocity() and p.momentum() already contain the values
    // corrected by Fermi motion, but those particles that have not yet
    // interacted are propagated along the beam-axis with v = (0, 0, beam_v)
    // (and not with p.velocity()).
    // To identify the corresponding hypersurface crossings the finding for
    // those paricles without prior interactions has to be performed with
    // v = vbeam instead of p.velcocity().
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

    if (hypersurface_is_crossed) {
      // Get exact coordinates where hypersurface is crossed
      FourVector crossing_position = coordinates_on_hypersurface(
          pdata_before_propagation, pdata_after_propagation, prop_time_);

      double time_until_crossing = crossing_position[0] - t0;

      ParticleData outgoing_particle(p);
      outgoing_particle.set_4position(crossing_position);
      ActionPtr action = make_unique<HypersurfacecrossingAction>(
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
      (fabs(pdata_before_propagation.position().x0()) >
               fabs(pdata_before_propagation.position().x3())
           ? 1
           : 0);
  const bool t_greater_z_after_prop =
      (fabs(pdata_after_propagation.position().x0()) >
               fabs(pdata_after_propagation.position().x3())
           ? 1
           : 0);

  if (t_greater_z_before_prop && t_greater_z_after_prop) {
    // proper time before and after propagation
    const double tau_before = pdata_before_propagation.position().tau();
    const double tau_after = pdata_after_propagation.position().tau();

    if (tau_before <= tau && tau <= tau_after) {
      hypersurface_is_crossed = true;
    }
  } else if (!t_greater_z_before_prop && t_greater_z_after_prop) {
    // proper time after propagation
    const double tau_after = pdata_after_propagation.position().tau();
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
  const double sol2 = n * m / (1 - m * m) -
                      std::sqrt((1 - m * m) * tau * tau + n * n) / (1 - m * m);

  SMASH_UNUSED(sol2);  // only used in DEBUG output
  assert((sol1 >= t1 && sol1 <= t2));
  assert(!(sol2 >= t1 && sol2 <= t2));

  // Propagate to point where hypersurface is crossed
  const ThreeVector v = pdata_before_propagation.velocity();
  const FourVector distance = FourVector(0.0, v * (sol1 - t1));
  FourVector crossing_position = pdata_before_propagation.position() + distance;
  crossing_position.set_x0(sol1);

  return crossing_position;
}

}  // namespace smash
