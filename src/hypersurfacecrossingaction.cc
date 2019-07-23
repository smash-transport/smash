/*
 *
 *    Copyright (c) 2017-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/hypersurfacecrossingaction.h"

#include "smash/logging.h"
#include "smash/particledata.h"
#include "smash/particles.h"

namespace smash {

void HypersurfacecrossingAction::generate_final_state() {
  const auto &log = logger<LogArea::HyperSurfaceCrossing>();
  log.debug("Process: Hypersurface Crossing. ");

  ParticleList empty_list;

  // check that there is only 1 incoming particle which lies on the hypersurface
  assert(incoming_particles_.size() == 1);
  assert(abs(incoming_particles[0].position().tau() - prop_time_) <=
         really_small);

  // Return empty list because we want to remove the incoming particle
  outgoing_particles_ = empty_list;
}

ActionList HyperSurfaceCrossActionsFinder::find_actions_in_cell(
    const ParticleList &plist, double t_max) const {
  std::vector<ActionPtr> actions;

  for (const ParticleData &p : plist) {
    ParticleData pdata_before_propagation = p;
    ParticleData pdata_after_propagation = p;  // Will receive updated position
    double t0 = p.position().x0();
    double t_end = t0 + t_max;  // Time at the end of timestep

    // We don't want to remove particles before the nuclei have interacted
    if (t_end < 0.0) {
      continue;
    }

    // propagate particles to position where they would be after t_max
    const ThreeVector &v = p.velocity();
    const FourVector distance = FourVector(0.0, v * t_max);
    FourVector position = p.position() + distance;
    position.set_x0(t_end);
    // update coordinates to the position corresponding to t_max
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

/**
 * Determine whether or not particle crosses hypersurface of given proper
 * time during timestepless propagation.
 *
 * \param[in] pdata_before_propagation data of the particle of interest before
 * it is propagated.
 * \param[in] pdata_after_propagation data of the particle of interest after it
 * was propagated.
 * \param[in] tau Proper time of hypersurface
 *
 * \return Does particle cross hypersurface?
 */
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

/**
 * Find the coordinates at which a particle crosses a hypersurface of constant
 * proper time during timestepless propagation.
 *
 * \param[in] pdata_before_propagation data of the particle of interest before
 * it is propagated.
 * \param[in] pdata_after_propagation data of the particle of interest after it
 * was propagated.
 * \param[in] tau Proper time of hypersurface
 *
 * \return Coordinates when crossing hypersurface
 */
FourVector HyperSurfaceCrossActionsFinder::coordinates_on_hypersurface(
    ParticleData &pdata_before_propagation,
    ParticleData &pdata_after_propagation, const double tau) const {
  // find t and z at start of propagation
  const double t1 = pdata_before_propagation.position().x0();
  const double z1 = pdata_before_propagation.position().x3();

  // find t and z after propagation
  const double t2 = pdata_after_propagation.position().x0();
  const double z2 = pdata_after_propagation.position().x3();

  // find slope and intercept of linear function
  const double m = (z2 - z1) / (t2 - t1);
  const double n = z1 - m * t1;

  // The equation to solve is a quadratic equation which provides two solutions,
  // the latter is usually out of the t-interval we are looking at.
  const double sol1 = n * m / (1 - m * m) +
                      std::sqrt((1 - m * m) * tau * tau + n * n) / (1 - m * m);
  const double sol2 = n * m / (1 - m * m) -
                      std::sqrt((1 - m * m) * tau * tau + n * n) / (1 - m * m);

  SMASH_UNUSED(sol2);  // only used in DEBUG output
  // std::cout << t1 << " " << sol1 << " " << t2 << '\n';
  // assert((sol1 >= t1 && sol1 <= t2));
  assert(!(sol2 >= t1 && sol2 <= t2));

  // Propagate to point where hypersurface is crossed
  const ThreeVector v = pdata_before_propagation.velocity();
  const FourVector distance = FourVector(0.0, v * (sol1 - t1));
  FourVector crossing_position = pdata_before_propagation.position() + distance;
  crossing_position.set_x0(sol1);

  return crossing_position;
}

}  // namespace smash
