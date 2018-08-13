/*
 *
 *    Copyright (c) 2013-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/propagation.h"

#include "smash/boxmodus.h"
#include "smash/collidermodus.h"
#include "smash/listmodus.h"
#include "smash/logging.h"
#include "smash/spheremodus.h"

namespace smash {

double calc_hubble(double time, const ExpansionProperties &metric) {
  double h;  // Hubble parameter

  switch (metric.mode_) {
    case ExpansionMode::NoExpansion:
      h = 0.;
      break;
    case ExpansionMode::MasslessFRW:
      h = metric.b_ / (2 * (metric.b_ * time + 1));
      break;
    case ExpansionMode::MassiveFRW:
      h = 2 * metric.b_ / (3 * (metric.b_ * time + 1));
      break;
    case ExpansionMode::Exponential:
      h = metric.b_ * time;
      break;
    default:
      h = 0.;
  }

  return h;
}

double propagate_straight_line(Particles *particles, double to_time,
                               const std::vector<FourVector> &beam_momentum) {
  const auto &log = logger<LogArea::Propagation>();
  bool negative_dt_error = false;
  double dt = 0.0;
  for (ParticleData &data : *particles) {
    const double t0 = data.position().x0();
    dt = to_time - t0;
    if (dt < 0.0 && !negative_dt_error) {
      // Print error message once, not for every particle
      negative_dt_error = true;
      log.error("propagate_straight_line - negative dt = ", dt);
    }
    assert(dt >= 0.0);
    /* "Frozen Fermi motion": Fermi momenta are only used for collisions,
     * but not for propagation. This is done to avoid nucleus flying apart
     * even if potentials are off. Initial nucleons before the first collision
     * are propagated only according to beam momentum.
     * Initial nucleons are distinguished by data.id() < the size of
     * beam_momentum, which is by default zero except for the collider modus
     * with the fermi motion == frozen.
     * todo(m. mayer): improve this condition (see comment #11 issue #4213)*/
    assert(data.id() >= 0);
    const bool avoid_fermi_motion =
        (static_cast<uint64_t>(data.id()) <
         static_cast<uint64_t>(beam_momentum.size())) &&
        (data.get_history().collisions_per_particle == 0);
    ThreeVector v;
    if (avoid_fermi_motion) {
      const FourVector vbeam = beam_momentum[data.id()];
      v = vbeam.velocity();
    } else {
      v = data.velocity();
    }
    const FourVector distance = FourVector(0.0, v * dt);
    log.debug("Particle ", data, " motion: ", distance);
    FourVector position = data.position() + distance;
    position.set_x0(to_time);
    data.set_4position(position);

    // If particle is formed reset cross_section_scaling_factor
    if (data.formation_time() < to_time) {
      data.set_cross_section_scaling_factor(1.0);
    }
  }
  return dt;
}

void expand_space_time(Particles *particles,
                       const ExperimentParameters &parameters,
                       const ExpansionProperties &metric) {
  const auto &log = logger<LogArea::Propagation>();
  const double dt = parameters.labclock.timestep_duration();
  for (ParticleData &data : *particles) {
    // Momentum and position modification to ensure appropriate expansion
    const double h = calc_hubble(parameters.labclock.current_time(), metric);
    FourVector delta_mom = FourVector(0.0, h * data.momentum().threevec() * dt);
    FourVector expan_dist =
        FourVector(0.0, h * data.position().threevec() * dt);

    log.debug("Particle ", data, " expansion motion: ", expan_dist);
    // New position and momentum
    FourVector position = data.position() + expan_dist;
    FourVector momentum = data.momentum() - delta_mom;

    // set the new momentum and position variables
    data.set_4position(position);
    data.set_4momentum(momentum);
    // force the on shell condition to ensure correct energy
    data.set_4momentum(data.pole_mass(), data.momentum().threevec());
  }
}

void update_momenta(Particles *particles, double dt, const Potentials &pot,
                    RectangularLattice<ThreeVector> *UB_grad_lat,
                    RectangularLattice<ThreeVector> *UI3_grad_lat) {
  // Copy particles before propagation to calculate potentials from them
  const ParticleList plist = particles->copy_to_vector();

  const auto &log = logger<LogArea::Propagation>();
  bool possibly_use_lattice =
      (pot.use_skyrme() ? (UB_grad_lat != nullptr) : true) &&
      (pot.use_symmetry() ? (UI3_grad_lat != nullptr) : true);
  ThreeVector dUB_dr, dUI3_dr;
  double min_time_scale = std::numeric_limits<double>::infinity();

  for (ParticleData &data : *particles) {
    // Only baryons will be affected by the potentials
    if (!data.is_baryon()) {
      continue;
    }
    const auto scale = pot.force_scale(data.type());
    const ThreeVector r = data.position().threevec();
    /* Lattices can be used for calculation if 1-2 are fulfilled:
     * 1) Required lattices are not nullptr - possibly_use_lattice
     * 2) r is not out of required lattices */
    const bool use_lattice =
        possibly_use_lattice &&
        (pot.use_skyrme() ? UB_grad_lat->value_at(r, dUB_dr) : true) &&
        (pot.use_symmetry() ? UI3_grad_lat->value_at(r, dUI3_dr) : true);
    if (!pot.use_skyrme()) {
      dUB_dr = ThreeVector(0.0, 0.0, 0.0);
    }
    if (!pot.use_symmetry()) {
      dUI3_dr = ThreeVector(0.0, 0.0, 0.0);
    }
    if (!use_lattice) {
      const auto tmp = pot.potential_gradient(r, plist);
      dUB_dr = tmp.first;
      dUI3_dr = tmp.second;
    }
    const ThreeVector dU_dr =
        scale.first * dUB_dr +
        scale.second * data.type().isospin3_rel() * dUI3_dr;
    log.debug("Update momenta: dU/dr [GeV/fm] = ", dU_dr);
    data.set_4momentum(data.effective_mass(),
                       data.momentum().threevec() - dU_dr * dt);

    // calculate the time scale of the change in momentum
    const double dU_dr_abs = dU_dr.abs();
    if (dU_dr_abs < really_small) {
      continue;
    }
    const double time_scale = data.momentum().x0() / dU_dr_abs;
    if (time_scale < min_time_scale) {
      min_time_scale = time_scale;
    }
  }

  // warn if the time step is too big
  constexpr double safety_factor = 0.1;
  if (dt > safety_factor * min_time_scale) {
    log.warn() << "The time step size is too large for an accurate propagation "
               << "with potentials. Maximum safe value: "
               << safety_factor * min_time_scale << " fm/c.";
  }
}

}  // namespace smash
