/*
 *
 *    Copyright (c) 2013-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/propagation.h"

#include "include/boxmodus.h"
#include "include/collidermodus.h"
#include "include/listmodus.h"
#include "include/logging.h"
#include "include/spheremodus.h"

namespace Smash {

/* Simple straight line propagation without potentials*/
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
    // "Frozen Fermi motion": Fermi momenta are only used for collisions,
    // but not for propagation. This is done to avoid nucleus flying apart
    // even if potentials are off. Initial nucleons before the first collision
    // are propagated only according to beam momentum.
    // Initial nucleons are distinguished by data.id() < the size of
    // beam_momentum, which is by default zero except for the collider modus
    // with the fermi motion == frozen.
    // todo(m. mayer): improve this condition (see comment #11 issue #4213)
    assert(data.id() >= 0);
    const bool avoid_fermi_motion =
      (static_cast<uint64_t>(data.id())
       < static_cast<uint64_t>(beam_momentum.size()))
      && (data.get_history().collisions_per_particle == 0);
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


void update_momenta(Particles *particles, double dt,
               const Potentials &pot,
               RectangularLattice<ThreeVector>* UB_grad_lat,
               RectangularLattice<ThreeVector>* UI3_grad_lat) {
  // Copy particles before propagation to calculate potentials from them
  const ParticleList plist = particles->copy_to_vector();

  const auto &log = logger<LogArea::Propagation>();
  bool possibly_use_lattice =
         (pot.use_skyrme() ? (UB_grad_lat != nullptr) : true) &&
         (pot.use_symmetry() ? (UI3_grad_lat != nullptr) : true);
  ThreeVector dUB_dr, dUI3_dr;
  double min_time_scale = std::numeric_limits<double>::infinity();

  for (ParticleData &data : *particles) {
    const ThreeVector r = data.position().threevec();
    /* Lattices can be used for calculation if 1-2 are fulfilled:
     * 1) Required lattices are not nullptr - possibly_use_lattice
     * 2) r is not out of required lattices
     */
    const bool use_lattice = possibly_use_lattice &&
              (pot.use_skyrme() ? UB_grad_lat->value_at(r, dUB_dr) : true) &&
              (pot.use_symmetry() ? UI3_grad_lat->value_at(r, dUI3_dr) : true);
    if (!pot.use_skyrme()) {
      dUB_dr = ThreeVector(0.0, 0.0, 0.0);
    }
    if (!pot.use_symmetry()) {
      dUI3_dr = ThreeVector(0.0, 0.0, 0.0);
    }
    // Compute potential gradient from lattice if possible
    const ThreeVector dU_dr = use_lattice ? (dUB_dr + dUI3_dr):
                              pot.potential_gradient(r, plist, data.type());
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

}  // namespace Smash
