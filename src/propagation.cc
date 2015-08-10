/*
 *
 *    Copyright (c) 2013-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/propagation.h"

#include "include/logging.h"
#include "include/boxmodus.h"
#include "include/collidermodus.h"
#include "include/listmodus.h"
#include "include/spheremodus.h"

namespace Smash {

/* Simple straight line propagation without potentials*/
void propagate_straight_line(Particles *particles,
                             const ExperimentParameters &parameters) {
  const auto &log = logger<LogArea::Propagation>();
  const double dt = parameters.timestep_duration();
  for (ParticleData &data : *particles) {
    FourVector distance = FourVector(0.0, data.velocity() * dt);
    log.debug("Particle ", data, " motion: ", distance);
    FourVector position = data.position() + distance;
    position.set_x0(parameters.new_particle_time());
    data.set_4position(position);
    /* Test, if particle is formed and reset cross_section_scaling_factor
       TODO: Is there a way to only do that for particles that were unformed 
       in the previous timestep? */
    if (data.formation_time() < data.position().x0()) {
      data.set_cross_section_scaling_factor(1.0);
    }
  }
}

void propagate(Particles *particles, const ExperimentParameters &parameters,
               const Potentials &pot,
               RectangularLattice<ThreeVector>* UB_grad_lat,
               RectangularLattice<ThreeVector>* UI3_grad_lat) {
  // Copy particles before propagation to calculate potentials from them
  const ParticleList plist = particles->copy_to_vector();
  const double dt = parameters.timestep_duration();
  const auto &log = logger<LogArea::Propagation>();
  bool possibly_use_lattice =
         (pot.use_skyrme() ? (UB_grad_lat != nullptr) : true) &&
         (pot.use_symmetry() ? (UI3_grad_lat != nullptr) : true);
  bool use_lattice;
  ThreeVector dUB_dr, dUI3_dr;
  float min_time_scale = std::numeric_limits<float>::infinity();

  for (ParticleData &data : *particles) {
    ThreeVector r = data.position().threevec();
    /* Lattices can be used for calculation if 1-2 are fulfilled:
     * 1) Required lattices are not nullptr - possibly_use_lattice
     * 2) r is not out of required lattices
     */
    use_lattice = possibly_use_lattice &&
              (pot.use_skyrme() ? UB_grad_lat->value_at(r, dUB_dr) : true) &&
              (pot.use_symmetry() ? UI3_grad_lat->value_at(r, dUI3_dr) : true);
    if (!pot.use_skyrme()) {
      dUB_dr = ThreeVector(0.0, 0.0, 0.0);
    }
    if (!pot.use_symmetry()) {
      dUI3_dr = ThreeVector(0.0, 0.0, 0.0);
    }
    // Compute potential gradient from lattice if possible
    ThreeVector dU_dr = use_lattice ? (dUB_dr + dUI3_dr):
                        pot.potential_gradient(r, plist, data.pdgcode());
    log.debug("Propagate: dU/dr = ", dU_dr);
    ThreeVector v = data.velocity();
    // predictor step assuming momentum-indep. potential, dU/dp = 0
    // then for momentum predictor = corrector
    data.set_4momentum(data.effective_mass(),
                       data.momentum().threevec() - dU_dr * dt);
    ThreeVector v_pred = data.velocity();
    // corrector step
    FourVector distance = FourVector(0.0, (v + v_pred) * (0.5 * dt));
    log.debug("Particle ", data, " motion: ", distance);
    FourVector position = data.position() + distance;
    position.set_x0(parameters.new_particle_time());
    data.set_4position(position);

    // calculate the time scale of the change in momentum
    const double dU_dr_abs = dU_dr.abs();
    if (dU_dr_abs < really_small) {
      continue;
    }
    const float time_scale = data.momentum().x0() / dU_dr_abs;
    if (time_scale < min_time_scale) {
      min_time_scale = time_scale;
    }
  }

  // warn if the time step is too big
  constexpr float safety_factor = 0.1f;
  if (parameters.timestep_duration() > safety_factor * min_time_scale) {
    log.warn() << "The time step size is too large for an accurate propagation "
               << "with potentials.";
  }
}

}  // namespace Smash
