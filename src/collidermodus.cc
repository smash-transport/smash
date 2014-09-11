/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "include/collidermodus.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <list>

#include "include/configuration.h"
#include "include/experimentparameters.h"
#include "include/logging.h"
#include "include/outputroutines.h"
#include "include/particles.h"
#include "include/random.h"

namespace Smash {
ColliderModus::ColliderModus(Configuration modus_config,
                             const ExperimentParameters &)
    : projectile_(modus_config.take({"Collider", "PROJECTILE"})
                      .convert_for(projectile_)),
      target_(modus_config.take({"Collider", "TARGET"}).convert_for(target_)),
      sqrts_(modus_config.take({"Collider", "SQRTS"})) {
  if (sqrts_ < ParticleType::find(projectile_).mass()
             + ParticleType::find(target_).mass()) {
    throw ModusDefault::InvalidEnergy(
          "Error in input: sqrt(s) is smaller than masses:\n"
          + std::to_string(sqrts_) + " GeV < "
          + std::to_string(ParticleType::find(projectile_).mass()) + " GeV + "
          + std::to_string(ParticleType::find(projectile_).mass()) + " GeV.");
  }
}

/* console output on startup of box specific parameters */
std::ostream &operator<<(std::ostream &out, const ColliderModus &m) {
  return out << "-- Collider Modus:\n"
                "Projectile PDG ID: " << m.projectile_
             << "\nTarget PDG ID: " << m.target_ << "\nCenter-of-mass energy "
             << format(m.sqrts_, "GeV", 10, 3);
}

/* initial_conditions - sets particle data for @particles */
float ColliderModus::initial_conditions(Particles *particles,
                                       const ExperimentParameters &) {
  const auto &log = logger<LogArea::Collider>();

  /* Create "projectile" particle */
  particles->create(1, projectile_);
  /* Pointer to "projectile" data */
  ParticleData &data_projectile = particles->data(particles->id_max());
  float mass_projectile = data_projectile.pole_mass();
  log.debug() << "Projectile: PDG code " << data_projectile.pdgcode()
              << ", mass: " << mass_projectile;
  /* Create "target" particle */
  particles->create(1, target_);
  /* Pointer to "target" data */
  ParticleData &data_target = particles->data(particles->id_max());
  float mass_target = data_target.pole_mass();
  log.debug() << "Target: PDG code " << data_target.pdgcode()
              << ", mass: " << mass_target;
  /* Projectile energy in CMS */
  double cms_energy_projectile = (sqrts_ * sqrts_
                                  + mass_projectile * mass_projectile
                                  - mass_target * mass_target)
                                 / (2 * sqrts_);
  /* CMS momentum, same in magnitude for both */
  double cms_momentum = sqrt(cms_energy_projectile * cms_energy_projectile
                             - mass_projectile * mass_projectile);
  /* Sample impact parameter */
  double impact_parameter = Random::uniform(0.0, 5.0);
  // collider start is hard-coded for now.
  const float start_time = -1.0f;
  /* Set positions and momenta */
  data_projectile.set_4position(FourVector(start_time, impact_parameter,
                                           0., -1.));
  data_projectile.set_4momentum(mass_projectile, 0.0, 0.0, cms_momentum);
  data_target.set_4position(FourVector(start_time, 0., 0., 1.));
  data_target.set_4momentum(mass_target, 0.0, 0.0, -cms_momentum);
  return start_time;
}

}  // namespace Smash
