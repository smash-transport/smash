/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <map>
#include <utility>
#include <vector>

#include "include/algorithms.h"
#include "include/angles.h"
#include "include/constants.h"
#include "include/configuration.h"
#include "include/distributions.h"
#include "include/experimentparameters.h"
#include "include/fourvector.h"
#include "include/logging.h"
#include "include/macros.h"
#include "include/outputroutines.h"
#include "include/particles.h"
#include "include/random.h"
#include "include/spheremodus.h"
#include "include/threevector.h"

namespace Smash {

SphereModus::SphereModus(Configuration modus_config,
                         const ExperimentParameters &)
    : radius_(modus_config.take({"Sphere", "RADIUS"})),
      number_of_particles_(modus_config.take({"Sphere", "NUMBEROFPARTICLES"})),
      sphere_temperature_(modus_config.take({"Sphere", "SPHERETEMPERATURE"})),
      start_time_(modus_config.take({"Sphere", "START_TIME"})) {
}

/* console output on startup of sphere specific parameters */
std::ostream &operator<<(std::ostream &out, const SphereModus &m) {
  return out << "-- Sphere Modus:\n"
                "Radius of the sphere: " << m.radius_ << " [fm]"
             << "\nTotal number of particles in sphere: "
             << m.number_of_particles_
             << "\nTemperature for momentum sampling: " << m.sphere_temperature_
             << "\nStarting time for Sphere calculation: " << m.start_time_
             << '\n';
}


/* initial_conditions - sets particle data for @particles */
float SphereModus::initial_conditions(Particles *particles,
  const ExperimentParameters& /*parameters*/) {
  const auto &log = logger<LogArea::Sphere>(); 
  FourVector momentum_total(0, 0, 0, 0);	  
  /* count number of stable types */
  int number_of_stable_types = 0;
  /* loop over all the particle types */
  for (const ParticleType &type : ParticleType::list_all()) {
    /* Particles with width > 0 (resonances) do not exist in the beginning */
    if (!type.is_stable()) {
      continue;
    }
    number_of_stable_types = number_of_stable_types + 1;
    printd("%s mass: %g [GeV]\n", type.name().c_str(), type.mass());
  }

  /* just produce equally many particles per type */
  int number_of_particles_per_type;
  number_of_particles_per_type = number_of_particles_/number_of_stable_types;
  for (const ParticleType &type : ParticleType::list_all()) {
  /* Particles with width > 0 (resonances) do not exist in the beginning */
    if (!type.is_stable()) {
      continue;
    }
  particles->create(number_of_particles_per_type, type.pdgcode());
  }
  /* loop over particle data to fill in momentum and position information */
  for (ParticleData &data : particles->data()) {
    Angles phitheta;
    /* thermal momentum according Maxwell-Boltzmann distribution */
    double momentum_radial;
    momentum_radial = sample_momenta(this->sphere_temperature_,
                                     data.pole_mass());
    phitheta.distribute_isotropically();
    log.debug() << data << ", radial mom:" << field << momentum_radial << ", "
                << phitheta;
    data.set_4momentum(data.pole_mass(), phitheta.threevec() * momentum_radial);
    momentum_total += data.momentum();
    /* uniform sampling in a sphere with radius r */
    double position_radial;
    position_radial = cbrt(Random::canonical()) * radius_;
    Angles pos_phitheta;
    pos_phitheta.distribute_isotropically();
    data.set_4position(FourVector(start_time_,
                                  pos_phitheta.threevec() * position_radial));
 }
  /* Make total 3-momentum 0 */
  for (ParticleData &data : particles->data()) {
    data.set_4momentum(data.pole_mass(), data.momentum().threevec() -
                       momentum_total.threevec()/particles->size());
  }

  /* Recalculate total momentum */
  momentum_total = FourVector(0, 0, 0, 0);
  for (ParticleData &data : particles->data()) {
    momentum_total += data.momentum();
    /* IC: debug checks */
    log.debug() << data;
  }
  /* allows to check energy conservation */
  log.info() << "Sphere initial total 4-momentum [GeV]: "
             << momentum_total;  
  return start_time_;
}
}  // namespace Smash
