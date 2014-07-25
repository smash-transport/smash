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

#include "include/angles.h"
#include "include/constants.h"
#include "include/configuration.h"
#include "include/crosssections.h"
#include "include/distributions.h"
#include "include/experimentparameters.h"
#include "include/fourvector.h"
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

/* print_startup - console output on startup of sphere specific parameters */
void SphereModus::print_startup() {
  printf("Radius of the sphere: %g [fm]\n", radius_);
  printf("Total number of particles in sphere: %i \n", number_of_particles_);
  printf("Temperature for momentum sampling: %f \n", sphere_temperature_);
  printf("Starting time for Sphere calculation: %f \n", start_time_);
}


/* initial_conditions - sets particle data for @particles */
float SphereModus::initial_conditions(Particles *particles,
  const ExperimentParameters& /*parameters*/) {
  /* count number of stable types */
  int number_of_stable_types = 0;
  /* loop over all the particle types */
  for (const ParticleType &type : ParticleType::list_all()) {
    /* Particles with width > 0 (resonances) do not exist in the beginning */
    if (!type.is_stable()) {
      continue;
    }
    number_of_stable_types = number_of_stable_types + 1;
    printd("%s mass: %g [GeV]\n", type.name().c_str(), type.pole_mass());
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
    printd("Particle %d radial momenta %g phi %g cos_theta %g\n", data.id(),
           momentum_radial, phitheta.phi(), phitheta.costheta());
    data.set_momentum(data.pole_mass(), phitheta.threevec() * momentum_radial);
    /* uniform sampling in a sphere with radius r */
    double position_radial;
    position_radial = cbrt(Random::canonical()) * radius_;
    Angles pos_phitheta;
    pos_phitheta.distribute_isotropically();
    data.set_position(FourVector(start_time_, pos_phitheta.threevec()
                                 * position_radial));
    /* IC Debug checks */
    printd_position(data);
    printd_momenta(data);
  }
  return start_time_;
}
}  // namespace Smash
