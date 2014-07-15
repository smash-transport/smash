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
#include "include/particledata.h"
#include "include/particles.h"
#include "include/random.h"
#include "include/spheremodus.h"


namespace Smash {

SphereModus::SphereModus(Configuration modus_config,
                         const ExperimentParameters &)
    : radius_(modus_config.take({"Sphere", "RADIUS"})),
      number_of_particles_(modus_config.take({"Sphere","NUMBEROFPARTICLES"})),
      sphere_temperature_(modus_config.take({"Sphere","SPHERETEMPERATURE"})),
      start_time_       (modus_config.take({"Sphere", "START_TIME"})){
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
                                     const ExperimentParameters &parameters){
/* count number of stable types */
  int number_of_stable_types = 0;
    /* loop over all the particle types */
  for (const ParticleType &type : ParticleType::list_all()) {
        /* Particles with width > 0 (resonances) do not exist in the beginning */
    if (!type.is_stable()) {
            continue;
    }
        number_of_stable_types = number_of_stable_types + 1;
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
    
auto uniform_radius = Random::make_uniform_distribution(0.0,
                                        static_cast<double>(this->radius_));
    /* loop over particle data to fill in momentum and position information */
  for (ParticleData &data : particles->data()) {
    double x, y, z, time_begin;
    Angles phitheta;
    /* thermal momentum according Maxwell-Boltzmann distribution */
      double momentum_radial;
      momentum_radial = sample_momenta(this->sphere_temperature_,
                                                 data.mass());
      phitheta.distribute_isotropically();
       data.set_momentum(data.mass(), phitheta.threevec() * momentum_radial);

      ThreeVector pos{uniform_radius(), uniform_radius(), uniform_radius()};
      data.set_position(FourVector(start_time_, pos));
      return start_time_;
  }

}
    
} // namespace Smash
