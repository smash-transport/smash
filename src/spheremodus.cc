/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <cmath>

#include "include/angles.h"
#include "include/constants.h"
#include "include/configuration.h"
#include "include/crosssections.h"
#include "include/distributions.h"
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
      sphere_temperature_(modus_config.take({"Sphere","SPHERETEMPERATURE"})){
}

/* print_startup - console output on startup of sphere specific parameters */
void SphereModus::print_startup() {
  printf("Radius of the sphere: %g [fm]\n", radius_);
  printf("Total number of particles in sphere: %i \n", number_of_particles_);
    printf("Temperature for momentum sampling: %f \n", sphere_temperature_);
}


/* initial_conditions - sets particle data for @particles */
void SphereModus::initial_conditions(Particles *particles,
                                     const ExperimentParameters &parameters){
/* count number of stable types */
  int number_of_stable_types = 0;
    /* loop over all the particle types */
  for (const ParticleType &type : particles->types()) {
        /* Particles with width > 0 (resonances) do not exist in the beginning */
    if (type.width() > 0.0) {
            continue;
    }
        number_of_stable_types = number_of_stable_types + 1;
  }

/* just produce equally many particles per type */
  int number_of_particles_per_type;
  number_of_particles_per_type = number_of_particles_/number_of_stable_types;
  for (const ParticleType &type : particles->types()) {
        /* Particles with width > 0 (resonances) do not exist in the beginning */
    if (type.width() > 0.0) {
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
    /* back to back pair creation with random momenta direction */
    if (unlikely(data.id() == particles->id_max() && !(data.id() % 2))) {
    /* poor last guy just sits around */
      data.set_momentum(particles->particle_type(data.pdgcode()).mass(), 0, 0, 0);
    } else if (!(data.id() % 2)) {
    /* thermal momentum according Maxwell-Boltzmann distribution */
      double momentum_radial;
      momentum_radial = sample_momenta(this->sphere_temperature_,
                                                 particles->particle_type(data.pdgcode()).mass());
      phitheta.distribute_isotropically();
      data.set_momentum(particles->particle_type(data.pdgcode()).mass(), momentum_radial * phitheta.x(),momentum_radial * phitheta.y(), momentum_radial * phitheta.z());
    }
    time_begin = 1.0;
    /* random position in a quadratic box */
    x = uniform_radius();
    y = uniform_radius();
    z = uniform_radius();
    data.set_position(time_begin, x, y, z);
  }

}
    
} // namespace Smash
