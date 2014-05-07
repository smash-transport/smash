/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <cmath>

#include "include/constants.h"
#include "include/configuration.h"
#include "include/crosssections.h"
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
      number_of_particles_(modus_config.take({"Sphere","NUMBEROFPARTICLES"})) {
}

/* print_startup - console output on startup of sphere specific parameters */
void SphereModus::print_startup() {
  printf("Radius of the sphere: %g [fm]\n", radius_);
  printf("Total number of particles in sphere: %i \n", number_of_particles_);
}


/* initial_conditions - sets particle data for @particles */
void SphereModus::initial_conditions(Particles *particles,
                                     const ExperimentParameters &parameters){
  FourVector momentum_total(0, 0, 0, 0);
  /* loop over all the particle types creating each particles */
//  for (auto i = particles->types_cbegin(); i != particles->types_cend(); ++i) {
    /* Particles with width > 0 (resonances) do not exist in the beginning */
//    if (data.width() > 0.0) continue;
//    printd("%s mass: %g [GeV]\n", data.name().c_str(), data.mass());
    /* Maxwell-Boltzmann statistics with temperature 0.3 GeV */
//    double number_density = number_density_maxwellboltzmann(data.mass(), 0.3);
//    printf("IC number density %.6g [fm^-3]\n", number_density);
//    double real_number = 4.0 / 3.0 * M_PI * radius_ * radius_ * radius_ *
//                         number_density * testparticles;
//    size_t int_number = static_cast<size_t>(real_number);
//    if (real_number - int_number > Random::canonical())
//      ++number;
    /* create bunch of particles */
//    printf("IC creating %zu particles\n", number);
//    particles->create(number, data.pdgcode());
//    number_total += number;
//  }
//  printf("IC contains %zu particles\n", number_total);
  /* now set position and momentum of the particles */
//  double momentum_radial;
//  Angles phitheta = Angles();
//  auto uniform_radius = Random::make_uniform_distribution(-radius_, +radius_);
//  for (ParticleData &data : particles->data()) {
//    if (unlikely(data.id() == particles->id_max() && !(data.id() % 2))) {
      /* poor last guy just sits around */
//      data.set_momentum(particles->particle_type(data.pdgcode()).mass(), 0, 0, 0);
//    } else if (!(data.id() % 2)) {
      /* thermal momentum according Maxwell-Boltzmann distribution */
//      momentum_radial = sample_momenta(0.3, particles->particle_type(data.pdgcode()).mass());
//      phitheta = Angles().distribute_isotropically();
//      printd("Particle %d radial momenta %g phi %g cos_theta %g\n", data.id(),
//             momentum_radial, phitheta.phi(), phitheta.costheta());
//      data.set_momentum(
//          particles->particle_type(data.pdgcode()).mass(), momentum_radial * phitheta.x(),
//          momentum_radial * phitheta.y(), momentum_radial * phitheta.z());
//    } else {
//      data.set_momentum(particles->particle_type(data.pdgcode()).mass(),
//                             -particles->data(data.id() - 1).momentum().x1(),
//                             -particles->data(data.id() - 1).momentum().x2(),
//                             -particles->data(data.id() - 1).momentum().x3());
//    }
//    momentum_total += data.momentum();
//    double x, y, z;
    /* ramdom position in a sphere
     * box length here has the meaning of the sphere radius
     */
//    x = uniform_radius();
//    y = uniform_radius();
//    z = uniform_radius();
    /* sampling points inside of the sphere, rejected if outside */
//    while (sqrt(x * x + y * y + z * z) > radius_) {
 //     x = uniform_radius();
 //     y = uniform_radius();
 //     z = uniform_radius();
  //  }
 //   data.set_position(time_start, x, y, z);
    /* IC: debug checks */
 //   printd_momenta(data);
 //   printd_position(data);
 // }
 // printf("IC total energy: %g [GeV]\n", momentum_total.x0());
}

    
} // namespace Smash
