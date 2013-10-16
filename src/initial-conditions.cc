/*
 *
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */
#include "include/initial-conditions.h"

#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>

#include "include/Box.h"
#include "include/constants.h"
#include "include/distributions.h"
#include "include/FourVector.h"
#include "include/Laboratory.h"
#include "include/macros.h"
#include "include/Particles.h"
#include "include/ParticleData.h"
#include "include/ParticleType.h"
#include "include/outputroutines.h"
#include "include/Sphere.h"

/* initial_conditions - sets particle type */
void initial_particles(Particles *particles) {
  /* XXX: use nosql table for particle type values */
  ParticleType piplus("pi+", 0.13957f, -1.0, 211, 1, 1, 0);
  particles->add_type(piplus, 211);
  ParticleType piminus("pi-", 0.13957f, -1.0, -211, 1, -1, 0);
  particles->add_type(piminus, -211);
  ParticleType pi0("pi0", 0.134977f, -1.0, 111, 1, 0, 0);
  particles->add_type(pi0, 111);
}

/* initial_conditions - sets particle data for @particles */
void initial_conditions(Particles *particles, Box *box) {
  double phi, cos_theta, sin_theta, momentum_radial, number_density_total = 0;
  FourVector momentum_total(0, 0, 0, 0);
  size_t number_total = 0, number = 0;

  /* loop over all the particle types */
  for (std::map<int, ParticleType>::const_iterator
       i = particles->types_cbegin(); i != particles->types_cend(); ++i) {
    /* Particles with width > 0 (resonances) do not exist in the beginning */
    if (i->second.width() > 0.0)
      continue;
    printd("%s mass: %g [GeV]\n", i->second.name().c_str(), i->second.mass());

    /* bose einstein distribution funtion */
    double number_density = number_density_bose(i->second.mass(),
      box->temperature());

    /* cast while reflecting probability of extra particle */
    number = box->length() * box->length() * box->length() * number_density
        * box->testparticles();
    if (box->length() * box->length() * box->length() * number_density
        - number > drand48())
      number++;

    printf("IC number density %.6g [fm^-3]\n", number_density);
    printf("IC %zu number of %s\n", number, i->second.name().c_str());
    number_density_total += number_density;
    /* create bunch of particles */
    printf("IC creating %zu particles\n", number);
    particles->create(number, i->second.pdgcode());
    number_total += number;
  }
  printf("IC total number density %.6g [fm^-3]\n", number_density_total);
  printf("IC contains %zu particles\n", number_total);

  /* Set paricles IC: */
  for (std::map<int, ParticleData>::iterator i = particles->begin();
       i != particles->end(); ++i) {
    double x, y, z, time_start;
    /* back to back pair creation with random momenta direction */
    if (unlikely(i->first == particles->id_max() && !(i->first % 2))) {
      /* poor last guy just sits around */
      i->second.set_momentum(particles->type(i->first).mass(), 0, 0, 0);
    } else if (!(i->first % 2)) {
      if (box->initial_condition() != 2) {
        /* thermal momentum according Maxwell-Boltzmann distribution */
        momentum_radial = sample_momenta(box->temperature(),
                                         particles->type(i->first).mass());
      } else {
        /* IC == 2 initial thermal momentum is the average 3T */
        momentum_radial = 3.0 * box->temperature();
      }
      /* phi in the range from [0, 2 * pi) */
      phi = 2.0 * M_PI * drand48();
      /* cos(theta) in the range from [-1.0, 1.0) */
      cos_theta = -1.0 + 2.0 * drand48();
      sin_theta = sqrt(1.0 - cos_theta * cos_theta);
      printd("Particle %zu radial momenta %g phi %g cos_theta %g\n", i->first,
          momentum_radial, phi, cos_theta);
      i->second.set_momentum(particles->type(i->first).mass(),
          momentum_radial * cos(phi) * sin_theta,
          momentum_radial * sin(phi) * sin_theta,
          momentum_radial * cos_theta);
    } else {
      i->second.set_momentum(particles->type(i->first).mass(),
          - particles->data(i->first - 1).momentum().x1(),
          - particles->data(i->first - 1).momentum().x2(),
          - particles->data(i->first - 1).momentum().x3());
    }
    momentum_total += i->second.momentum();

    time_start = 1.0;
    /* ramdom position in a quadratic box */
    x = drand48() * box->length();
    y = drand48() * box->length();
    z = drand48() * box->length();
    i->second.set_position(time_start, x, y, z);

    /* IC: debug checks */
    printd_momenta(i->second);
    printd_position(i->second);
  }

  /* Display on startup if pseudo grid is used */
  number = number_total;
  int const grid_number = round(box->length()
              / sqrt(box->cross_section() * fm2_mb * M_1_PI) * 0.5);
  /* pseudo grid not used for 3^3 or extremely small particle numbers */
  if (grid_number >= 4 && number > 10)
    printf("Simulation with pseudo grid: %d^3\n", grid_number);
  else
    printf("W: Not using pseudo grid: %d^3\n", grid_number);

  /* allows to check energy conservation */
  printf("IC total energy: %g [GeV]\n", momentum_total.x0());
  box->set_energy_initial(momentum_total.x0());
  box->set_number_density_inital(number_density_total);
}

/* initial_conditions - sets particle data for @particles */
void initial_conditions(Particles *particles, Sphere *ball) {
  size_t number_total = 0;
  double time_start = 1.0;
  FourVector momentum_total(0, 0, 0, 0);

  /* loop over all the particle types creating each particles */
  for (std::map<int, ParticleType>::const_iterator
       i = particles->types_cbegin(); i != particles->types_cend(); ++i) {
    /* Particles with width > 0 (resonances) do not exist in the beginning */
    if (i->second.width() > 0.0)
      continue;
    printd("%s mass: %g [GeV]\n", i->second.name().c_str(),
           particles->type(i->first).mass());

    /* bose einstein distribution funtion with temperature 0.3 GeV */
    double number_density = number_density_bose(
                            particles->type(i->first).mass(), 0.3);
    printf("IC number density %.6g [fm^-3]\n", number_density);

    /* cast while reflecting probability of extra particle */
    size_t number = 4.0 / 3.0 * M_PI * ball->radius() * ball->radius()
      * ball->radius() * number_density * ball->testparticles();
    if (4.0 / 3.0 * M_PI * ball->radius() * ball->radius() * ball->radius()
        * number_density - number > drand48())
        number++;

    /* create bunch of particles */
    printf("IC creating %zu particles\n", number);
    particles->create(number, i->second.pdgcode());
    number_total += number;
  }
  printf("IC contains %zu particles\n", number_total);

  /* now set position and momentum of the particles */
  double momentum_radial, phi, cos_theta, sin_theta;
  for (std::map<int, ParticleData>::iterator i = particles->begin();
       i != particles->end(); ++i) {
    if (unlikely(i->first == particles->id_max() && !(i->first % 2))) {
      /* poor last guy just sits around */
      i->second.set_momentum(particles->type(i->first).mass(), 0, 0, 0);
    } else if (!(i->first % 2)) {
      /* thermal momentum according Maxwell-Boltzmann distribution */
      momentum_radial = sample_momenta(0.3, particles->type(i->first).mass());
      /* phi in the range from [0, 2 * pi) */
      phi = 2.0 * M_PI * drand48();
      /* cos(theta) in the range from [-1.0, 1.0) */
      cos_theta = -1.0 + 2.0 * drand48();
      sin_theta = sqrt(1.0 - cos_theta * cos_theta);
      printd("Particle %d radial momenta %g phi %g cos_theta %g\n", i->first,
        momentum_radial, phi, cos_theta);
      i->second.set_momentum(particles->type(i->first).mass(),
        momentum_radial * cos(phi) * sin_theta,
        momentum_radial * sin(phi) * sin_theta,
        momentum_radial * cos_theta);
    } else {
      i->second.set_momentum(particles->type(i->first).mass(),
        - particles->data(i->first - 1).momentum().x1(),
        - particles->data(i->first - 1).momentum().x2(),
        - particles->data(i->first - 1).momentum().x3());
    }
    momentum_total += i->second.momentum();

    double x, y, z;
    /* ramdom position in a sphere
     * box length here has the meaning of the sphere radius
     */
    x = -ball->radius() + 2.0 * drand48() * ball->radius();
    y = -ball->radius() + 2.0 * drand48() * ball->radius();
    z = -ball->radius() + 2.0 * drand48() * ball->radius();
    /* sampling points inside of the sphere, rejected if outside */
    while (sqrt(x * x + y * y + z * z) > ball->radius()) {
      x = -ball->radius() + 2.0 * drand48() * ball->radius();
      y = -ball->radius() + 2.0 * drand48() * ball->radius();
      z = -ball->radius() + 2.0 * drand48() * ball->radius();
    }
    i->second.set_position(time_start, x, y, z);
    /* IC: debug checks */
    printd_momenta(i->second);
    printd_position(i->second);
  }
  printf("IC total energy: %g [GeV]\n", momentum_total.x0());
}
