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

#include <stdio.h>
#include <stdint.h>
#include <gsl/gsl_sf_bessel.h>

#include "include/box.h"
#include "include/constants.h"
#include "include/ParticleData.h"
#include "include/ParticleType.h"
#include "include/outputroutines.h"

/* density_integrand - Maxwell-Boltzmann distribution */
static double inline density_integrand(double momentum, double temp,
  double mass) {
  return 4 * M_PI * momentum * momentum
    * exp(- sqrt(momentum * momentum + mass * mass) / temp);
}

/* sample_momenta - return thermal momenta */
static double sample_momenta(box *box, ParticleType type) {
  double momentum_radial, momentum_average, momentum_min, momentum_max;
  double probability = 0, probability_max, probability_random = 1;

  /* massless particles peak would be at <E>=3T */
  momentum_average = sqrt((3 * box->temperature()) * (3 * box->temperature())
    - type.mass() * type.mass());
  momentum_min = type.mass();
  momentum_max = 50.0 * box->temperature();
  /* double the massless peak value to be above maximum of the distribution */
  probability_max = 2 * density_integrand(momentum_average, box->temperature(),
    type.mass());

  /* sample by rejection method: (see numerical recipes for more efficient)
   * random momenta and random probability need to be below the distribution
   */
  while (probability_random > probability) {
    momentum_radial = (momentum_max - momentum_min) * drand48() + momentum_min;
    momentum_radial = sqrt(momentum_radial * momentum_radial
      - type.mass() * type.mass());
    probability = density_integrand(momentum_radial, box->temperature(),
      type.mass());
    probability_random = probability_max * drand48();
  }

  return momentum_radial;
}

/* initial_conditions - sets particle type */
ParticleType* initial_particles(ParticleType *type) {
  /* XXX: use nosql table for particle type values */
  type = new ParticleType[3];
  type[0].set("pi+", 0.13957, 211);
  type[1].set("pi-", 0.13957, -211);
  type[2].set("pi0", 0.134977, 111);

  return type;
}


/* initial_conditions - sets particle data for @particles */
ParticleData* initial_conditions(ParticleData *particles,
  ParticleType *type, int &number, box *box) {
  double phi, theta, momentum_radial, number_density_total = 0;
  FourVector momentum_total(0, 0, 0, 0);
  int number_total = 0;

  /* initialize random seed */
  srand48(box->seed());

  /* loop over all the particle types */
  for (int i = 0; i < 3; i++) {
    printd("%s mass: %g [GeV]\n", type[i].name().c_str(), type[i].mass());
    /* 
     * The particle number depends on distribution function
     * (assumes Bose-Einstein):
     * Volume m^2 T BesselK[2, m/T] / (2\pi^2)
     */
    double number_density = type[i].mass() * type[i].mass() * box->temperature()
      * gsl_sf_bessel_Knu(2, type[i].mass() / box->temperature())
      * 0.5 * M_1_PI * M_1_PI / hbarc / hbarc / hbarc;
    /* cast while reflecting probability of extra particle */
    number = box->a() * box->a() * box->a() * number_density
      * box->testparticle();
    if (box->a() * box->a() * box->a() * number_density - number > drand48())
      number++;
    printf("IC number density %.6g [fm^-3]\n", number_density);
    printf("IC %d number of %s\n", number, type[i].name().c_str());
    number_density_total += number_density;

    /* Set random IC:
     * particles at random position in the box with thermal momentum
     */
    if (!particles) {
      particles = new ParticleData[number];
    } else {
      /* XXX: hack realloc, maybe better to use resize from vector class */
      ParticleData *particles_tmp = new ParticleData[number_total + number];
      for (int id = 0; id < number_total; id++)
        particles_tmp[id] = particles[id];
      delete[] particles;
      particles = particles_tmp;
    }
    for (int id = number_total; id < number_total + number; id++) {
      double x_pos, y_pos, z_pos, time_start;
      /* set id and particle type */
      particles[id].set_id(id);

      /* back to back pair creation with random momenta direction */
      if (unlikely(id == number + number_total - 1 && !(id % 2) && i == 2)) {
        /* poor last guy just sits around */
        particles[id].set_momentum(type[i].mass(), 0, 0, 0);
      } else if (!(id % 2)) {
        /* thermal momentum according Maxwell-Boltzmann distribution */
        momentum_radial = sample_momenta(box, type[i]);
        phi =  2 * M_PI * drand48();
        theta = M_PI * drand48();
        printd("Particle %d radial momenta %g phi %g theta %g\n", id,
          momentum_radial, phi, theta);
        particles[id].set_momentum(type[i].mass(),
          momentum_radial * cos(phi) * sin(theta),
          momentum_radial * sin(phi) * sin(theta),
          momentum_radial * cos(theta));
      } else {
        particles[id].set_momentum(type[i].mass(),
          - particles[id - 1].momentum().x1(),
          - particles[id - 1].momentum().x2(),
          - particles[id - 1].momentum().x3());
      }
      momentum_total += particles[id].momentum();

      /* ramdom position in the box */
      time_start = 1.0;
      x_pos = drand48() * box->a();
      y_pos = drand48() * box->a();
      z_pos = drand48() * box->a();
      particles[id].set_position(time_start, z_pos, x_pos, y_pos);

      /* no collision yet hence zero time and unexisting id */
      particles[id].set_collision(0, -1);

      /* IC: debug checks */
      printd_momenta(particles[id]);
      printd_position(particles[id]);
    }
    number_total += number;
  }
  printf("IC total number density %.6g [fm^-3]\n", number_density_total);
  printf("IC contains %i particles\n", number_total);
  /* loop over all particles */
  number = number_total;
  /* reducing cross section according to number of test particle */
  if (box->testparticle() > 1) {
    printf("IC test particle: %i\n", box->testparticle());
    box->set_cross_section(box->cross_section() / box->testparticle());
    printf("Elastic cross section: %g [mb]\n", box->cross_section());
  }

  /* allows to check energy conservation */
  printf("IC total energy: %g [GeV]\n", momentum_total.x0());
  box->set_energy_initial(momentum_total.x0());
  box->set_number_density_inital(number_density_total);
  print_header();

  return particles;
}
