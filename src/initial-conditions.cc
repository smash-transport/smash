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
#include "include/random.h"
#include "include/outputroutines.h"

/* Set default IC for the box */
void init_box(box *box) {
  int steps = 10000, update = 10;
  int64_t seed = 1;
  float A = 10.0, EPS = 0.01, temperature = 0.3, sigma = 10.0;

  box->set(A, EPS, seed, sigma, steps, temperature, update);
}

/* density_integrand - Maxwell-Boltzmann distribution */
static double inline density_integrand(double momentum, double temp,
  double mass) {
  return 4 * M_PI * momentum * momentum
    * exp(- sqrt(momentum * momentum + mass * mass) / temp);
}

/* sample_momenta - return thermal momenta */
static double sample_momenta(box *box, ParticleType pi) {
  double momentum_radial, momentum_average, momentum_min, momentum_max;
  double probability = 0, probability_max, probability_random = 1;

  /* massless particles peak would be at <E>=3T */
  momentum_average = sqrt((3 * box->temperature()) * (3 * box->temperature())
    - pi.mass() * pi.mass());
  momentum_min = pi.mass();
  momentum_max = 50.0 * box->temperature();
  /* double the massless peak value to be above maximum of the distribution */
  probability_max = 2 * density_integrand(momentum_average, box->temperature(),
    pi.mass());

  /* random momenta and random probability need to be below the distribution */
  while (probability_random > probability) {
    momentum_radial = (momentum_max - momentum_min) * drand48() + momentum_min;
    momentum_radial = sqrt(momentum_radial * momentum_radial
      - pi.mass() * pi.mass());
    probability = density_integrand(momentum_radial, box->temperature(),
      pi.mass());
    probability_random = probability_max * drand48();
  }

  return momentum_radial;
}

/* initial_conditions - sets partilce data for @particles */
ParticleData* initial_conditions(ParticleData *particles, int &number,
      box *box) {
  double x_pos, y_pos, z_pos, time_start, number_density;
  double phi, theta, momentum_radial;
  FourVector momentum_total(0, 0, 0, 0);
  ParticleType pi("pi", 0.13957);
  ParticleType pi0("pi0", 0.134977);

  /* initialize random seed */
  srand48(box->seed());

  /* XXX: move to proper startup */
  printd("Pi^± mass: %g [GeV]\n", pi.mass());
  printd("Pi^0 mass: %g [GeV]\n", pi0.mass());

  /* 
   * The particle number depends on distribution function
   * (assumes Bose-Einstein):
   * Volume m^2 T BesselK[2, m/T] / (2\pi^2)
   */
  number_density = pi.mass() * pi.mass() * box->temperature()
    * gsl_sf_bessel_Knu(2, pi.mass() / box->temperature())
    / 2 / M_PI / M_PI / hbarc / hbarc / hbarc;
  /* cast while reflecting probability of extra particle */
  number = box->a() * box->a() * box->a() * number_density;
  if (box->a() * box->a() * box->a() * number_density - number > drand48())
    number++;
  printf("IC number density %.6g [fm^-3]\n", number_density);
  printf("IC %d number of %s\n", number, pi.name().c_str());

  /* Set random IC:
   * particles at random position in the box with thermal momentum
   */
  particles = new ParticleData[number];
  for (int i = 0; i < number; i++) {
    particles[i].set_id(i);

    /* thermal momentum according Maxwell-Boltzmann distribution */
    momentum_radial = sample_momenta(box, pi);
    /* back to back pair creation with random momenta direction */
    if (unlikely(i == number - 1 && !(i % 2))) {
      /* poor last guy just sits around */
      particles[i].set_momentum(pi.mass(), 0, 0, 0);
    } else if (!(i % 2)) {
      phi =  2 * M_PI * drand48();
      theta = M_PI * drand48();
      printd("Particle %d radial momenta %g phi %g theta %g\n", i,
        momentum_radial, phi, theta);
      particles[i].set_momentum(pi.mass(),
        momentum_radial * cos(phi) * sin(theta),
        momentum_radial * sin(phi) * sin(theta),
        momentum_radial * cos(theta));
    } else {
      particles[i].set_momentum(pi.mass(),
        - particles[i - 1].momentum().x1(),
        - particles[i - 1].momentum().x2(),
        - particles[i - 1].momentum().x3());
    }
    momentum_total += particles[i].momentum();

    /* ramdom position in the box */
    time_start = 1.0;
    x_pos = drand48() * box->a();
    y_pos = drand48() * box->a();
    z_pos = drand48() * box->a();
    particles[i].set_position(time_start, z_pos, x_pos, y_pos);

    /* no collision yet to happen */
    particles[i].set_collision(0, 0);

    /* IC: debug checks */
    printd_momenta(particles[i]);
    printd_position(particles[i]);
  }
  /* allows to check energy conservation */
  printf("IC total energy: %g [GeV]\n", momentum_total.x0());
  box->set_energy_initial(momentum_total.x0());
  print_header();

  return particles;
}
