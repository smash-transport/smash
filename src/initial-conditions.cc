/*
 *
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */
#include "include/initial-conditions.h"

#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>

#include "include/box.h"
#include "include/constants.h"
#include "include/ParticleData.h"
#include "include/ParticleType.h"
#include "include/random.h"
#include "include/outputroutines.h"

/* Set default IC for the box */
box init_box(box box) {
  int steps = 10000, update = 10;
  float A = 10.0, EPS = 0.01, temperature = 0.3, sigma = 10.0;

  box.set(steps, update, A, EPS, temperature, sigma);
  return box;
}

/* initial_conditions - sets partilce data for @particles */
ParticleData* initial_conditions(ParticleData *particles, int &number,
      box box) {
  double x_pos, y_pos, z_pos, time_start, number_density;
  double phi, theta, momentum_radial;
  ParticleType pi("pi", 0.13957);
  ParticleType pi0("pi0", 0.134977);

  /* XXX: move to proper startup */
  printd("Pi^± mass: %g [GeV]\n", pi.mass());
  printd("Pi^0 mass: %g [GeV]\n", pi0.mass());

  /* 
   * The particle number depends on distribution function
   * (assumes Bose-Einstein):
   * Volume m^2 T BesselK[2, m/T] / (2\pi^2)
   */
  number_density = pi.mass() * pi.mass() * box.temperature()
    * gsl_sf_bessel_Knu(2, pi.mass() / box.temperature())
    / 2 / M_PI / M_PI / hbarc / hbarc / hbarc;
  /* cast while reflecting probability of extra particle */
  number = box.a() * box.a() * box.a() * number_density;
  srand48(time(NULL));
  if (box.a() * box.a() * box.a() * number_density - number > drand48())
    number++;
  printf("IC number density %.6g [fm^-3]\n", number_density);
  printf("IC %d number of %s\n", number, pi.name().c_str());
  print_header();

  /* Set random IC:
   * particles at random position in the box with thermal momentum
   */
  particles = new ParticleData[number];
  for (int i = 0; i < number; i++) {
    particles[i].set_id(i);

    /* XXX: replace with full distribution */
    momentum_radial = box_muller(sqrt((3 * box.temperature())
      * (3 * box.temperature()) - pi.mass() * pi.mass()), 0.01);
    /* back to back pair creation with random momenta direction */
    if (unlikely(i == number - 1 && !(i % 2))) {
      /* poor last guy just sits around */
      particles[i].set_momentum(pi.mass(), 0, 0, 0);
    } else if (!(i % 2)) {
      phi =  2 * M_PI * rand_r(&seedp) / RAND_MAX;
      theta = M_PI * rand_r(&seedp) / RAND_MAX;
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

    time_start = 1.0;
    x_pos = 1.0 * rand_r(&seedp) / RAND_MAX * box.a();
    y_pos = 1.0 * rand_r(&seedp) / RAND_MAX * box.a();
    z_pos = 1.0 * rand_r(&seedp) / RAND_MAX * box.a();
    particles[i].set_position(time_start, z_pos, x_pos, y_pos);

    printd_momenta(particles[i]);
    printd_position(particles[i]);
  }

  return particles;
}
