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

#include "include/box.h"
#include "include/ParticleData.h"
#include "include/ParticleType.h"
#include "include/random.h"
#include "include/outputroutines.h"

/* Set default IC for the box */
box init_box(box box) {
  int steps = 10000;
  float A = 5.0, EPS = 0.05, temperature = 0.1;

  box.set(steps, A, EPS, temperature);
  return box;
}

/* initial_conditions - sets partilce data for @particles */
ParticleData* initial_conditions(ParticleData *particles, int *number,
      box box) {
  double x_pos, y_pos, z_pos, time_start;
  ParticleType pi("pi", 0.13957);
  ParticleType pi0("pi0", 0.134977);

  printd("Pi^± mass: %g [GeV]\n", pi.mass());
  printd("Pi^0 mass: %g [GeV]\n", pi0.mass());

  /* Set random IC:
   * particles at random position in the box with random momentum
   */
  *number = 50;
  particles = new ParticleData[*number];
  for (int i = 0; i < *number; i++) {
    particles[i].set_id(i);
    particles[i].set_momentum(pi.mass(), randGauss(1.0), randGauss(1.0),
      randGauss(1.0));
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
