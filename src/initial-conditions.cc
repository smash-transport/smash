/*
 *
 *    Copyright (c) 2012
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

/* initial_conditions - sets partilce data for @particles */
ParticleData* initial_conditions(ParticleData *particles) {
  int num = 5;
  float x_pos, y_pos, z_pos;
  ParticleType pi("pi", 0.13957);
  ParticleType pi0("pi0", 0.134977);

  printd("Pi^± mass: %g [GeV]\n", pi.mass());
  printd("Pi^0 mass: %g [GeV]\n", pi0.mass());

  /* Set random IC:
   * particles at random position in the box with random momenta
   */
  particles = new ParticleData[num];
  for (int i = 0; i < num; i++) {
    particles[i].set_id(i);
    particles[i].set_momenta(randGauss(1.0), randGauss(1.0));
    x_pos = 1.0 * rand_r(&seedp) / RAND_MAX * A;
    y_pos = 1.0 * rand_r(&seedp) / RAND_MAX * A;
    z_pos = 1.0 * rand_r(&seedp) / RAND_MAX * A;
    particles[i].set_position(x_pos, y_pos, z_pos);
    printd("Particle %d momenta: %g %g [GeV]\n", particles[i].id(),
      particles[i].momenta_l(), particles[i].momenta_t());
    printd("Particle %d position: %g %g %g\n", particles[i].id(),
      particles[i].x(), particles[i].y(), particles[i].z());
  }

  return particles;
}
