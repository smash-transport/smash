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

/* initial_conditions - sets partilce data for @particles */
ParticleData* initial_conditions(ParticleData *particles) {
  int num = 5;
  double x_pos, y_pos, z_pos, time_start;
  FourVector position;
  ParticleType pi("pi", 0.13957);
  ParticleType pi0("pi0", 0.134977);

  printd("Pi^± mass: %g [GeV]\n", pi.mass());
  printd("Pi^0 mass: %g [GeV]\n", pi0.mass());

  /* Set random IC:
   * particles at random position in the box with random momentum
   */
  particles = new ParticleData[num];
  for (int i = 0; i < num; i++) {
    particles[i].set_id(i);
    particles[i].set_momentum(pi.mass(), randGauss(1.0), randGauss(1.0),
      randGauss(1.0));
    time_start = 1.0;
    x_pos = 1.0 * rand_r(&seedp) / RAND_MAX * A;
    y_pos = 1.0 * rand_r(&seedp) / RAND_MAX * A;
    z_pos = 1.0 * rand_r(&seedp) / RAND_MAX * A;
    position.set_FourVector(time_start, z_pos, x_pos, y_pos);
    particles[i].set_position(position);
    printd("Particle %d momentum: %g %g %g %g [GeV]\n", particles[i].id(),
      particles[i].momentum().x0(), particles[i].momentum().x1(),
      particles[i].momentum().x2(), particles[i].momentum().x3());
    printd("Particle %d position: %g %g %g %g\n", particles[i].id(),
      particles[i].x().x0(), particles[i].x().x1(),
      particles[i].x().x2(), particles[i].x().x3());
  }

  return particles;
}
