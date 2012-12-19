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
void initial_conditions(ParticleData *particles) {
  int num = 5;
  ParticleType pi("pi", 0.13957);
  ParticleType pi0("pi0", 0.134977);

  printd("Pi^± mass: %g [GeV]\n", pi.mass());
  printd("Pi^0 mass: %g [GeV]\n", pi0.mass());

  /* Set random IC */
  particles = new ParticleData[num];
  for (int i = 0; i < num; i++) {
    particles[i].set_id(i);
    particles[i].set_momenta(randGauss(1.0), randGauss(1.0));
    particles[i].set_position(randGauss(0.2), randGauss(0.2), randGauss(0.2));
    printd("Particle %d momenta: %g %g [GeV]\n", particles[i].id(),
      particles[i].momenta_l(), particles[i].momenta_t());
    printd("Particle %d position: %g %g %g\n", particles[i].id(),
      particles[i].x(), particles[i].y(), particles[i].z());
  }
}
