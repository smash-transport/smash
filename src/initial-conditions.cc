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

void initial_conditions(ParticleData *particles) {
  int num = 5;
  ParticleType pi("pi", 0.13957);
  ParticleType pi0("pi0", 0.134977);

  printd("Pi^± mass: %g [GeV]\n", pi.mass());
  printd("Pi^0 mass: %g [GeV]\n", pi0.mass());

  particles = new ParticleData[num];
  for (int i = 0; i < num; i++) {
    particles[i].set_id(i);
    particles[i].set_momenta(randGauss(1.0), randGauss(1.0));
    printd("Particle %d momenta: %g %g [GeV]\n", particles[i].id(),
      particles[i].momenta_l(), particles[i].momenta_t());
  }
}
