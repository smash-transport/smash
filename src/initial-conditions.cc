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

void initial_conditions() {
  ParticleType pi("pi", 0.13957);
  ParticleType pi0("pi0", 0.134977);

  printd("Pi^± mass: %f [GeV]\n", pi.mass());
  printd("Pi^0 mass: %f [GeV]\n", pi0.mass());
}
