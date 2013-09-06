/*
 *    Copyright (c) 2012 maximilian attems <attems@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_INPUT_PARTICLES_H_
#define SRC_INCLUDE_INPUT_PARTICLES_H_

#include <vector>

/* forward declarations */
class Particles;

extern const char *sep;

/* read input file particle types */
void input_particles(Particles *particles, char *path);

#endif  // SRC_INCLUDE_INPUT_PARTICLES_H_
