/*
 *    Copyright (c) 2012 
 *      SMASH Team
 * 
 *    GNU General Public License (GPLv3 or later)
 * 
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
