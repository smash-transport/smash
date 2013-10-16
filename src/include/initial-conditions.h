/*
 *    Copyright (c) 2012 maximilian attems <attems@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_INITIAL_CONDITIONS_H_
#define SRC_INCLUDE_INITIAL_CONDITIONS_H_

/* forward declarations */
class Box;
class Laboratory;
class Particles;

/* initialisation functions */
void initial_particles(Particles *particles);
void initial_conditions(Particles *particles,
  const Laboratory &parameters, Box *box);

#endif  // SRC_INCLUDE_INITIAL_CONDITIONS_H_
