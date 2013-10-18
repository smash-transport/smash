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
class Sphere;

/* initialisation functions */
void initial_conditions(Particles *particles, Box *box);
void initial_conditions(Particles *particles, Sphere *ball);

#endif  // SRC_INCLUDE_INITIAL_CONDITIONS_H_
