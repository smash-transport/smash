/*
 *    Copyright (c) 2012 maximilian attems <attems@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_OUTPUTROUTINES_H_
#define SRC_INCLUDE_OUTPUTROUTINES_H_

/* forward declaration */
class ParticleData;

/* console output */
void print_startup(void);

/* data files */
void print_particles(ParticleData *particles, const int number);

#endif  // SRC_INCLUDE_OUTPUTROUTINES_H_
