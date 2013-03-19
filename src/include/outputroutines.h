/*
 *    Copyright (c) 2012 maximilian attems <attems@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_OUTPUTROUTINES_H_
#define SRC_INCLUDE_OUTPUTROUTINES_H_

/* forward declarations */
class box;
class ParticleData;

/* console output */
void print_startup(box box);

/* data directory */
void mkdir_data(void);

/* Compile time debug info */
#ifdef DEBUG
# define printd printf
#else
# define printd(...) ((void)0)
#endif

/* console debug output */
void printd_position(ParticleData particle);
void printd_momenta(ParticleData particle);

/* output data files */
void write_particles(ParticleData *particles, const int number);

#endif  // SRC_INCLUDE_OUTPUTROUTINES_H_
