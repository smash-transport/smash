/*
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_OUTPUTROUTINES_H_
#define SRC_INCLUDE_OUTPUTROUTINES_H_

#include <cstdlib>

#include "../include/Particles.h"

/* forward declarations */
class Box;
class Parameters;
class ParticleData;
class ParticleType;

/* console output */
void print_startup(const Parameters &parameters);
void print_startup(const Box &box);
void print_header(void);
void print_measurements(const Particles &particles,
                        const size_t &scatterings_total,
                        const size_t &scatterings_this_interval,
                        const Box &box);
void print_tail(const Box &box, const double &scattering_rate);

/* data directory */
void mkdir_data(void);

/* Compile time debug info */
#ifdef DEBUG
# define printd printf
#else
# define printd(...) ((void)0)
#endif

/* console debug output */
void printd_position(const ParticleData &particle);
void printd_position(const char *message, const ParticleData &particle);
void printd_momenta(const ParticleData &particle);
void printd_momenta(const char *message, const ParticleData &particle);

/* output data files */
void write_particles(const Particles &particles);
void write_measurements_header(const Particles &particles);
void write_measurements(const Particles &particles,
  int interactions_total, int interactions_this_interval, int decays,
  int resonances, const size_t &rejection_conflict);
void write_oscar_header(void);
void write_oscar(const ParticleData &particle_data,
                 const ParticleType &particle_type, int initial, int final);
void write_oscar(const ParticleData &particle_data,
                 const ParticleType &particle_type);
void write_vtk(const Particles &particles);

/* timing measure */
double measure_timediff(const Box &box);

#endif  // SRC_INCLUDE_OUTPUTROUTINES_H_
