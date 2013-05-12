/*
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_OUTPUTROUTINES_H_
#define SRC_INCLUDE_OUTPUTROUTINES_H_

#include <cstdlib>
#include <map>

/* forward declarations */
class box;
class ParticleData;
class ParticleType;

/* console output */
void print_startup(const box &box);
void print_header(void);
void print_measurements(const ParticleData *particle,
  const int &num, const size_t &scatterings_total, const box &box);
void print_tail(const box &box, const double &scattering_rate);

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
void printd_momenta(const ParticleData &particle);

/* output data files */
void write_particles(ParticleData *particles, const int number);
void write_oscar_header(void);
void write_oscar(const ParticleData &particle1, const ParticleData &particle2,
  const ParticleType &type1, const ParticleType &type2, int flag);
void write_vtk(const ParticleData *particles, int number);

/* timing measure */
double measure_timediff(const box &box);

#endif  // SRC_INCLUDE_OUTPUTROUTINES_H_
