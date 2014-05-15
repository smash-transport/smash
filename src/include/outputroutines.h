/*
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_OUTPUTROUTINES_H_
#define SRC_INCLUDE_OUTPUTROUTINES_H_

#include <list>

#include "include/chrono.h"
#include "include/particles.h"
#include "forwarddeclarations.h"

namespace Smash {

/* console output */
void print_startup(const ModusDefault &parameters);
void print_header(void);
void print_measurements(const Particles &particles,
                        const size_t &scatterings_total,
                        const size_t &scatterings_this_interval,
                        float energy_ini,
                 SystemTimePoint time_start);
void print_tail(const
                SystemTimePoint time_start,
                const double &scattering_rate);

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
void printd_list(const std::list<int> &collision_list);

}  // namespace Smash

#endif  // SRC_INCLUDE_OUTPUTROUTINES_H_
