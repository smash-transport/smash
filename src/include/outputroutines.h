/*
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_OUTPUTROUTINES_H_
#define SRC_INCLUDE_OUTPUTROUTINES_H_

#include <list>

#include "chrono.h"
#include "forwarddeclarations.h"
#include "quantumnumbers.h"

namespace Smash {

/* console output */
void print_startup(const ModusDefault &parameters);
void print_header(void);
void print_measurements(const Particles &particles,
                        const size_t scatterings_total,
                        const size_t scatterings_this_interval,
                        const QuantumNumbers& conserved_initial,
                        const SystemTimePoint time_start,
                        const float time);
void print_tail(const SystemTimePoint time_start,
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

}  // namespace Smash

#endif  // SRC_INCLUDE_OUTPUTROUTINES_H_
