/*
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_OUTPUTROUTINES_H_
#define SRC_INCLUDE_OUTPUTROUTINES_H_

#include "forwarddeclarations.h"
#include "logging.h"

namespace Smash {

/* Compile time debug info */
#ifndef NDEBUG
template <typename... Ts>
void printd(const char *format, Ts &&... args) {
  char tmp[512];
  snprintf(tmp, 512, format, std::forward<Ts>(args)...);
  logger<LogArea::Legacy>().debug(tmp);
}
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
