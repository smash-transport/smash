/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "include/oscarfullhistoryoutput.h"

#include <boost/filesystem.hpp>
#include <string>

#include <include/config.h>
#include "include/clock.h"
#include "include/forwarddeclarations.h"
#include "include/particles.h"
#include "include/outputroutines.h"
#include "include/configuration.h"

namespace Smash {

OscarFullHistoryOutput::OscarFullHistoryOutput(bf::path path,
                                               Configuration &&conf)
  : OscarOutput(path, "full_event_history.oscar", std::move(conf)),
    print_start_end_(conf.has_value({"Print_start_end"})
                     ? conf.take({"Print_start_end"}) : false) {
  if (modern_format_) {
    fprintf(file_.get(), "#!OSCAR2013 full_event_history %s ", VERSION_MAJOR);
  } else {
    fprintf(file_.get(), "# OSC1999A\n");
    fprintf(file_.get(), "# full_event_history\n");
    fprintf(file_.get(), "# %s\n", VERSION_MAJOR);
  }
  write_format_description();
}


OscarFullHistoryOutput::~OscarFullHistoryOutput() {}

void OscarFullHistoryOutput::at_eventstart(const Particles &particles,
                                const int event_number) {
  if (print_start_end_) {
    if (modern_format_) {
        fprintf(file_.get(), "# event %i in %zu\n", event_number + 1,
                particles.size());
    } else {
      /* OSCAR line prefix : initial particles; final particles; event id
       * First block of an event: initial = 0, final = number of particles
       */
      const size_t zero = 0;
      fprintf(file_.get(), "%zu %zu %i\n", zero, particles.size(),
              event_number + 1);
    }
    write(particles);
  }
}

void OscarFullHistoryOutput::at_eventend(const Particles &particles,
                              const int event_number) {
  if (modern_format_) {
    if (print_start_end_) {
      fprintf(file_.get(), "# event %i out %zu\n", event_number + 1,
              particles.size());
      write(particles);
    }
    /* Comment end of an event */
    fprintf(file_.get(), "# event %i end\n", event_number + 1);
  } else {
    /* OSCAR line prefix : initial particles; final particles; event id
     * Last block of an event: initial = number of particles, final = 0
     * Block ends with null interaction
     */
    const size_t zero = 0;
    if (print_start_end_) {
      fprintf(file_.get(), "%zu %zu %i\n", particles.size(), zero,
              event_number + 1);
      write(particles);
    }
    /* Null interaction marks the end of an event */
    fprintf(file_.get(), "%zu %zu %i\n", zero, zero, event_number + 1);
  }
  /* Flush to disk */
  std::fflush(file_.get());
}

void OscarFullHistoryOutput::at_interaction(
  const ParticleList &incoming_particles,
  const ParticleList &outgoing_particles) {
  if (modern_format_) {
    fprintf(file_.get(), "# interaction in %zu out %zu\n",
            incoming_particles.size(), outgoing_particles.size());
  } else {
    /* OSCAR line prefix : initial final
     * particle creation: 0 1
     * particle 2<->2 collision: 2 2
     * resonance formation: 2 1
     * resonance decay: 1 2
     * etc.
     */
    fprintf(file_.get(), "%zu %zu\n", incoming_particles.size(),
            outgoing_particles.size());
  }
  for (const auto &p : incoming_particles) {
    write_particledata(p);
  }
  for (const auto &p : outgoing_particles) {
    write_particledata(p);
  }
}

}  // namespace Smash
