/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/oscarparticlelistoutput.h"

#include <boost/filesystem.hpp>

#include "include/clock.h"
#include "include/particles.h"
#include "include/outputroutines.h"
#include "include/configuration.h"


namespace Smash {

OscarParticleListOutput::OscarParticleListOutput(bf::path path,
                                                 Configuration&& conf)
  : OscarOutput(path, "particle_lists.oscar", std::move(conf)),
    only_final_(conf.has_value({"Only_final"})
              ? conf.take({"Only_final"}) : true)  {
  if (modern_format_) {
    fprintf(file_.get(), "#!OSCAR2013 particle_lists ");
    fprintf(file_.get(), "t x y z mass p0 px py pz pdg ID\n");
    fprintf(file_.get(),
            "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none\n");
  } else {
    fprintf(file_.get(), "# OSC1999A\n");
    fprintf(file_.get(), "# final_id_p_x\n");
    fprintf(file_.get(), "# smash\n");
    fprintf(file_.get(), "# Block format:\n");
    fprintf(file_.get(), "# nin nout event_number\n");
    fprintf(file_.get(), "# id pdg 0 px py pz p0 mass x y z t\n");
    fprintf(file_.get(), "# End of event: 0 0 event_number\n");
    fprintf(file_.get(), "#\n");
  }
}

OscarParticleListOutput::~OscarParticleListOutput() {}

void OscarParticleListOutput::at_eventstart(const Particles &particles,
                                            const int event_number) {
  if (!only_final_) {
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

void OscarParticleListOutput::at_eventend(const Particles &particles,
                              const int event_number) {
  if (modern_format_) {
    if (only_final_) {
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
    if (only_final_) {
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

void OscarParticleListOutput::after_Nth_timestep(const Particles &particles,
                                                 const int event_number,
                                                 const Clock&/*clock*/) {
  if (!only_final_) {
    if (modern_format_) {
      fprintf(file_.get(), "# event %i out %zu\n", event_number + 1,
              particles.size());
    } else {
      const size_t zero = 0;
      fprintf(file_.get(), "%zu %zu %i\n", particles.size(), zero,
              event_number + 1);
    }
    write(particles);
  }
}

}  // namespace Smash
