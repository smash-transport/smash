/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/clock.h"
#include "include/oscarfullhistoryoutput.h"
#include "include/oscarparticlelistoutput.h"
#include "include/particles.h"
#include "include/outputroutines.h"

#include <boost/filesystem.hpp>

namespace Smash {

OscarParticleListOutput::OscarParticleListOutput(bf::path path,
                                                 Options op)
  : OscarFullHistoryOutput(path / "final_id_p_x.oscar", "# final_id_p_x\n",
                                                                    op),
  only_final_(true) {
  std::string opt_str = op["only_final"];
  for (auto &c : opt_str) {
    c = tolower(c);
  }
  if (opt_str == "false") {
    only_final_=false;
  }
}

OscarParticleListOutput::~OscarParticleListOutput() {}

void OscarParticleListOutput::at_eventstart(const Particles &particles,
                                            const int event_number) {
  /* OSCAR line prefix : initial particles; final particles; event id
   * First block of an event: initial = 0, final = number of particles
   */
  if (!only_final_) {
    const size_t zero = 0;
    fprintf(file_.get(), "%zu %zu %i\n", zero, particles.size(),
            event_number + 1);
    write(particles);
  }
}

void OscarParticleListOutput::at_eventend(const Particles &particles,
                              const int event_number) {
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

  /* Flush to disk */
  std::fflush(file_.get());
}

void OscarParticleListOutput::write_interaction(
  const ParticleList &/*incoming_particles*/,
  const ParticleList &/*outgoing_particles*/) {

  /* No interaction output */
}

void OscarParticleListOutput::after_Nth_timestep(const Particles &particles,
                                                 const int event_number,
                                                 const Clock&/*clock*/) {
  if (!only_final_) {
    const size_t zero = 0;
    fprintf(file_.get(), "%zu %zu %i\n", particles.size(), zero,
            event_number + 1);
    write(particles);
  }
}

}  // namespace Smash
