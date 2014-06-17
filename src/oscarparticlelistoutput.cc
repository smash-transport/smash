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
                                                 std::string option)
  : OscarFullHistoryOutput(path / "final_id_p_x.oscar", "# final_id_p_x\n",
                           option) {}

OscarParticleListOutput::~OscarParticleListOutput() {}

void OscarParticleListOutput::at_eventstart(const Particles &particles,
                                            const int event_number) {
  /* OSCAR line prefix : initial particles; final particles; event id
   * First block of an event: initial = 0, final = number of particles
   */
  if (config_option_ != "Final") {
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
  if (config_option_ == "Final") {
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
  if (config_option_ != "Final") {
    const size_t zero = 0;
    fprintf(file_.get(), "%zu %zu %i\n", particles.size(), zero,
            event_number + 1);
    write(particles);
  }
}

}  // namespace Smash
