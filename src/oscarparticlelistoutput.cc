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
#include "include/configuration.h"

#include <boost/filesystem.hpp>

namespace Smash {

OscarParticleListOutput::OscarParticleListOutput(bf::path path,
                                                 Configuration&& conf)
  : file_{std::fopen((path / "final_id_p_x.oscar").native().c_str(), "w")},
  only_final_(static_cast<bool>(conf.take({"only_final"}))) {
  fprintf(file_.get(), "# OSC1999A\n");
  fprintf(file_.get(), "# final_id_p_x\n");
  fprintf(file_.get(), "# smash\n");
  fprintf(file_.get(), "# Block format:\n");
  fprintf(file_.get(), "# nin nout event_number\n");
  fprintf(file_.get(), "# id pdg 0 px py pz p0 mass x y z t\n");
  fprintf(file_.get(), "# End of event: 0 0 event_number\n");
  fprintf(file_.get(), "#\n");
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

void OscarParticleListOutput::write(const Particles &particles) {
  for (const ParticleData &data : particles.data()) {
    fprintf(file_.get(), "%i %s %i %g %g %g %g %g %g %g %g %g\n", data.id(),
            data.pdgcode().string().c_str(), 0, data.momentum().x1(),
            data.momentum().x2(), data.momentum().x3(), data.momentum().x0(),
            sqrt(data.momentum().Dot(data.momentum())), data.position().x1(),
            data.position().x2(), data.position().x3(),
            data.position().x0());
  }
}

}  // namespace Smash
