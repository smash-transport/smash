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

OscarParticleListOutput::OscarParticleListOutput(bf::path path)
  : OscarFullHistoryOutput(path / "final_id_p_x.oscar", "# final_id_p_x\n") {}

OscarParticleListOutput::~OscarParticleListOutput() {}

void OscarParticleListOutput::at_eventstart(const Particles &/*particles*/,
                                            const int /*event_number*/) {
  /* No initial state output */
}

void OscarParticleListOutput::write_interaction(
  const ParticleList &/*incoming_particles*/,
  const ParticleList &/*outgoing_particles*/) {

  /* No interaction output */
}

}  // namespace Smash
