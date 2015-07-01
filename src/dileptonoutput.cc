/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "include/dileptonoutput.h"

#include <boost/filesystem.hpp>
#include <string>

#include "include/config.h"
#include "include/configuration.h"
#include "include/cxx14compat.h"
#include "include/forwarddeclarations.h"
#include "include/particles.h"
#include "include/pdgcode.h"

namespace Smash {

DileptonOutput::DileptonOutput(bf::path path)
: file_{std::fopen((path / "dilepton_out.txt").native().c_str(), "w")} {
  std::fprintf(file_.get(), "# DILEPTONS pdg p0 p1 p2 p3 \n" );
  std::fprintf(file_.get(), "# Units: none GeV GeV GeV GeV\n" );
  std::fprintf(file_.get(), "# parent pdg mass weight\n" ); // TODO the more
}

void DileptonOutput::at_eventstart(const Particles &, const int event_number) {
  std::fprintf(file_.get(), "# event %i in\n", event_number + 1);
  // write(particles); ???
}

void DileptonOutput::at_eventend(const Particles &, const int event_number) {
  std::fprintf(file_.get(), "# event %i out\n", event_number + 1);

  // write(particles); ???
  // Flush to disk
  std::fflush(file_.get());  // why?
}

void DileptonOutput::dileptons(const ParticleList &incoming_particles,
                       const ParticleList &outgoing_particles,
                       float shining_weight){
    std::fprintf(file_.get(), "parent particle %s %g %g\n",
                incoming_particles[0].pdgcode().string().c_str(),
                std::sqrt(incoming_particles[0].momentum().Dot(incoming_particles[0].momentum())),
                shining_weight);
    std::fprintf(file_.get(), "%s %g %g %g %g\n",
            outgoing_particles[0].pdgcode().string().c_str(),
            outgoing_particles[0].momentum().x0(),
            outgoing_particles[0].momentum().x1(),
            outgoing_particles[0].momentum().x2(),
            outgoing_particles[0].momentum().x3());
    std::fprintf(file_.get(), "%s %g %g %g %g\n",
            outgoing_particles[1].pdgcode().string().c_str(),
            outgoing_particles[1].momentum().x0(),
            outgoing_particles[1].momentum().x1(),
            outgoing_particles[1].momentum().x2(),
            outgoing_particles[1].momentum().x3());
}








}  // namespace Smash
