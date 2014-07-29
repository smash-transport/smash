/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "include/oscaroutput.h"

#include <boost/filesystem.hpp>
#include <string>

#include "include/clock.h"
#include "include/forwarddeclarations.h"
#include "include/particles.h"
#include "include/outputroutines.h"
#include "include/configuration.h"

namespace Smash {

OscarOutput::OscarOutput(bf::path path, std::string filename,
                           Configuration &&conf)
  : file_{std::fopen((path / filename).native().c_str(), "w")},
    modern_format_(conf.has_value({"2013_format"})
                   ? conf.take({"2013_format"}) : false) {}

OscarOutput::~OscarOutput() {}

void OscarOutput::write_format_description() {
  if (modern_format_) {
    fprintf(file_.get(), "t x y z mass p0 px py pz pdg ID\n");
    fprintf(file_.get(),
            "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none\n");
  } else {
    fprintf(file_.get(), "# Block format:\n");
    fprintf(file_.get(), "# nin nout event_number\n");
    fprintf(file_.get(), "# id pdg 0 px py pz p0 mass x y z t\n");
    fprintf(file_.get(), "# End of event: 0 0 event_number\n");
    fprintf(file_.get(), "#\n");
  }
}

void OscarOutput::write_particledata(const ParticleData &data) {
  if (modern_format_) {
    fprintf(file_.get(), "%g %g %g %g %g %g %g %g %g %s %i\n",
            data.position().x0(), data.position().x1(),
            data.position().x2(), data.position().x3(),
            std::sqrt(data.momentum().Dot(data.momentum())),
            data.momentum().x0(), data.momentum().x1(),
            data.momentum().x2(), data.momentum().x3(),
            data.pdgcode().string().c_str(), data.id());
  } else {
    fprintf(file_.get(), "%i %s %i %g %g %g %g %g %g %g %g %g\n", data.id(),
            data.pdgcode().string().c_str(), 0, data.momentum().x1(),
            data.momentum().x2(), data.momentum().x3(), data.momentum().x0(),
            std::sqrt(data.momentum().Dot(data.momentum())),
            data.position().x1(), data.position().x2(), data.position().x3(),
            data.position().x0());
  }
}

void OscarOutput::write(const Particles &particles) {
  for (const ParticleData &data : particles.data()) {
    write_particledata(data);
  }
}

void OscarOutput::at_eventstart(const Particles& /*particles*/,
                                const int /*event_number*/) {
  /* Behavior of this function depends on subclass */
}

void OscarOutput::at_eventend(const Particles& /*particles*/,
                              const int /*event_number*/) {
  /* Behavior of this function depends on subclass */
}

void OscarOutput::at_interaction(
  const ParticleList& /*incoming_particles*/,
  const ParticleList& /*outgoing_particles*/) {
  /* Behavior of this function depends on subclass */
}

void OscarOutput::at_intermediate_time(const Particles & /*particles*/,
                                                const int /*event_number*/,
                                                const Clock& /*clock*/) {
  /* Behavior of this function depends on subclass */
}

}  // namespace Smash
