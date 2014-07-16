/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <string>

#include "include/clock.h"
#include "include/forwarddeclarations.h"
#include "include/oscarfullhistoryoutput.h"
#include "include/particles.h"
#include "include/outputroutines.h"

#include <boost/filesystem.hpp>

namespace Smash {

OscarFullHistoryOutput::OscarFullHistoryOutput(bf::path path, Options op)
  : OscarFullHistoryOutput(path / "full_event_history.oscar",
                           "full_event_history", op){}

OscarFullHistoryOutput::OscarFullHistoryOutput(bf::path path,
                                               const char* format_specifier,
                                               Options op)
  : file_{std::fopen(path.native().c_str(), "w")},
    modern_format_(false), write_initial_final_lists_(false) {

  std::string opt_str;
  if (op.count("Print_start_end") != 0) {
    opt_str = op.at("Print_start_end");
    for (auto &c : opt_str) {
      c = tolower(c);
    }
    if (opt_str == "true") {
      write_initial_final_lists_ = true;
    }
  }
  if (op.count("Format") != 0) {
    if (op.at("Format") == "2013") {
      modern_format_ = true;
    }
  }
  if (modern_format_) {
    fprintf(file_.get(), "#!OSCAR2013 %s ", format_specifier);
    fprintf(file_.get(), "t x y z mass p0 px py pz pdg ID\n");
    fprintf(file_.get(),
            "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none\n");
  } else {
    fprintf(file_.get(), "# OSC1999A\n");
    fprintf(file_.get(), "# %s\n", format_specifier);
    fprintf(file_.get(), "# smash\n");
    fprintf(file_.get(), "# Block format:\n");
    fprintf(file_.get(), "# nin nout event_number\n");
    fprintf(file_.get(), "# id pdg 0 px py pz p0 mass x y z t\n");
    fprintf(file_.get(), "# End of event: 0 0 event_number\n");
    fprintf(file_.get(), "#\n");
  }
}


OscarFullHistoryOutput::~OscarFullHistoryOutput() {}

void OscarFullHistoryOutput::at_eventstart(const Particles &particles,
                                const int event_number) {
  if (write_initial_final_lists_) {
    if (modern_format_) {
        fprintf(file_.get(), "# event %i in %zu\n", event_number + 1,
                particles.size());
        write_2013(particles);
    } else {
      /* OSCAR line prefix : initial particles; final particles; event id
       * First block of an event: initial = 0, final = number of particles
       */
      const size_t zero = 0;
      fprintf(file_.get(), "%zu %zu %i\n", zero, particles.size(),
              event_number + 1);
      write(particles);
    }
  }
}

void OscarFullHistoryOutput::at_eventend(const Particles &particles,
                              const int event_number) {
  if (modern_format_) {
    if (write_initial_final_lists_) {
      fprintf(file_.get(), "# event %i out %zu\n", event_number + 1,
              particles.size());
      write_2013(particles);
    }
    /* Comment end of an event */
    fprintf(file_.get(), "# event %i end\n", event_number + 1);
  } else {
    /* OSCAR line prefix : initial particles; final particles; event id
     * Last block of an event: initial = number of particles, final = 0
     * Block ends with null interaction
     */
    const size_t zero = 0;
    if (write_initial_final_lists_) {
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

void OscarFullHistoryOutput::write(const Particles &particles) {
  for (const ParticleData &data : particles.data()) {
    fprintf(file_.get(), "%i %s %i %g %g %g %g %g %g %g %g %g\n", data.id(),
            data.pdgcode().string().c_str(), 0, data.momentum().x1(),
            data.momentum().x2(), data.momentum().x3(), data.momentum().x0(),
            sqrt(data.momentum().Dot(data.momentum())), data.position().x1(),
            data.position().x2(), data.position().x3(),
            data.position().x0());
  }
}

void OscarFullHistoryOutput::write_2013(const Particles &particles) {
  for (const ParticleData &data : particles.data()) {
    fprintf(file_.get(), "%g %g %g %g %g %g %g %g %g %s %i\n",
            data.position().x0(), data.position().x1(),
            data.position().x2(), data.position().x3(),
            std::sqrt(data.momentum().Dot(data.momentum())),
            data.momentum().x0(), data.momentum().x1(),
            data.momentum().x2(), data.momentum().x3(),
            data.pdgcode().string().c_str(), data.id());
  }
}

void OscarFullHistoryOutput::write_interaction(
  const ParticleList &incoming_particles,
  const ParticleList &outgoing_particles) {

  if (modern_format_) {
    fprintf(file_.get(), "# interaction in %zu out %zu\n",
            incoming_particles.size(), outgoing_particles.size());
    const auto print = [&](const ParticleData &p) {
      const float mass = std::sqrt(p.momentum().Dot(p.momentum()));
      fprintf(file_.get(), "%g %g %g %g %g %g %g %g %g %s %i\n",
              p.position().x0(), p.position().x1(),
              p.position().x2(), p.position().x3(),
              mass, p.momentum().x0(), p.momentum().x1(),
              p.momentum().x2(), p.momentum().x3(),
              p.pdgcode().string().c_str(), p.id());
    };
    for (const auto &p : incoming_particles) {
      print(p);
    }
    for (const auto &p : outgoing_particles) {
      print(p);
    }
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
    const auto print = [&](const ParticleData &p) {
      const float mass = std::sqrt(p.momentum().Dot(p.momentum()));
      fprintf(file_.get(), "%i %s %i %g %g %g %g %g %g %g %g %g\n", p.id(),
              p.pdgcode().string().c_str(), 0, p.momentum().x1(),
              p.momentum().x2(), p.momentum().x3(), p.momentum().x0(), mass,
              p.position().x1(), p.position().x2(), p.position().x3(),
              p.position().x0());
    };
    for (const auto &p : incoming_particles) {
      print(p);
    }
    for (const auto &p : outgoing_particles) {
      print(p);
    }
  }
}

void OscarFullHistoryOutput::after_Nth_timestep(const Particles & /*particles*/,
                                                const int /*event_number*/,
                                     const Clock& /*clock*/) {
  /* No time interval output for interaction history */
}

}  // namespace Smash
