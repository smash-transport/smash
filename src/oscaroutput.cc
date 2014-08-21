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

#include <include/config.h>
#include "include/clock.h"
#include "include/forwarddeclarations.h"
#include "include/particles.h"
#include "include/outputroutines.h"
#include "include/configuration.h"

namespace Smash {

template <OscarOutputFormat Format, int Contents>
OscarOutput<Format, Contents>::OscarOutput(bf::path path, std::string name)
    : file_{std::fopen((path / (name + ".oscar")).native().c_str(), "w")} {
  if (Format == OscarFormat2013) {
    fprintf(file_.get(), "#!OSCAR2013 %s %s ", name.c_str(), VERSION_MAJOR);
    fprintf(file_.get(), "t x y z mass p0 px py pz pdg ID\n");
    fprintf(file_.get(),
            "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none\n");
  } else {
    if (name == "particle_lists") {
      name = "final_id_p_x";  // FIXME: why is this necessary? I.e. what does
                              // the string on the second line tell, and why
                              // does it have to be this specific string?
    }
    fprintf(file_.get(), "# OSC1999A\n# %s\n# %s\n", name.c_str(), VERSION_MAJOR);
    fprintf(file_.get(), "# Block format:\n");
    fprintf(file_.get(), "# nin nout event_number\n");
    fprintf(file_.get(), "# id pdg 0 px py pz p0 mass x y z t\n");
    fprintf(file_.get(), "# End of event: 0 0 event_number\n");
    fprintf(file_.get(), "#\n");
  }
}

template <OscarOutputFormat Format, int Contents>
inline void OscarOutput<Format, Contents>::write(const Particles &particles) {
  for (const ParticleData &data : particles.data()) {
    write_particledata(data);
  }
}

template <OscarOutputFormat Format, int Contents>
void OscarOutput<Format, Contents>::at_eventstart(const Particles &particles,
                                                  const int event_number) {
  if (Contents & OscarAtEventstart) {
    if (Format == OscarFormat2013) {
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

template <OscarOutputFormat Format, int Contents>
void OscarOutput<Format, Contents>::at_eventend(const Particles &particles,
                                                const int event_number) {
  if (Format == OscarFormat2013) {
    if (Contents & OscarParticlesAtEventend) {
      fprintf(file_.get(), "# event %i out %zu\n", event_number + 1,
              particles.size());
      write(particles);
    }
    // Comment end of an event
    fprintf(file_.get(), "# event %i end\n", event_number + 1);
  } else {
    // OSCAR line prefix : initial particles; final particles; event id
    // Last block of an event: initial = number of particles, final = 0
    // Block ends with null interaction
    const size_t zero = 0;
    if (Contents & OscarParticlesAtEventend) {
      fprintf(file_.get(), "%zu %zu %i\n", particles.size(), zero,
              event_number + 1);
      write(particles);
    }
    // Null interaction marks the end of an event
    fprintf(file_.get(), "%zu %zu %i\n", zero, zero, event_number + 1);
  }
  // Flush to disk
  std::fflush(file_.get());
}

template <OscarOutputFormat Format, int Contents>
void OscarOutput<Format, Contents>::at_interaction(
    const ParticleList &incoming_particles,
    const ParticleList &outgoing_particles) {
  if (Contents & OscarInteractions) {
    if (Format == OscarFormat2013) {
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
}

template <OscarOutputFormat Format, int Contents>
void OscarOutput<Format, Contents>::at_intermediate_time(
    const Particles &particles, const int event_number,
    const Clock & /*clock*/) {
  if (Contents & OscarTimesteps) {
    if (Format == OscarFormat2013) {
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

template <OscarOutputFormat Format, int Contents>
void OscarOutput<Format, Contents>::write_particledata(
    const ParticleData &data) {
  if (Format == OscarFormat2013) {
    fprintf(file_.get(), "%g %g %g %g %g %g %g %g %g %s %i\n",
            data.position().x0(), data.position().x1(), data.position().x2(),
            data.position().x3(),
            std::sqrt(data.momentum().Dot(data.momentum())),
            data.momentum().x0(), data.momentum().x1(), data.momentum().x2(),
            data.momentum().x3(), data.pdgcode().string().c_str(), data.id());
  } else {
    fprintf(file_.get(), "%i %s %i %g %g %g %g %g %g %g %g %g\n", data.id(),
            data.pdgcode().string().c_str(), 0, data.momentum().x1(),
            data.momentum().x2(), data.momentum().x3(), data.momentum().x0(),
            std::sqrt(data.momentum().Dot(data.momentum())),
            data.position().x1(), data.position().x2(), data.position().x3(),
            data.position().x0());
  }
}

/**
 * Will be in C++14's standard library
 */
template <typename T, typename... Args>
std::unique_ptr<T> make_unique(Args &&... args) {
  return std::unique_ptr<T>{new T{std::forward<Args>(args)...}};
}

std::unique_ptr<OutputInterface> create_oscar_output(bf::path path,
                                                     Configuration config) {
  if (config.has_value({"OSCAR_PARTICLELIST", "Enable"})) {
    auto subconfig = config["OSCAR_PARTICLELIST"];
    const bool enabled = subconfig.take({"Enable"});
    if (!enabled) {
      config.take({"OSCAR_PARTICLELIST"});
    } else {
      const bool modern_format = subconfig.has_value({"2013_format"})
                                     ? subconfig.take({"2013_format"})
                                     : false;
      const bool only_final = subconfig.has_value({"Only_final"})
                                  ? subconfig.take({"Only_final"})
                                  : true;
      // combine the two bools to get four possible states
      switch (modern_format * 2 + only_final) {
        case 0:  // modern_format == only_final == false
          return make_unique<
              OscarOutput<OscarFormat1999, OscarTimesteps | OscarAtEventstart>>(
              std::move(path), "particle_lists");
        case 1:  // modern_format == false && only_final == true
          return make_unique<
              OscarOutput<OscarFormat1999, OscarParticlesAtEventend>>(
              std::move(path), "particle_lists");
        case 2:  // modern_format == true  && only_final == false
          return make_unique<
              OscarOutput<OscarFormat2013, OscarTimesteps | OscarAtEventstart>>(
              std::move(path), "particle_lists");
        case 3:  // modern_format == only_final == true
          return make_unique<
              OscarOutput<OscarFormat2013, OscarParticlesAtEventend>>(
              std::move(path), "particle_lists");
        default:
          throw std::runtime_error(
              "impossible Oscar output configuration requested");
      }
    }
  }
  if (config.has_value({"OSCAR_COLLISIONS", "Enable"})) {
    auto subconfig = config["OSCAR_COLLISIONS"];
    const bool enabled = subconfig.take({"Enable"});
    if (!enabled) {
      config.take({"OSCAR_COLLISIONS"});
    } else {
      const bool modern_format = subconfig.has_value({"2013_format"})
                                     ? subconfig.take({"2013_format"})
                                     : false;
      const bool print_start_end = subconfig.has_value({"Print_start_end"})
                                  ? subconfig.take({"Print_start_end"})
                                  : false;
      // combine the two bools to get four possible states
      switch (modern_format * 2 + print_start_end) {
        case 0:  // modern_format == print_start_end == false
          return make_unique<OscarOutput<OscarFormat1999, OscarInteractions>>(
              std::move(path), "full_event_history");
        case 1:  // modern_format == false && print_start_end == true
          return make_unique<OscarOutput<OscarFormat1999,
                                         OscarInteractions | OscarAtEventstart |
                                             OscarParticlesAtEventend>>(
              std::move(path), "full_event_history");
        case 2:  // modern_format == true  && print_start_end == false
          return make_unique<OscarOutput<OscarFormat2013, OscarInteractions>>(
              std::move(path), "full_event_history");
        case 3:  // modern_format == print_start_end == true
          return make_unique<OscarOutput<OscarFormat2013,
                                         OscarInteractions | OscarAtEventstart |
                                             OscarParticlesAtEventend>>(
              std::move(path), "full_event_history");
        default:
          throw std::runtime_error(
              "impossible Oscar output configuration requested");
      }
    }
  }
  return {};  // return a nullptr to signify the end of OSCAR outputs in the
              // config file
}

}  // namespace Smash
