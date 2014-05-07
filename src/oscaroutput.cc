/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/oscaroutput.h"
#include "include/particles.h"
#include "include/outputroutines.h"

namespace Smash {

OscarOutput::OscarOutput(boost::filesystem::path path)
    : file_{std::fopen((path / "collision.dat").native().c_str(), "w")} {
  fprintf(file_.get(), "# OSC1999A\n");
  fprintf(file_.get(), "# Interaction history\n");
  fprintf(file_.get(), "# smash\n");
  fprintf(file_.get(), "#\n");
}

OscarOutput::~OscarOutput() {}

void OscarOutput::at_eventstart(const Particles &particles,
                                const int event_number) {
  /* OSCAR line prefix : initial particles; final particles; event id
   * First block of an event: initial = 0, final = number of particles
   * Vice versa for the last block
   */
  const size_t zero = 0;
  fprintf(file_.get(), "%zu %zu %i\n", zero, particles.size(), event_number + 1);
  write(particles);
}

void OscarOutput::at_eventend(const Particles &particles,
                              const int event_number) {
  /* OSCAR line prefix : initial particles; final particles; event id
   * First block of an event: initial = 0, final = number of particles
   * Vice versa for the last block
   */
  const size_t zero = 0;
  fprintf(file_.get(), "%zu %zu %i\n", particles.size(), zero, event_number + 1);
  write(particles);
}

void OscarOutput::write(const Particles &particles) {
  for (const ParticleData &data : particles.data()) {
    fprintf(file_.get(), "%i %s %i %g %g %g %g %g %g %g %g %g\n", data.id(),
            data.pdgcode().string().c_str(), 0, data.momentum().x1(),
            data.momentum().x2(), data.momentum().x3(), data.momentum().x0(),
            sqrt(data.momentum().Dot(data.momentum())), data.position().x1(),
            data.position().x2(), data.position().x3(),
            data.position().x0() - 1.0);
  }
}

void OscarOutput::before_collision() {}

void OscarOutput::after_collision() {}

void OscarOutput::write_interaction(const ParticleList &incoming_particles,
                                    const ParticleList &outgoing_particles) {
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
            p.position().x0() - 1.0);
  };
  for (const auto &p : incoming_particles) {
    print(p);
  }
  for (const auto &p : outgoing_particles) {
    print(p);
  }
}

void OscarOutput::after_Nth_timestep(const Particles & /*particles*/,
                                     const int /*event_number*/,
                                     const int /*timestep*/) {
  /*
  char filename[64];

  snprintf(filename, sizeof(filename), "momenta_%.5f.dat",
           particles.time());
  std::unique_ptr<FILE> momenta_file{
      fopen((base_path_ / filename).native().c_str(), "w")};
  for (const ParticleData &data : particles.data()) {
    fprintf(momenta_file.get(), "%g %g %g %g %i %s\n",
            data.momentum().x0(),
            data.momentum().x1(), data.momentum().x2(),
            data.momentum().x3(), data.id(), data.pdgcode().string().c_str());
  }
  snprintf(filename, sizeof(filename), "position_%.5f.dat",
           particles.time());
  std::unique_ptr<FILE> position_file{
      fopen((base_path_ / filename).native().c_str(), "w")};
  for (const ParticleData &data : particles.data()) {
    fprintf(position_file.get(), "%g %g %g %g %i %s\n",
            data.position().x0(),
            data.position().x1(), data.position().x2(),
            data.position().x3(), data.id(), data.pdgcode().string().c_str());
  }
  */
}

}  // namespace Smash
