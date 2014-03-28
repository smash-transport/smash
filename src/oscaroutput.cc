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
    : base_path_(std::move(path)) {
  write_oscar_header();
}

OscarOutput::~OscarOutput() {}

void OscarOutput::at_eventstart(const Particles &particles,
                                const int event_number) {
  /* Write the initial data block of the event */
  write_oscar_event_block(&particles, 0, particles.size(), event_number + 1);
}

void OscarOutput::at_eventend(const Particles &particles,
                              const int event_number) {
  /* Write the final data block of the event */
  write_oscar_event_block(&particles, particles.size(), 0, event_number + 1);
}

void OscarOutput::before_collision() {}

void OscarOutput::after_collision() {}

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
