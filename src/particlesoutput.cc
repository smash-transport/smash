/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/particlesoutput.h"
#include "include/particles.h"
#include "include/filedeleter.h"
#include <memory>

namespace Smash {

ParticlesOutput::ParticlesOutput(boost::filesystem::path path)
    : base_path_(std::move(path)) {
}

ParticlesOutput::~ParticlesOutput() {
}

void ParticlesOutput::write_state(const Particles &particles) {
  char filename[64];

  snprintf(filename, sizeof(filename), "momenta_%.5f.dat",
           particles.time());
  std::unique_ptr<FILE> momenta_file{
      fopen((base_path_ / filename).native().c_str(), "w")};
  for (auto i = particles.cbegin(); i != particles.cend(); ++i) {
    fprintf(momenta_file.get(), "%g %g %g %g %i %i\n", i->second.momentum().x0(),
            i->second.momentum().x1(), i->second.momentum().x2(),
            i->second.momentum().x3(), i->second.id(), i->second.pdgcode());
  }
  snprintf(filename, sizeof(filename), "position_%.5f.dat",
           particles.time());
  std::unique_ptr<FILE> position_file{
      fopen((base_path_ / filename).native().c_str(), "w")};
  for (auto i = particles.cbegin(); i != particles.cend(); ++i) {
    fprintf(position_file.get(), "%g %g %g %g %i %i\n", i->second.position().x0(),
            i->second.position().x1(), i->second.position().x2(),
            i->second.position().x3(), i->second.id(), i->second.pdgcode());
  }
}

}  // namespace Smash
