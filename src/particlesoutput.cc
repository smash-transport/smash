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
}

}  // namespace Smash
