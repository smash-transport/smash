/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_PARTICLESOUTPUT_H_
#define SRC_INCLUDE_PARTICLESOUTPUT_H_

#include "outputinterface.h"
#include "forwarddeclarations.h"
#include <boost/filesystem.hpp>

namespace Smash {
class Particles;

class ParticlesOutput : public OutputInterface {
 public:
  ParticlesOutput(bf::path path);
  ~ParticlesOutput();

  void at_eventstart(const Particles &particles, const int) override {
    write_state(particles);
  }

  void at_eventend(const Particles &particles, const int) override {
    write_state(particles);
  }

  void after_Nth_timestep(const Particles &particles, const int,
                          const int) override {
    write_state(particles);
  }

 private:
  void write_state(const Particles &particles);

  const bf::path base_path_;
};
}  // namespace Smash

#endif  // SRC_INCLUDE_PARTICLESOUTPUT_H_
