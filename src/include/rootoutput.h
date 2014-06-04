/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_ROOTOUTPUT_H_
#define SRC_INCLUDE_ROOTOUTPUT_H_

#include <boost/filesystem.hpp>
#include <vector>
#include "outputinterface.h"
#include "TFile.h"
#include "TTree.h"

namespace Smash {
class Particles;

class RootOutput : public OutputInterface {
 public:
  explicit RootOutput(boost::filesystem::path path);
  ~RootOutput();

  void at_eventstart(const Particles &particles,
                     const int event_number) override;
  void at_eventend(const Particles &particles,
                   const int event_number) override;
  void after_Nth_timestep(const Particles &particles,
                          const int event_number,
                          const Clock &) override;
  void write_interaction(const ParticleList &incoming_particles,
                         const ParticleList &outgoing_particles) override;

 private:
  const boost::filesystem::path base_path_;
  std::unique_ptr<TFile> root_out_file_;
  // TFile takes ownership of all TTrees.
  // That's why TTrees are not unique pointers.
  std::vector<TTree*> tree_list_;
  void particles_to_tree(const char* treename,
                         const char* treedescr,
                         const Particles &particles,
                         const int event_number);
};
}  // namespace Smash

#endif  // SRC_INCLUDE_ROOTOUTPUT_H_
