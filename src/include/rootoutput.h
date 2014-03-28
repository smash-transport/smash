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

#include "outputinterface.h"
#include <boost/filesystem.hpp>
#include "TFile.h"
#include "TTree.h"

namespace Smash {
class Particles;

class RootOutput : public OutputInterface {
 public:
  RootOutput(boost::filesystem::path path);
  ~RootOutput();

  void at_runstart() override;
  void at_eventstart(const Particles &particles, const int evt_num) override;
  void at_eventend(const Particles &particles, const int evt_num) override;
  //void at_collision(const Collisions &collisions) override;
  void at_outtime(const Particles &particles, const int evt_num, const int timestep) override;
  void at_runend() override;
  void at_crash() override;
  void particles_to_tree(const char* treename, const char* treedescr, const Particles &particles, const int evt_num);

 private:
   const boost::filesystem::path base_path_;
   std::unique_ptr<TFile> root_out_file;
   std::vector<std::unique_ptr<TTree>> tree_list_;
};
}  // namespace Smash

#endif  // SRC_INCLUDE_ROOTOUTPUT_H_
