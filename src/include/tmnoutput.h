/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_TMNOUTPUT_H_
#define SRC_INCLUDE_TMNOUTPUT_H_

#include "outputinterface.h"
#include <boost/filesystem.hpp>
#include "TFile.h"
#include "TTree.h"

namespace Smash {
class Particles;

class TmnOutput : public OutputInterface {
 public:
  TmnOutput(boost::filesystem::path path);
  ~TmnOutput();

  void at_runstart() override;
  void at_eventstart(const Particles &particles, const int evt_num) override;
  void at_eventend(const Particles &particles, const int evt_num) override;
  //void at_collision(const Collisions &collisions) override;
  void at_outtime(const Particles &particles, const int evt_num, const int timestep) override;
  void at_runend() override;
  void at_crash() override;
  void Tmn_to_tree(const char* treename, const char* treedescr, const Particles &particles, const int evt_num);

 private:
   const boost::filesystem::path base_path_;
   TFile* root_out_file;
   TTree* curr_tree;
   double p0, px,py,pz;
   double t00,t01,t02,t03,t11,t12,t13,t22,t23,t33;
   int    ev;
};
}  // namespace Smash

#endif  // SRC_INCLUDE_ROOTOUTPUT_H_
