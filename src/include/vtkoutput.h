/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_VTKOUTPUT_H_
#define SRC_INCLUDE_VTKOUTPUT_H_

#include "outputinterface.h"
#include <boost/filesystem.hpp>
class Particles;

namespace Smash {

class VtkOutput : public OutputInterface {
 public:
   VtkOutput(boost::filesystem::path path);
   ~VtkOutput();

   void at_eventstart(const Particles &particles, const int event_number) override;
   void at_eventend(const Particles &particles, const int event_number) override;
   void after_collision() override;
   void before_collision() override;
   void after_Nth_timestep(const Particles &particles, const int event_number, const int timestep) override;


 private:
   const boost::filesystem::path base_path_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_VTKOUTPUT_H_
