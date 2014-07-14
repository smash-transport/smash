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
#include "forwarddeclarations.h"

#include <boost/filesystem.hpp>

namespace Smash {

class VtkOutput : public OutputInterface {
 public:
  VtkOutput(bf::path path);
  ~VtkOutput();

  void at_eventstart(const Particles &particles,
                     const int event_number) override;
  void at_eventend(const Particles &particles, const int event_number) override;
  void after_Nth_timestep(const Particles &particles, const int event_number,
                          const Clock &clock) override;

 private:
  const bf::path base_path_;
  /// Number of vtk output in current event
  int vtk_output_counter_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_VTKOUTPUT_H_
