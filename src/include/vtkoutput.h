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


#include <boost/filesystem.hpp>

#include "configuration.h"
#include "density.h"
#include "forwarddeclarations.h"
#include "outputinterface.h"

namespace Smash {

/**
 * \ingroup output
 * SMASH output in a paraview format,
 * intended for simple visualization.
 */
class VtkOutput : public OutputInterface {
 public:
  VtkOutput(bf::path path, Configuration&& conf);
  ~VtkOutput();

  void at_eventstart(const Particles &particles,
                     const int event_number) override;
  void at_eventend(const Particles &particles, const int event_number) override;
  void at_intermediate_time(const Particles &particles, const int event_number,
                          const Clock &clock) override;

  /// Prints 3D Lattice in vtk format on a grid
  void thermodynamics_output(const std::string varname,
                             RectangularLattice<DensityOnLattice> &lattice,
                             const int event_number) override;

 private:
  void write(const Particles &particles, const int event_number);

  /// filesystem path for output
  const bf::path base_path_;

  /// Number of vtk output in current event
  int vtk_output_counter_;

  /// Number of thermodynamical vtk output in current event
  int vtk_thermodynamics_output_counter_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_VTKOUTPUT_H_
