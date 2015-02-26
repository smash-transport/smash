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

  /** Prints 3D density map in vtk format on a grid [-nx;nx]x[-ny;ny]x[-nz;nz]
   *  with steps dx, dy, dz. This allows to look at density profiles and
   *  make easy plots.
   */
  void vtk_density_map(const char * file_name, const ParticleList &plist,
                      double gs_sigma, Density_type dens_type, int ntest,
                      int nx, int ny, int nz, double dx, double dy, double dz);


 private:
  void write(const Particles &particles, const int event_number);

  /// filesystem path for output
  const bf::path base_path_;

  /// Number of vtk output in current event
  int vtk_output_counter_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_VTKOUTPUT_H_
