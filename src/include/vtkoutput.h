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

#include <string>

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
  VtkOutput(const bf::path &path, Configuration &&conf);
  ~VtkOutput();

  void at_eventstart(const Particles &particles,
                     const int event_number) override;
  void at_eventend(const Particles &particles, const int event_number) override;
  void at_intermediate_time(const Particles &particles, const Clock &clock,
                            const DensityParameters &dens_param) override;

  /// Prints 3D Lattice in vtk format on a grid
  void thermodynamics_output(
      const ThermodynamicQuantity tq, const DensityType dt,
      RectangularLattice<DensityOnLattice> &lattice) override;

  /// Prints 3D Lattice in vtk format on a grid
  void thermodynamics_output(
      const ThermodynamicQuantity tq, const DensityType dt,
      RectangularLattice<EnergyMomentumTensor> &lattice) override;

  /// Printout of the thermodynamic quantities from thethermalizer class
  void thermodynamics_output(const GrandCanThermalizer &gct) override;

 private:
  void write(const Particles &particles);

  std::string make_filename(const std::string &description, int counter);

  std::string make_varname(const ThermodynamicQuantity tq,
                           const DensityType dens_type);

  template <typename T>
  void write_vtk_header(std::ofstream &file, RectangularLattice<T> &lat,
                        const std::string &description);

  template <typename T, typename F>
  void write_vtk_scalar(std::ofstream &file, RectangularLattice<T> &lat,
                        const std::string &varname, F &&function);

  template <typename T, typename F>
  void write_vtk_vector(std::ofstream &file, RectangularLattice<T> &lat,
                        const std::string &varname, F &&function);

  /// filesystem path for output
  const bf::path base_path_;

  /// Event number
  int current_event_ = 0;
  /// Number of vtk output in current event
  int vtk_output_counter_ = 0;

  /// Number of density lattice vtk output in current event
  int vtk_density_output_counter_ = 0;
  /// Number of energy-momentum tensor lattice vtk output in current event
  int vtk_tmn_output_counter_ = 0;
  /// Number of Landau frame energy-momentum tensor vtk output in current event
  int vtk_tmn_landau_output_counter_ = 0;
  /// Number of Landau rest frame velocity vtk output in current event
  int vtk_v_landau_output_counter_ = 0;
  /// Number of fluidization output
  int vtk_fluidization_counter_ = 0;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_VTKOUTPUT_H_
