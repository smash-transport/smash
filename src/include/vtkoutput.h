/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_VTKOUTPUT_H_
#define SRC_INCLUDE_VTKOUTPUT_H_

#include <string>

#include <boost/filesystem.hpp>

#include "density.h"
#include "forwarddeclarations.h"
#include "outputinterface.h"

namespace smash {

/**
 * \ingroup output
 * SMASH output in a paraview format, intended for simple visualization.
 */
class VtkOutput : public OutputInterface {
 public:
  /**
   * Create a new VTK output.
   *
   * \param path Path to the output file.
   * \param name Name of the output.
   */
  VtkOutput(const bf::path &path, const std::string &name);
  ~VtkOutput();

  /**
   * Writes the initial particle information list of an event to the VTK
   * output.
   *
   * \param particles Current list of all particles.
   * \param event_number Number of the current event.
   */
  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  /**
   * Writes the final particle information list of an event to the VTK
   * output. This currently does not do anything, because it is not
   * required for the VTK output.
   *
   * \param particles Unused. Current list of particles.
   * \param event_number Unused. Number of event.
   * \param impact_parameter Unused. Impact parameter of this event.
   */
  void at_eventend(const Particles &particles, const int event_number,
                   double impact_parameter) override;

  /**
   * Writes out all current particles.
   *
   * \param particles Current list of particles.
   * \param clock Unused, needed since inherited.
   * \param dens_param Unused, needed since inherited.
   */
  void at_intermediate_time(const Particles &particles, const Clock &clock,
                            const DensityParameters &dens_param) override;

  /**
   * Prints the density lattice in VTK format on a grid.
   *
   * \param tq The quantity whose density should be written.
   * \param dt The type of the density.
   * \param lattice The lattice from which the quantity is taken.
   */
  void thermodynamics_output(
      const ThermodynamicQuantity tq, const DensityType dt,
      RectangularLattice<DensityOnLattice> &lattice) override;

  /**
   * Prints the energy-momentum-tensor lattice in VTK format on a grid.
   *
   * \param tq The quantity whose density should be written.
   * \param dt The type of the density.
   * \param lattice The lattice from which the quantity is taken.
   */
  void thermodynamics_output(
      const ThermodynamicQuantity tq, const DensityType dt,
      RectangularLattice<EnergyMomentumTensor> &lattice) override;

  /**
   * Printout of all thermodynamic quantities from the thermalizer class.
   *
   * \param gct Grand-canonical thermalizer from which the quantities are
   *            taken.
   */
  void thermodynamics_output(const GrandCanThermalizer &gct) override;

 private:
  /**
   * Write the given particles to the output.
   *
   * \param particles The particles.
   */
  void write(const Particles &particles);

  /**
   * Make a file name given a description and a counter.
   *
   * \param description The description.
   * \param counter The counter enumerating the outputs.
   */
  std::string make_filename(const std::string &description, int counter);

  /**
   * Make a variable name given quantity and density type.
   *
   * \param tq The quantity.
   * \param dens_type The density type.
   */
  std::string make_varname(const ThermodynamicQuantity tq,
                           const DensityType dens_type);

  /**
   * Write the VTK header.
   *
   * \param file Output file.
   * \param lat Lattice corresponding to output.
   * \param description Description of the output.
   */
  template <typename T>
  void write_vtk_header(std::ofstream &file, RectangularLattice<T> &lat,
                        const std::string &description);

  /**
   * Write a VTK scalar.
   *
   * \param file Output file.
   * \param lat Lattice corresponding to output.
   * \param varname Name of the output variable.
   * \param function Function that gets the scalar given a lattice node.
   */
  template <typename T, typename F>
  void write_vtk_scalar(std::ofstream &file, RectangularLattice<T> &lat,
                        const std::string &varname, F &&function);

  /**
   * Write a VTK vector.
   *
   * \param file Output file.
   * \param lat Lattice corresponding to output.
   * \param varname Name of the output variable.
   * \param function Function that gets the vector given a lattice node.
   */
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
  /// Is the VTK output a thermodynamics output
  bool is_thermodynamics_output_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_VTKOUTPUT_H_
