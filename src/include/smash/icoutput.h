/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_ICOUTPUT_H_
#define SRC_INCLUDE_ICOUTPUT_H_

#include <string>

#include <boost/filesystem.hpp>

#include "file.h"
#include "outputinterface.h"
#include "outputparameters.h"
#include "smash/config.h"

namespace smash {

/**
 * \ingroup output
 * SMASH output in ASCII format containing initial conditions for hydrodynamic
 * codes.
 */
class ICOutput : public OutputInterface {
 public:
  /**
   * Create a new IC output.
   *
   * \param path Path to the output file.
   * \param name Name of the output.
   * \param out_par Additional information on the configured output.
   */
  ICOutput(const bf::path &path, const std::string &name,
           const OutputParameters &out_par);
  ~ICOutput();

  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  void at_eventend(const Particles &particles, const int event_number,
                   double impact_parameter, bool empty_event) override;
  /**
   * Writes out all current particles at correct proper time and removes
   * them from the total list of particles.
   *
   * \param particles Current list of particles.
   * \param clock Current time in
   * \param dens_param Unused, needed since inherited.
   */
  void at_intermediate_time(const Particles &particles, const Clock &clock,
                            const DensityParameters &dens_param) override;

  void at_interaction(const Action &action, const double density) override;
  // void at_hypersurface_crossing(const Particles &particles);

  /**
   * Prints the density lattice in VTK format on a grid.
   *
   * \param tq The quantity whose density should be written,
   *           see ThermodynamicQuantity.
   * \param dt The type of the density, see DensityType.
   * \param lattice The lattice from which the quantity is taken.
   */
  // void thermodynamics_output(
  //     const ThermodynamicQuantity tq, const DensityType dt,
  //     RectangularLattice<DensityOnLattice> &lattice) override;

  /**
   * Prints the energy-momentum-tensor lattice in VTK format on a grid.
   *
   * \param tq The quantity whose density should be written,
               see ThermodynamicQuantity.
   * \param dt The type of the density, see DensityType
   * \param lattice The lattice from which the quantity is taken.
   */
  // void thermodynamics_output(
  //     const ThermodynamicQuantity tq, const DensityType dt,
  //     RectangularLattice<EnergyMomentumTensor> &lattice) override;
  //
  /**
   * Printout of all thermodynamic quantities from the thermalizer class.
   *
   *  \param gct Grand-canonical thermalizer from which the quantities are
   *            taken.
   */
  // void thermodynamics_output(const GrandCanThermalizer &gct) override;

 private:
  /**
   * Write the given particles to the output.
   *
   * \param particles The particles.
   */
  // void write(const Particles &particles);

  /**
   * Make a file name given the type (particles/Tmunu)
   *
   * \param type The type.
   */
  // std::string make_filename(const std::string &type);

  /**
   * Make a variable name given quantity and density type.
   *
   * \param tq The quantity.
   * \param dens_type The density type.
   */
  // std::string make_varname(const ThermodynamicQuantity tq,
  //                          const DensityType dens_type);

  /**
   * Write the VTK header.
   *
   * \param file Output file.
   * \param lat Lattice corresponding to output.
   * \param description Description of the output.
   */
  // template <typename T>
  // void write_vtk_header(std::ofstream &file, RectangularLattice<T> &lat,
  //                       const std::string &description);

  // /// filesystem path for output
  // const bf::path base_path_;
  //
  // /// Event number
  // int current_event_ = 0;
  // /// Number of IC output in current event
  // //int IC_output_counter_ = 0;
  //
  // /// Is Ic output?
  // bool is_IC_output_;

  /// Pointer to output file
  RenamingFilePtr file_;
  /// Structure that holds all the information about what to printout
  const OutputParameters out_par_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_ICOUTPUT_H_
