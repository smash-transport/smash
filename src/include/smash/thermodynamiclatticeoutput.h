/*
 *
 *    Copyright (c) 2014-2021
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_THERMODYNAMICLATTICEOUTPUT_H_
#define SRC_INCLUDE_SMASH_THERMODYNAMICLATTICEOUTPUT_H_

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>

#include "density.h"
#include "experimentparameters.h"
#include "file.h"
#include "forwarddeclarations.h"
#include "outputinterface.h"
#include "outputparameters.h"
#include "threevector.h"

namespace smash {

/**
 * \ingroup output
 *
 * \brief Writes the thermodynamic quantities at lattice points versus time
 *
 * This class is a temporary solution to write thermodynamic
 * quantities out. Calculations are called directly inside the
 * output functions. In future it should be substituted by some
 * more general output.
 *
 **/
class ThermodynamicLatticeOutput : public OutputInterface {
 public:
  /// Version of the thermodynamic lattice output
  static const char* version;

  /// aliases for direction indexes
  static const int xdir, ydir, zdir;

  /**
   * Construct Output
   * param[in] path Path to output
   * param[in] name Filename
   * param[in] out_par Parameters of output
   */
  ThermodynamicLatticeOutput(const bf::path &path, const std::string &name,
                      const OutputParameters &out_par);
  /// Default destructor
  ~ThermodynamicLatticeOutput();
   /**
   * Output launched at event start after initialization, when particles are
   * generated but not yet propagated.
   * \param[in] event_number Number of the current event.
   * \param[in] tq Thermodynamic quantity to deal with
   * \param[in] dens_type Density type for the reference frame.
   * \param[in] lattice Specialized Lattice for DensityOnLattice
   */
  void at_eventstart(const int event_number, const ThermodynamicQuantity tq,
      const DensityType dens_type, const RectangularLattice<DensityOnLattice>
      lattice) override;
/**
   * Output launched at event start after initialization, when particles are
   * generated but not yet propagated.
   * \param[in] event_number Number of the current event.
   * \param[in] tq Thermodynamic quantity to deal with.
   * \param[in] dens_type Density type for the reference frame.
   * \param[in] lattice Specialized Lattice for EnergyMomentumTensor
   */
  void at_eventstart(const int event_number, const ThermodynamicQuantity tq,
      const DensityType dens_type,
      const RectangularLattice<EnergyMomentumTensor>
      lattice) override;

  /**
   *  Closes the output files
   *  \param[in] tq The quantity that has been written in the output file,
   *           see ThermodynamicQuantity.
   */
  void at_eventend(const ThermodynamicQuantity tq) override;

  /**
   * Prints the density lattice on a grid.
   *
   * \param[in] lattice DensityOnLattice lattice to use.
   * \param[in] current_time The output time in the computational frame
   * 
   */
  void thermodynamics_lattice_output(
      RectangularLattice<DensityOnLattice> &lattice,
      const double current_time);

  /**
   * Prints the energy-momentum-tensor lattice on a grid.
   *
   * \param[in] tq The quantity whose energy-momentum tensor should be written,
   *           see ThermodynamicQuantity.
   * \param[in] lattice EnergyMomentumTensor lattice to use.
   * \param[in] current_time The output time in the computational frame
   * 
   */
  void thermodynamics_lattice_output(
      const ThermodynamicQuantity tq,
      RectangularLattice<EnergyMomentumTensor> &lattice,
      const double current_time);


 private:
  /// Structure that holds all the information about what to printout
  const OutputParameters out_par_;

  /**
   * Make a file name given a description and a counter.
   *
   * \param[in] description The description.
   * \param[in] event_number The event number.
   */
  std::string make_filename(const std::string &description,
                            const int event_number);

  /**
   * Make a variable name given quantity and density type.
   *
   * \param tq The quantity.
   * \param dens_type The density type.
   */
  std::string make_varname(const ThermodynamicQuantity tq,
                           const DensityType dens_type);

  /**
   * Write the header.
   *
   * \param file Output file.
   * \param[in] tq The quantity to be written,
   *           see ThermodynamicQuantity.
   */
  void write_therm_lattice_header(std::shared_ptr<std::ofstream> file,
      const ThermodynamicQuantity &tq);

  /**
   * Write a scalar.
   *
   * \param file Output file.
   * \param lat Lattice corresponding to output.
   * \param varname Name of the output variable.
   * \param function Function that gets the scalar given a lattice node.
   */
  template <typename T, typename F>
  void write_therm_lattice_scalar(std::ofstream &file,
      RectangularLattice<T> &lat, const std::string &varname, F &&function);

 /**
   * Write a vector.
   *
   * \param file Output file.
   * \param lat Lattice corresponding to output.
   * \param varname Name of the output variable.
   * \param function Function that gets the vector given a lattice node.
   */
  template <typename T, typename F>
  void write_therm_lattice_vector(std::ofstream &file,
      RectangularLattice<T> &lat, const std::string &varname, F &&function);

  /// filesystem path for output
  const bf::path base_path_;

  /// map of output file handlers
  std::map <ThermodynamicQuantity, std::shared_ptr<std::ofstream>>
       output_files_;

  /// number of nodes in the lattice along the three axes
  int nx_, ny_, nz_;

  /// lattice resolution along the three axes
  double dx_, dy_, dz_;

  /// lattice origin
  /// orientation: if 0,0,0 is the origin of a cube with face widths 10,
  /// the center is at 5,5,5
  double x0_, y0_, z0_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_THERMODYNAMICLATTICEOUTPUT_H_
