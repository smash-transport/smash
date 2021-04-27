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
  static const double_t version;

  /**
   * Construct Output
   * \param[in] path Path to output
   * \param[in] name Filename
   * \param[in] out_par Parameters of output
   * \param[in] enable_ascii Bool (True or False) to enable ASCII format
   * \param[in] enable_binary Bool (True or False) to enable binary format
   */
  ThermodynamicLatticeOutput(const bf::path &path, const std::string &name,
                             const OutputParameters &out_par,
                             const bool enable_ascii, const bool enable_binary);
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
  void at_eventstart(
      const int event_number, const ThermodynamicQuantity tq,
      const DensityType dens_type,
      const RectangularLattice<DensityOnLattice> lattice) override;
  /**
   * Output launched at event start after initialization, when particles are
   * generated but not yet propagated.
   * \param[in] event_number Number of the current event.
   * \param[in] tq Thermodynamic quantity to deal with.
   * \param[in] dens_type Density type for the reference frame.
   * \param[in] lattice Specialized Lattice for EnergyMomentumTensor
   */
  void at_eventstart(
      const int event_number, const ThermodynamicQuantity tq,
      const DensityType dens_type,
      const RectangularLattice<EnergyMomentumTensor> lattice) override;

  /**
   *  Final actions at the end of each event (it closes the output files).
   *  \param[in] tq The quantity that has been written in the output file,
   *           see ThermodynamicQuantity.
   */
  void at_eventend(const ThermodynamicQuantity tq) override;

  /**
   * Prints the density lattice on a grid.   *
   * \param[in] lattice DensityOnLattice lattice to use.
   * \param[in] current_time The output time in the computational frame
   */
  void thermodynamics_lattice_output(
      RectangularLattice<DensityOnLattice> &lattice, double current_time);

  /**
   * Prints the density lattice on a grid.
   *
   * \param[in] lattice DensityOnLattice lattice to use.
   * \param[in] current_time The output time in the computational frame
   * \param[in] ensembles Particles, from which the 4-currents j_{Q,B,S} are
   *            computed
   * * \param[in] dens_param set of parameters, defining smearing.
   *            For more info about
   *            smearing see \ref thermodyn_output_user_guide_.
   */
  void thermodynamics_lattice_output(
      RectangularLattice<DensityOnLattice> &lattice, const double current_time,
      const std::vector<Particles> &ensembles,
      const DensityParameters &dens_param);

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
      RectangularLattice<EnergyMomentumTensor> &lattice, double current_time);

 private:
  /// Structure that holds all the information about what to printout
  const OutputParameters out_par_;

  /**
   * Makes a file name given a description and a counter.
   *
   * \param[in] description The description.
   * \param[in] event_number The event number.
   * \param[in] type Flag for the file type: 'a' for ASCII, 'b' for Binary
   */
  std::string make_filename(const std::string &description,
                            const int event_number, const char type);

  /**
   * Makes a variable name given quantity and density type.
   *
   * \param[in] tq The quantity.
   * \param[in] dens_type The density type.
   */
  std::string make_varname(const ThermodynamicQuantity tq,
                           const DensityType dens_type);

  /**
   * Writes the header for the ASCII output files
   *
   * \param file Output file.
   * \param[in] tq The quantity to be written,
   *           see ThermodynamicQuantity.
   */
  void write_therm_lattice_ascii_header(std::shared_ptr<std::ofstream> file,
                                        const ThermodynamicQuantity &tq);

  /**
   * Writes the header for the binary output files
   *
   * \param[in] file Output file.
   * \param[in] tq The quantity to be written,
   *           see ThermodynamicQuantity.
   */
  void write_therm_lattice_binary_header(std::shared_ptr<std::ofstream> file,
                                         const ThermodynamicQuantity &tq);

  /**
   * Convert a ThermodynamicQuantity into an int
   * \param[in] tq The quantity to be converted, see ThermodynamicQuantity.
   * \return An int corresponding to ThermodynamicQuantity.
   */
  int to_int(const ThermodynamicQuantity &tq);

  /// filesystem path for output
  const bf::path base_path_;

  /// map of output file handlers for ASCII format
  std::map<ThermodynamicQuantity, std::shared_ptr<std::ofstream>>
      output_ascii_files_;

  /// map of output file handlers for binary format
  std::map<ThermodynamicQuantity, std::shared_ptr<std::ofstream>>
      output_binary_files_;

  /// number of nodes in the lattice along the three axes
  int nodes_[3];

  /// lattice resolution along the three axes
  double sizes_[3];

  /// lattice origin
  /// orientation: if 0,0,0 is the origin of a cube with face widths 10,
  /// the center is at 5,5,5
  double origin_[3];

  /// enable output type ASCII
  bool enable_ascii_;

  /// enable output type Binary
  bool enable_binary_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_THERMODYNAMICLATTICEOUTPUT_H_
