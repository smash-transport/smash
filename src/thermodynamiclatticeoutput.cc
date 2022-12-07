/*
 *
 *    Copyright (c) 2014-2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/thermodynamiclatticeoutput.h"

#include <filesystem>
#include <fstream>
#include <memory>

#include "smash/clock.h"
#include "smash/config.h"
#include "smash/density.h"
#include "smash/energymomentumtensor.h"
#include "smash/experimentparameters.h"
#include "smash/thermodynamicoutput.h"
#include "smash/vtkoutput.h"

namespace smash {

/*!\Userguide
 * \page doxypage_output_thermodyn_lattice
 *
 * The thermodynamics lattice output prints "smeared" thermodynamic quantities
 *  evaluated at the nodes of a Lattice, defined as defined in
 * \ref doxypage_input_conf_lattice.
 *
 * The calculated quantities can include: the Eckart density, the
 * energy-momentum tensor in the lab and/or Landau frame, the Landau
 * velocity and the electric/baryonic/strange currents.
 * Which of these quantities are outputted needs to be specified in the config
 * file. They will be calculated using **one** of the following density types:
 * hadron, baryon, baryonic isospin, pion, or none.
 *
 * A separate output file is printed for each quantity and each event.
 *
 * It is possible to print the output in:
 * - ASCII format (option "Lattice_ASCII")
 * - Binary format (option "Lattice_Binary")
 *
 * For example:
 *\verbatim
 Output:
   Output_Interval: 4.0
   Thermodynamics:
     Format: ["Lattice_ASCII","Lattice_Binary"]
     Type: "baryon"
     Quantities: ["rho_eckart", "tmn", "tmn_landau","landau_velocity","j_QBS"]
 \endverbatim
 *
 * Please, note that the Thermodynamic output is not printed at the end of an
   event, even if the end time of the simulation is a multiple of
   Output_Interval (if this final output is important, just increase a little
   the parameter End_Time).
 *
 * **Output files**
 * Each file has a header and a payload. The content is the same both in
 * Lattice_ASCII and Lattice_Binary formats, but, of course, the details are a
 * bit different.
 *
 *
 * **Header of the output files**
 *
 * - version of the output (ASCII: fixed 2 digits float, Binary: double)
 * - thermodynamic quantity (ASCII: string, Binary: int)
 * - nx, ny, nz:  cells of the lattice along x, y, z, respectively
 *   (set in: Lattice->Cell_Number, see \ref doxypage_input_conf_lattice) (3 ints)
 * - x0, y0, z0: coordinates of the origin of the lattice
 *   (set in Lattice->Origin)
 *   (ASCII: 3 fixed 6 digits precision floats, Binary: 3 doubles)
 * - dx, dy, dz: size of the lattice (set in Lattice->Sizes)
 *   (ASCII: 3 fixed 6 digits precision floats, Binary: 3 doubles)
 *
 * In "Lattice_ASCII" the data entries are separated by 1 space and each line of
 * the header starts with "#", followed by the description of the entry.
 *
 * The numbers corresponding to the various thermodynamic quantities are:
 * -0 EckartDensity
 * -1 Tmn (energy momentum tensor in the lab frame)
 * -2 TmnLandau (energy momentum tensor in Landau's frame)
 * -3 LandauVelocity
 * -4 j_QBS (electric charge, baryon and strangeness four currents)
 * -99 Unknown quantity
 *
 * **Payload of the output files**
 *
 * All the data in the payload are represented as:
 * - ASCII files: scientific format, 14 precision digits
 * - Binary files: double type
 *
 * In the case of the energy-momentum tensor, the data payload consists
 * in the values of the quantity in the following order:
 *
 \code
  for (h=0;h<{number of timesteps};h++) {
    output time
    for (l=0;l<{number of components};l++) {
     for (k=0;k<nz;k++) {
      for (j=0;j<ny;j++) {
       for (i=0;i<nx;i++) {
            energy momentum tensor component T[l]
       }
       \newline (only in ASCII files)
      }
     }
    }
  }
 \endcode
 * (see SMASH documentation of the function tmn_index for more details about the
 * components)
 *
 *
 * In the case of densities, the data payload consists
 * in the values of the quantity in the following order:
 *
 \code
  for (h=0;h<{number of timesteps};h++) {
    output time
    for (k=0;k<nz;k++) {
     for (j=0;j<ny;j++) {
      for (i=0;i<nx;i++) {
        density value
      }
      \newline (only in ASCII files)
     }
    }
  }
 \endcode
 *
 * In the case of Landau velocity the data payload
 * consists in the values of the components
 * Vx, Vy, Vz in the following order:

\code
 for (h=0;h<{number of timesteps};h++) {
    output time
    for (k=0;k<nz;k++) {
     for (j=0;j<ny;j++) {
      for (i=0;i<nx;i++) {
        Vx  Vy  Vz (+ \newline in ASCII files)
      }
     }
    }
  }
 \endcode
 *
 * In the case of the electric charge, baryon and strangeness currents
 * jQ, jB and jS the data payload consists in the values of their components
 * in following order:

\code
  for (h=0;h<{number of timesteps};h++) {
    output time
    for (k=0;k<nz;k++) {
     for (j=0;j<ny;j++) {
      for (i=0;i<nx;i++) {
        jQ[0] jQ[1] jQ[2] jQ[3] jB[0] jB[1] jB[2] jB[3] jS[0] jS[1] jS[2] jS[3]
        (+ \newline in ASCII files)
      }
     }
    }
  }
 \endcode
 * (see SMASH documentation of the function current_eckart for more details
 * about
 * the definition of the currents)
 * Please, note that note that all currents are given in units of
   "number of charges"; multiply the electric current by the
   elementary charge \f$\sqrt{4 \pi \alpha_{EM}} \f$ for charge units.
 *
 * Please, have a look also at \ref input_output_thermodynamics_ for additional
 * information about the computation of the various Thermodynamics quantities.
 */

/* initialization of the static member version */
const double_t ThermodynamicLatticeOutput::version = 1.0;

ThermodynamicLatticeOutput::ThermodynamicLatticeOutput(
    const std::filesystem::path &path, const std::string &name,
    const OutputParameters &out_par, const bool enable_ascii,
    const bool enable_binary)
    : OutputInterface(name),
      out_par_(out_par),
      base_path_(std::move(path)),
      enable_ascii_(enable_ascii),
      enable_binary_(enable_binary) {
  if (enable_ascii_ || enable_binary_) {
    enable_output_ = true;
  } else {
    enable_output_ = false;
  }
  if (enable_ascii_) {
    if (out_par_.td_rho_eckart) {
      output_ascii_files_[ThermodynamicQuantity::EckartDensity] =
          std::make_shared<std::ofstream>(nullptr);
    }
    if (out_par_.td_tmn_landau) {
      output_ascii_files_[ThermodynamicQuantity::TmnLandau] =
          std::make_shared<std::ofstream>(nullptr);
    }
    if (out_par_.td_tmn) {
      output_ascii_files_[ThermodynamicQuantity::Tmn] =
          std::make_shared<std::ofstream>(nullptr);
    }
    if (out_par_.td_v_landau) {
      output_ascii_files_[ThermodynamicQuantity::LandauVelocity] =
          std::make_shared<std::ofstream>(nullptr);
    }
    if (out_par_.td_jQBS) {
      output_ascii_files_[ThermodynamicQuantity::j_QBS] =
          std::make_shared<std::ofstream>(nullptr);
    }
  }
  if (enable_binary_) {
    if (out_par_.td_rho_eckart) {
      output_binary_files_[ThermodynamicQuantity::EckartDensity] =
          std::make_shared<std::ofstream>(nullptr);
    }
    if (out_par_.td_tmn_landau) {
      output_binary_files_[ThermodynamicQuantity::TmnLandau] =
          std::make_shared<std::ofstream>(nullptr);
    }
    if (out_par_.td_tmn) {
      output_binary_files_[ThermodynamicQuantity::Tmn] =
          std::make_shared<std::ofstream>(nullptr);
    }
    if (out_par_.td_v_landau) {
      output_binary_files_[ThermodynamicQuantity::LandauVelocity] =
          std::make_shared<std::ofstream>(nullptr);
    }
    if (out_par_.td_jQBS) {
      output_binary_files_[ThermodynamicQuantity::j_QBS] =
          std::make_shared<std::ofstream>(nullptr);
    }
  }
}

ThermodynamicLatticeOutput::~ThermodynamicLatticeOutput() {}

void ThermodynamicLatticeOutput::at_eventstart(
    const int event_number, const ThermodynamicQuantity tq,
    const DensityType dens_type, RectangularLattice<DensityOnLattice> lattice) {
  if (!enable_output_) {
    return;
  }
  assert((tq == ThermodynamicQuantity::EckartDensity) ||
         (tq == ThermodynamicQuantity::j_QBS));
  // at the next refactoring of the code,
  // this piece should go in the constructor
  const auto dim = lattice.n_cells();
  const auto cs = lattice.cell_sizes();
  const auto orig = lattice.origin();
  for (int l = 0; l < 3; l++) {
    nodes_[l] = dim[l];
    sizes_[l] = cs[l];
    origin_[l] = orig[l];
  }
  std::shared_ptr<std::ofstream> fp(nullptr);
  std::string varname;
  std::string filename;
  varname = make_varname(tq, dens_type);
  if (tq == ThermodynamicQuantity::EckartDensity) {
    if (enable_ascii_) {
      filename = make_filename(varname, event_number, 'a');
      try {
        output_ascii_files_[ThermodynamicQuantity::EckartDensity]->open(
            filename, std::ios::out);
      } catch (std::ofstream::failure &e) {
        logg[LogArea::Main::id].fatal()
            << "Error in opening " << filename << std::endl;
        throw std::runtime_error(
            "Not possible to write thermodynamic "
            "lattice output to file.");
      }
      fp = output_ascii_files_[ThermodynamicQuantity::EckartDensity];
      write_therm_lattice_ascii_header(fp, tq);
    }
    if (enable_binary_) {
      filename = make_filename(varname, event_number, 'b');
      try {
        output_binary_files_[ThermodynamicQuantity::EckartDensity]->open(
            filename, std::ios::out | std::ios::binary);
      } catch (std::ofstream::failure &e) {
        logg[LogArea::Main::id].fatal()
            << "Error in opening " << filename << std::endl;
        throw std::runtime_error(
            "Not possible to write thermodynamic "
            "lattice output to file.");
      }
      fp = output_binary_files_[ThermodynamicQuantity::EckartDensity];
      write_therm_lattice_binary_header(fp, tq);
    }
  } else {
    if (enable_ascii_) {
      filename = make_filename(varname, event_number, 'a');
      try {
        output_ascii_files_[ThermodynamicQuantity::j_QBS]->open(filename,
                                                                std::ios::out);
      } catch (std::ofstream::failure &e) {
        logg[LogArea::Main::id].fatal()
            << "Error in opening " << filename << std::endl;
        throw std::runtime_error(
            "Not possible to write thermodynamic "
            "lattice output to file.");
      }
      fp = output_ascii_files_[ThermodynamicQuantity::j_QBS];
      write_therm_lattice_ascii_header(fp, tq);
    }
    if (enable_binary_) {
      filename = make_filename(varname, event_number, 'b');
      try {
        output_binary_files_[ThermodynamicQuantity::j_QBS]->open(
            filename, std::ios::out | std::ios::binary);
      } catch (std::ofstream::failure &e) {
        logg[LogArea::Main::id].fatal()
            << "Error in opening " << filename << std::endl;
        throw std::runtime_error(
            "Not possible to write thermodynamic "
            "lattice output to file.");
      }
      fp = output_binary_files_[ThermodynamicQuantity::j_QBS];
      write_therm_lattice_binary_header(fp, tq);
    }
  }
}

void ThermodynamicLatticeOutput::at_eventstart(
    const int event_number, const ThermodynamicQuantity tq,
    const DensityType dens_type,
    RectangularLattice<EnergyMomentumTensor> lattice) {
  if (!enable_output_) {
    return;
  }
  const auto dim = lattice.n_cells();
  const auto cs = lattice.cell_sizes();
  const auto orig = lattice.origin();
  for (int l = 0; l < 3; l++) {
    nodes_[l] = dim[l];
    sizes_[l] = cs[l];
    origin_[l] = orig[l];
  }
  std::shared_ptr<std::ofstream> fp(nullptr);
  std::string varname;
  std::string filename;
  varname = make_varname(tq, dens_type);
  if (enable_ascii_) {
    filename = make_filename(varname, event_number, 'a');
    if (tq == ThermodynamicQuantity::TmnLandau) {
      try {
        output_ascii_files_[ThermodynamicQuantity::TmnLandau]->open(
            filename, std::ios::out);
      } catch (std::ofstream::failure &e) {
        logg[LogArea::Main::id].fatal()
            << "Error in opening " << filename << std::endl;
        throw std::runtime_error(
            "Not possible to write thermodynamic "
            "lattice output to file.");
      }
      fp = output_ascii_files_[ThermodynamicQuantity::TmnLandau];
    } else if (tq == ThermodynamicQuantity::Tmn) {
      try {
        output_ascii_files_[ThermodynamicQuantity::Tmn]->open(filename,
                                                              std::ios::out);
      } catch (std::ofstream::failure &e) {
        logg[LogArea::Main::id].fatal()
            << "Error in opening " << filename << std::endl;
        throw std::runtime_error(
            "Not possible to write thermodynamic "
            "lattice output to file.");
      }
      fp = output_ascii_files_[ThermodynamicQuantity::Tmn];
    } else if (tq == ThermodynamicQuantity::LandauVelocity) {
      try {
        output_ascii_files_[ThermodynamicQuantity::LandauVelocity]->open(
            filename, std::ios::out);
      } catch (std::ofstream::failure &e) {
        logg[LogArea::Main::id].fatal()
            << "Error in opening " << filename << std::endl;
        throw std::runtime_error(
            "Not possible to write thermodynamic "
            "lattice output to file.");
      }
      fp = output_ascii_files_[ThermodynamicQuantity::LandauVelocity];
    } else {
      try {
        output_ascii_files_[ThermodynamicQuantity::j_QBS]->open(filename,
                                                                std::ios::out);
      } catch (std::ofstream::failure &e) {
        logg[LogArea::Main::id].fatal()
            << "Error in opening " << filename << std::endl;
        throw std::runtime_error(
            "Not possible to write thermodynamic "
            "lattice output to file.");
      }
    }
    write_therm_lattice_ascii_header(fp, tq);
  }
  if (enable_binary_) {
    filename = make_filename(varname, event_number, 'b');
    if (tq == ThermodynamicQuantity::TmnLandau) {
      try {
        output_binary_files_[ThermodynamicQuantity::TmnLandau]->open(
            filename, std::ios::out | std::ios::binary);
      } catch (std::ofstream::failure &e) {
        logg[LogArea::Main::id].fatal()
            << "Error in opening " << filename << std::endl;
        throw std::runtime_error(
            "Not possible to write thermodynamic "
            "lattice output to file.");
      }
      fp = output_binary_files_[ThermodynamicQuantity::TmnLandau];
    } else if (tq == ThermodynamicQuantity::Tmn) {
      try {
        output_binary_files_[ThermodynamicQuantity::Tmn]->open(
            filename, std::ios::out | std::ios::binary);
      } catch (std::ofstream::failure &e) {
        logg[LogArea::Main::id].fatal()
            << "Error in opening " << filename << std::endl;
        throw std::runtime_error(
            "Not possible to write thermodynamic "
            "lattice output to file.");
      }
      fp = output_binary_files_[ThermodynamicQuantity::Tmn];
    } else if (tq == ThermodynamicQuantity::LandauVelocity) {
      try {
        output_binary_files_[ThermodynamicQuantity::LandauVelocity]->open(
            filename, std::ios::out | std::ios::binary);
      } catch (std::ofstream::failure &e) {
        logg[LogArea::Main::id].fatal()
            << "Error in opening " << filename << std::endl;
        throw std::runtime_error(
            "Not possible to write thermodynamic "
            "lattice output to file.");
      }
      fp = output_binary_files_[ThermodynamicQuantity::LandauVelocity];
    } else {
      try {
        output_binary_files_[ThermodynamicQuantity::j_QBS]->open(filename,
                                                                 std::ios::out);
      } catch (std::ofstream::failure &e) {
        logg[LogArea::Main::id].fatal()
            << "Error in opening " << filename << std::endl;
        throw std::runtime_error(
            "Not possible to write thermodynamic "
            "lattice output to file.");
      }
    }
    write_therm_lattice_binary_header(fp, tq);
  }
}

void ThermodynamicLatticeOutput::at_eventend(const ThermodynamicQuantity tq) {
  if (!enable_output_) {
    return;
  }
  if (tq == ThermodynamicQuantity::EckartDensity) {
    if (enable_ascii_) {
      output_ascii_files_[ThermodynamicQuantity::EckartDensity]->close();
    }
    if (enable_binary_) {
      output_binary_files_[ThermodynamicQuantity::EckartDensity]->close();
    }
    return;
  }
  if (tq == ThermodynamicQuantity::Tmn) {
    if (enable_ascii_) {
      output_ascii_files_[ThermodynamicQuantity::Tmn]->close();
    }
    if (enable_binary_) {
      output_binary_files_[ThermodynamicQuantity::Tmn]->close();
    }
    return;
  }
  if (tq == ThermodynamicQuantity::TmnLandau) {
    if (enable_ascii_) {
      output_ascii_files_[ThermodynamicQuantity::TmnLandau]->close();
    }
    if (enable_binary_) {
      output_binary_files_[ThermodynamicQuantity::TmnLandau]->close();
    }
    return;
  }
  if (tq == ThermodynamicQuantity::LandauVelocity) {
    if (enable_ascii_) {
      output_ascii_files_[ThermodynamicQuantity::LandauVelocity]->close();
    }
    if (enable_binary_) {
      output_binary_files_[ThermodynamicQuantity::LandauVelocity]->close();
    }
    return;
  }
  if (tq == ThermodynamicQuantity::j_QBS) {
    if (enable_ascii_) {
      output_ascii_files_[ThermodynamicQuantity::j_QBS]->close();
    }
    if (enable_binary_) {
      output_binary_files_[ThermodynamicQuantity::j_QBS]->close();
    }
    return;
  }
}

void ThermodynamicLatticeOutput::thermodynamics_lattice_output(
    RectangularLattice<DensityOnLattice> &lattice, double ctime) {
  double result;
  const auto dim = lattice.n_cells();
  std::shared_ptr<std::ofstream> fp(nullptr);
  if (enable_ascii_) {
    fp = output_ascii_files_[ThermodynamicQuantity::EckartDensity];
    *fp << std::setprecision(14);
    *fp << std::scientific;
    *fp << ctime << std::endl;
  }
  if (enable_binary_) {
    fp = output_binary_files_[ThermodynamicQuantity::EckartDensity];
    assert(sizeof(ctime) == sizeof(double));
    fp->write(reinterpret_cast<char *>(&ctime), sizeof(ctime));
  }
  lattice.iterate_sublattice(
      {0, 0, 0}, dim, [&](DensityOnLattice &node, int ix, int, int) {
        if (enable_ascii_) {
          *fp << node.rho() << " ";
          if (ix == dim[0] - 1) {
            *fp << "\n";
          }
        }
        if (enable_binary_) {
          result = node.rho();
          fp->write(reinterpret_cast<char *>(&result), sizeof(double));
        }
      });
}

void ThermodynamicLatticeOutput::thermodynamics_lattice_output(
    RectangularLattice<DensityOnLattice> &lattice, double ctime,
    const std::vector<Particles> &ensembles,
    const DensityParameters &dens_param) {
  if (!enable_output_) {
    return;
  }
  double result;
  const auto dim = lattice.n_cells();
  std::shared_ptr<std::ofstream> fp(nullptr);
  FourVector jQ = FourVector(), jB = FourVector(), jS = FourVector();
  constexpr bool compute_gradient = false;
  if (enable_ascii_) {
    fp = output_ascii_files_[ThermodynamicQuantity::j_QBS];
    *fp << std::setprecision(14);
    *fp << std::scientific;
    *fp << ctime << std::endl;
  }
  if (enable_binary_) {
    fp = output_binary_files_[ThermodynamicQuantity::j_QBS];
    assert(sizeof(ctime) == sizeof(double));
    fp->write(reinterpret_cast<char *>(&ctime), sizeof(ctime));
  }
  lattice.iterate_sublattice(
      {0, 0, 0}, dim, [&](DensityOnLattice &, int ix, int iy, int iz) {
        const ThreeVector position = lattice.cell_center(ix, iy, iz);
        jQ.reset();
        jB.reset();
        jS.reset();
        for (const Particles &particles : ensembles) {
          jQ += std::get<1>(current_eckart(
              position, particles, dens_param, DensityType::Charge,
              compute_gradient, out_par_.td_smearing));
          jB += std::get<1>(current_eckart(
              position, particles, dens_param, DensityType::Baryon,
              compute_gradient, out_par_.td_smearing));
          jS += std::get<1>(current_eckart(
              position, particles, dens_param, DensityType::Strangeness,
              compute_gradient, out_par_.td_smearing));
        }
        if (enable_ascii_) {
          *fp << jQ[0];
          for (int l = 1; l < 4; l++) {
            *fp << " " << jQ[l];
          }
          for (int l = 0; l < 4; l++) {
            *fp << " " << jB[l];
          }
          for (int l = 0; l < 4; l++) {
            *fp << " " << jS[l];
          }
          *fp << "\n";
        }
        if (enable_binary_) {
          for (int l = 0; l < 4; l++) {
            result = jQ[l];
            fp->write(reinterpret_cast<char *>(&result), sizeof(double));
          }
          for (int l = 0; l < 4; l++) {
            result = jB[l];
            fp->write(reinterpret_cast<char *>(&result), sizeof(double));
          }
          for (int l = 0; l < 4; l++) {
            result = jS[l];
            fp->write(reinterpret_cast<char *>(&result), sizeof(double));
          }
        }
      });
}

void ThermodynamicLatticeOutput::thermodynamics_lattice_output(
    const ThermodynamicQuantity tq,
    RectangularLattice<EnergyMomentumTensor> &lattice, double ctime) {
  if (!enable_output_) {
    return;
  }
  double result;
  const auto dim = lattice.n_cells();
  std::shared_ptr<std::ofstream> fp(nullptr);
  if (enable_ascii_) {
    switch (tq) {
      case ThermodynamicQuantity::Tmn:
        fp = output_ascii_files_[ThermodynamicQuantity::Tmn];
        break;
      case ThermodynamicQuantity::TmnLandau:
        fp = output_ascii_files_[ThermodynamicQuantity::TmnLandau];
        break;
      case ThermodynamicQuantity::LandauVelocity:
        fp = output_ascii_files_[ThermodynamicQuantity::LandauVelocity];
        break;
      default:
        return;
    }
    *fp << std::setprecision(14);
    *fp << std::scientific;
    *fp << ctime << std::endl;
  }
  if (enable_binary_) {
    switch (tq) {
      case ThermodynamicQuantity::Tmn:
        fp = output_binary_files_[ThermodynamicQuantity::Tmn];
        break;
      case ThermodynamicQuantity::TmnLandau:
        fp = output_binary_files_[ThermodynamicQuantity::TmnLandau];
        break;
      case ThermodynamicQuantity::LandauVelocity:
        fp = output_binary_files_[ThermodynamicQuantity::LandauVelocity];
        break;
      default:
        return;
    }
    assert(sizeof(ctime) == sizeof(double));
    fp->write(reinterpret_cast<char *>(&ctime), sizeof(double));
  }
  switch (tq) {
    case ThermodynamicQuantity::Tmn:
      for (int i = 0; i < 4; i++) {
        for (int j = i; j < 4; j++) {
          lattice.iterate_sublattice(
              {0, 0, 0}, dim,
              [&](EnergyMomentumTensor &node, int ix, int, int) {
                if (enable_ascii_) {
                  *fp << node[EnergyMomentumTensor::tmn_index(i, j)] << " ";
                  if (ix == dim[0] - 1) {
                    *fp << "\n";
                  }
                }
                if (enable_binary_) {
                  result = node[EnergyMomentumTensor::tmn_index(i, j)];
                  fp->write(reinterpret_cast<char *>(&result), sizeof(double));
                }
              });
        }
      }
      break;
    case ThermodynamicQuantity::TmnLandau:
      for (int i = 0; i < 4; i++) {
        for (int j = i; j < 4; j++) {
          lattice.iterate_sublattice(
              {0, 0, 0}, dim,
              [&](EnergyMomentumTensor &node, int ix, int, int) {
                if (enable_ascii_) {
                  const FourVector u = node.landau_frame_4velocity();
                  const EnergyMomentumTensor Tmn_L = node.boosted(u);
                  *fp << Tmn_L[EnergyMomentumTensor::tmn_index(i, j)] << " ";
                  if (ix == dim[0] - 1) {
                    *fp << "\n";
                  }
                }
                if (enable_binary_) {
                  const FourVector u = node.landau_frame_4velocity();
                  const EnergyMomentumTensor Tmn_L = node.boosted(u);
                  result = Tmn_L[EnergyMomentumTensor::tmn_index(i, j)];
                  fp->write(reinterpret_cast<char *>(&result), sizeof(double));
                }
              });
        }
      }
      break;
    case ThermodynamicQuantity::LandauVelocity:
      lattice.iterate_sublattice(
          {0, 0, 0}, dim, [&](EnergyMomentumTensor &node, int, int, int) {
            if (enable_ascii_) {
              const FourVector u = node.landau_frame_4velocity();
              const ThreeVector v = -u.velocity();
              *fp << v.x1() << " " << v.x2() << " " << v.x3() << "\n";
            }
            if (enable_binary_) {
              const FourVector u = node.landau_frame_4velocity();
              ThreeVector v = -u.velocity();
              fp->write(reinterpret_cast<char *>(&v), 3 * sizeof(double));
            }
          });
      break;
    default:
      return;
  }
}

int ThermodynamicLatticeOutput::to_int(const ThermodynamicQuantity &tq) {
  switch (tq) {
    case ThermodynamicQuantity::EckartDensity:
      return 0;
    case ThermodynamicQuantity::Tmn:
      return 1;
    case ThermodynamicQuantity::TmnLandau:
      return 2;
    case ThermodynamicQuantity::LandauVelocity:
      return 3;
    case ThermodynamicQuantity::j_QBS:
      return 4;
    default:
      throw std::runtime_error(
          "Error when converting a thermodynamic quantity "
          "to an int, unknown quantity.");
  }
}

std::string ThermodynamicLatticeOutput::make_filename(const std::string &descr,
                                                      const int event_number,
                                                      const char type) {
  char suffix[13];
  assert((type == 'a') || (type == 'b'));
  if (type == 'a') {
    snprintf(suffix, sizeof(suffix), "_%07i.dat", event_number);
  } else {
    snprintf(suffix, sizeof(suffix), "_%07i.bin", event_number);
  }
  return base_path_.string() + std::string("/") + descr + std::string(suffix);
}

std::string ThermodynamicLatticeOutput::make_varname(
    const ThermodynamicQuantity tq, const DensityType dens_type) {
  return std::string(to_string(dens_type)) + std::string("_") +
         std::string(to_string(tq));
}

void ThermodynamicLatticeOutput::write_therm_lattice_ascii_header(
    std::shared_ptr<std::ofstream> fp, const ThermodynamicQuantity &tq) {
  *fp << std::setprecision(2);
  *fp << std::fixed;
  *fp << "#Thermodynamic Lattice Output version: "
      << ThermodynamicLatticeOutput::version << std::endl;
  *fp << std::setprecision(6);
  *fp << "#Quantity:"
      << " " << std::string(to_string(tq)) << std::endl;
  *fp << "#Grid dimensions: " << nodes_[0] << " " << nodes_[1] << " "
      << nodes_[2] << std::endl;
  *fp << "#Grid spacing: " << sizes_[0] << " " << sizes_[1] << " " << sizes_[2]
      << std::endl;
  *fp << "#Grid origin: " << origin_[0] << " " << origin_[1] << " "
      << origin_[2] << std::endl;
}

void ThermodynamicLatticeOutput::write_therm_lattice_binary_header(
    std::shared_ptr<std::ofstream> fp, const ThermodynamicQuantity &tq) {
  auto variable_id = to_int(tq);
  fp->write(
      reinterpret_cast<const char *>(&ThermodynamicLatticeOutput::version),
      sizeof(double));
  fp->write(reinterpret_cast<char *>(&variable_id), sizeof(int));
  fp->write(reinterpret_cast<char *>(&nodes_), sizeof(nodes_));
  fp->write(reinterpret_cast<char *>(&sizes_), sizeof(sizes_));
  fp->write(reinterpret_cast<char *>(&origin_), sizeof(origin_));
}
}  // namespace smash
