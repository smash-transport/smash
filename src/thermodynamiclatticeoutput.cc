/*
 *
 *    Copyright (c) 2014-2021
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/thermodynamiclatticeoutput.h"

#include <fstream>
#include <memory>

#include <boost/filesystem.hpp>

#include "smash/clock.h"
#include "smash/config.h"
#include "smash/density.h"
#include "smash/energymomentumtensor.h"
#include "smash/experimentparameters.h"
#include "smash/thermodynamicoutput.h"
#include "smash/vtkoutput.h"

namespace smash {

/*!\Userguide
 * \page thermodyn_lattice_output_ Thermodynamics Lattice Output
 *
 * The thermodynamics lattice output prints "smeared" thermodynamic quantities
 *  evaluated at the nodes of a Lattice, defined as defined in \ref input_lattice_.
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
 * The format of the file is the following: \n
 *
 * - string containing the version of the output
 * - nx, ny, nz (int) :  cells of the lattice along x, y, z, respectively (set in: Lattice->Cell_Number)
 * - 0, 0, 0 (int) : three zeroes for padding purposes and free slots for future uses
 * - x0, y0, z0 (double) : coordinates of the origin of the lattice (set in Lattice->Origin)
 * - dx, dy, dz (double) : size of the lattice (set in Lattice->Sizes)
 * - the data payload
 *
 * In the case of the energy-momentum tensor, the data payload consists
 * in the values of the quantity in the following order:
 *
 * <\code>
 * for (h=0;h<{number of timesteps};h++) {
 *   output time
 *   for (l=0;l<{number of quantity components};l++) {
 *    for (k=0;k<nz;k++) {
 *     for (j=0;j<ny;j++) {
 *      for (i=0;i<nx;i++) {
 *           quantity value
 *      }
 *      \newline
 *     }
 *    }
 *   }
 * }           
 * <\endcode>          
 *
 * In the case of densities, the data payload consists
 * in the values of the quantity in the following order:
 *
 * <\code>
 * for (h=0;h<{number of timesteps};h++) {
 *   output time
 *   for (k=0;k<nz;k++) {
 *    for (j=0;j<ny;j++) {
 *     for (i=0;i<nx;i++) {
 *         quantity value
 *     }
 *     \newline
 *    }
 *   }
 * }
 * <\endcode>   
 *
 * In the case of Landau velocity the data payload
 * consists in the values of the components
 * Vx, Vy, Vz in the following order:
 *
 * <\code>
 * for (h=0;h<{number of timesteps};h++) {
 *   output time
 *   for (k=0;k<nz;k++) {
 *    for (j=0;j<ny;j++) {
 *     for (i=0;i<nx;i++) {
 *         Vx  Vy  Vz \newline
 *     }
 *    }
 *   }
 * }
 * <\endcode>   
 */


/* initialization of the static member version */
const char* ThermodynamicLatticeOutput::version="v0.5-ASCII";

ThermodynamicLatticeOutput::ThermodynamicLatticeOutput(const bf::path &path,
                                         const std::string &name,
                                         const OutputParameters &out_par)
    : OutputInterface(name),
      base_path_(std::move(path)),
      out_par_(out_par) {
          if (out_par_.td_rho_eckart) {
           output_files_[ThermodynamicQuantity::EckartDensity]=
           std::make_shared<std::ofstream>(nullptr);
      }
          if (out_par_.td_tmn_landau) {
           output_files_[ThermodynamicQuantity::TmnLandau]=
           std::make_shared<std::ofstream>(nullptr);
      }
          if (out_par_.td_tmn) {
           output_files_[ThermodynamicQuantity::Tmn]=
           std::make_shared<std::ofstream>(nullptr);
      }
          if (out_par_.td_v_landau) {
           output_files_[ThermodynamicQuantity::LandauVelocity]=
           std::make_shared<std::ofstream>(nullptr);
      }
          if (out_par_.td_jQBS) {
           std::cout << "Sorry, the Thermodynamic Lattice Output for j_{Q,B,S}"
          << " is not yet implemented.\n";
      }
      }


ThermodynamicLatticeOutput::~ThermodynamicLatticeOutput() {}


void ThermodynamicLatticeOutput::at_eventstart(
    const int event_number, const ThermodynamicQuantity tq,
    const DensityType dens_type,
    RectangularLattice<DensityOnLattice> lattice) {
    // at the next refactoring of the code,
    // this piece should go in the constructor
    const auto dim = lattice.dimensions();
    const auto cs = lattice.cell_sizes();
    const auto orig = lattice.origin();
    nx_ = dim[0];
    ny_ = dim[1];
    nz_ = dim[2];
    dx_ = cs[0];
    dy_ = cs[1];
    dz_ = cs[2];
    x0_ = orig[0];
    y0_ = orig[1];
    z0_ = orig[2];
    std::shared_ptr<std::ofstream> fp(nullptr);
    std::string varname;
    std::string filename;
    varname = make_varname(tq, dens_type);
    filename = make_filename(varname, event_number);
    if (tq==ThermodynamicQuantity::EckartDensity) {
        try {
            output_files_[ThermodynamicQuantity::EckartDensity]->
            open(filename, std::ios::out);
            }
        catch (std::ofstream::failure &e) {
            std::cout << "Error in opening " << filename << std::endl;
        exit(1);
        }
        fp=output_files_[ThermodynamicQuantity::EckartDensity];
    } else {
        std::cout << "There is a problem in the implementation " <<
          "of ThermodynamicLatticeOutput::at_eventstart :\n";
        std::cout << "this specialization expects only " <<
          "tq==ThermodynamicQuantity::EckartDensity, while\n";
        std::cout << "it got " << std::string(to_string(tq)) << " .\n";
        exit(1);
    }
        write_therm_lattice_header(fp, tq);
}

void ThermodynamicLatticeOutput::at_eventstart(
    const int event_number, const ThermodynamicQuantity tq,
    const DensityType dens_type,
    RectangularLattice<EnergyMomentumTensor> lattice ) {
    const auto dim = lattice.dimensions();
    const auto cs = lattice.cell_sizes();
    const auto orig = lattice.origin();
    nx_ = dim[0];
    ny_ = dim[1];
    nz_ = dim[2];
    dx_ = cs[0];
    dy_ = cs[1];
    dz_ = cs[2];
    x0_ = orig[0];
    y0_ = orig[1];
    z0_ = orig[2];

    std::shared_ptr<std::ofstream> fp(nullptr);
    std::string varname;
    std::string filename;
    varname = make_varname(tq, dens_type);
    filename = make_filename(varname, event_number);
    if (tq==ThermodynamicQuantity::TmnLandau) {
        try {
            output_files_[ThermodynamicQuantity::TmnLandau]->
            open(filename, std::ios::out);
        }
        catch (std::ofstream::failure &e) {
            std::cout << "Error in opening " << filename << std::endl;
            exit(1);
        }

        fp=output_files_[ThermodynamicQuantity::TmnLandau];
    } else if (tq==ThermodynamicQuantity::Tmn) {
        try {
            output_files_[ThermodynamicQuantity::Tmn]->
            open(filename, std::ios::out);
        }
        catch (std::ofstream::failure &e) {
            std::cout << "Error in opening " << filename << std::endl;
            exit(1);
        }
        fp=output_files_[ThermodynamicQuantity::Tmn];
    } else {
        try {
            output_files_[ThermodynamicQuantity::LandauVelocity]->
            open(filename, std::ios::out);
        }
        catch (std::ofstream::failure &e) {
            std::cout << "Error in opening " << filename << std::endl;
            exit(1);
        }
        fp=output_files_[ThermodynamicQuantity::LandauVelocity];
    }
        write_therm_lattice_header(fp, tq);
}

void ThermodynamicLatticeOutput::at_eventend(const ThermodynamicQuantity tq) {
    if (tq==ThermodynamicQuantity::EckartDensity) {
        output_files_[ThermodynamicQuantity::EckartDensity]->close();
        return;
    }
    if (tq==ThermodynamicQuantity::TmnLandau) {
        output_files_[ThermodynamicQuantity::TmnLandau]->close();
        return;
    }
    if (tq==ThermodynamicQuantity::Tmn) {
        output_files_[ThermodynamicQuantity::Tmn]->close();
        return;
    }
    if (tq==ThermodynamicQuantity::LandauVelocity) {
        output_files_[ThermodynamicQuantity::LandauVelocity]->close();
        return;
    }
}

void ThermodynamicLatticeOutput::thermodynamics_lattice_output(
    RectangularLattice<DensityOnLattice> &lattice, const double ctime) {
    const auto dim = lattice.dimensions();
    std::shared_ptr<std::ofstream> fp(nullptr);
    fp=output_files_[ThermodynamicQuantity::EckartDensity];
    *fp << std::setprecision(14);
    *fp << std::scientific;
    *fp << ctime << std::endl;
    lattice.iterate_sublattice({0, 0, 0}, dim,
    [&](DensityOnLattice &node, int ix, int, int) {
        *fp << node.density() << " ";
        if (ix == dim[0] - 1) {
            *fp << "\n";
        }
    });
}

void ThermodynamicLatticeOutput::thermodynamics_lattice_output(
    const ThermodynamicQuantity tq,
    RectangularLattice<EnergyMomentumTensor> &lattice, const double ctime) {
    const auto dim = lattice.dimensions();
    std::shared_ptr<std::ofstream> fp(nullptr);
    switch (tq) {
        case ThermodynamicQuantity::Tmn:
            fp=output_files_[ThermodynamicQuantity::Tmn];
            break;
        case ThermodynamicQuantity::TmnLandau:
            fp=output_files_[ThermodynamicQuantity::TmnLandau];
            break;
        case ThermodynamicQuantity::LandauVelocity:
            fp=output_files_[ThermodynamicQuantity::LandauVelocity];
            break;
        default:
            return;
    }
    *fp << std::setprecision(14);
    *fp << std::scientific;
    *fp << ctime << std::endl;
    switch (tq) {
         case ThermodynamicQuantity::Tmn:
             for (int i = 0; i < 4; i++) {
               for (int j = i; j < 4; j++) {
                 lattice.iterate_sublattice({0, 0, 0}, dim,
                 [&](EnergyMomentumTensor &node, int ix, int, int) {
                   *fp << node[EnergyMomentumTensor::tmn_index(i, j)] << " ";
                   if (ix == dim[0] - 1) {
                     *fp << "\n";
                   }
                 });
               }
             }
             break;
         case ThermodynamicQuantity::TmnLandau:
             for (int i = 0; i < 4; i++) {
               for (int j = i; j < 4; j++) {
                 lattice.iterate_sublattice({0, 0, 0}, dim,
                 [&](EnergyMomentumTensor &node, int ix, int, int) {
                   const FourVector u = node.landau_frame_4velocity();
                   const EnergyMomentumTensor Tmn_L = node.boosted(u);
                   *fp << Tmn_L[EnergyMomentumTensor::tmn_index(i, j)] << " ";
                   if (ix == dim[0] - 1) {
                     *fp << "\n";
                   }
                 });
               }
             }
             break;
         case ThermodynamicQuantity::LandauVelocity:
             lattice.iterate_sublattice({0, 0, 0}, dim,
             [&](EnergyMomentumTensor &node, int, int, int) {
               const FourVector u = node.landau_frame_4velocity();
               const ThreeVector v = -u.velocity();
               *fp << v.x1() << " " << v.x2() << " " << v.x3() << "\n";
                });
             break;
         default:
             return;
    }
}

std::string ThermodynamicLatticeOutput::make_filename(const std::string &descr,
    const int event_number) {
    char suffix[13];
    snprintf(suffix, sizeof(suffix), "_%07i.dat", event_number);
    return base_path_.string() + std::string("/") + descr + std::string(suffix);
}

std::string ThermodynamicLatticeOutput::make_varname(const ThermodynamicQuantity
  tq, const DensityType dens_type) {
    return std::string(to_string(dens_type)) + std::string("_") +
    std::string(to_string(tq));
}

void ThermodynamicLatticeOutput::write_therm_lattice_header
  (std::shared_ptr<std::ofstream> fp, const ThermodynamicQuantity &tq) {
    *fp << std::setprecision(5);
    *fp << std::fixed;
    *fp << "#Thermodynamic Lattice Output version: " <<
      ThermodynamicLatticeOutput::version << std::endl;
    *fp << "#Quantity:" << std::string(to_string(tq)) << std::endl;
    *fp << "#Grid information:" << std::endl;
    *fp << "#Dimensions: " << nx_ << " " << ny_ << " " << nz_ << std::endl;
    *fp << "#Spacing: " << dx_ << " " << dy_ << " " << dz_ << std::endl;
    *fp << "#Origin: " << x0_ << " " << y0_ << " " << z0_ << std::endl;
}
}  // namespace smash
