/*
 *
 *    Copyright (c) 2014-2021
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/thermodynamiclatticeoutput.h"
#include "smash/thermodynamicoutput.h"

#include <fstream>
#include <memory>

#include <boost/filesystem.hpp>

#include "smash/clock.h"
#include "smash/config.h"
#include "smash/density.h"
#include "smash/energymomentumtensor.h"
#include "smash/experimentparameters.h"
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
 * - variable id (int) : see the correspondence list below
 * - id of the density type (int ) : see the correspondence list below 
 * - nx, ny, nz (int) :  cells of the lattice along x, y, z, respectively (set in: Lattice->Cell_Number)
 * - 0, 0, 0 (int) : three zeroes for padding purposes and free slot for future uses
 * - t (double) : time 
 * - x0, y0, z0 (double) : coordinates of the origin of the lattice (set in Lattice->Origin)
 * - dx, dy, dz (double) : size of the lattice (set in Lattice->Sizes)
 * - smash version number (double)
 * - the data payload
 *
 * The data payload consists in the values of the quantity in the following order:
 * <\code>
 * for (h=0;h<{number of timesteps};h++) {
 *   for (i=0;i<nx;i++) {
 *     for (j=0;j<ny;j++) {
 *       for (k=0;k<nz;k++) {
 *         for (l=0;l<{number of quantity components};l++) {
 *           quantity (double)
 * <\endcode>          
 *
 * Here is the list of ids of the quantities and their number of components
 * 
 * \li \key density: The density specified in the configuration file - variable id: 0, components: 1
 * \li \key Tmunu_Lab: Energy-momentum tensor in the lab frame - variable id: 1, components: 10
 * \li \key Tmunu_Landau: Energy-momentum tensor in the Landau frame - variable id: 2, components: 10
 * \li \key v_Landau: The velocity in Landau frame - variable id: 3, components: 3
 * \li \key el_current: The electric current in the lab frame - variable id: 4, components: 4
 * \li \key bar_current: The baryonic current in the lab frame - variable id: 5, components: 4
 * \li \key str_current: The strange current in the lab frame - variable id: 6, components: 4
 *
 * The id of the density types are:
 *
 * \li "baryon": 0
 * \li "hadron": 1
 * \li "baryonic isospin": 2
 * \li "pion": 3
 * \li "total isospin": 4
 * \li "none": 5
 */

ThermodynamicLatticeOutput::ThermodynamicLatticeOutput(const bf::path &path,
                                         const std::string &name,
                                         const OutputParameters &out_par)
    : OutputInterface(name),
      file_{path / "thermodynamics.dat", "w"},
      out_par_(out_par) {
  std::fprintf(file_.get(), "# %s thermodynamics output\n", VERSION_MAJOR);
  const ThreeVector r = out_par.td_position;
  if (out_par_.td_smearing) {
    std::fprintf(file_.get(), "# @ point (%6.2f, %6.2f, %6.2f) [fm]\n", r.x1(),
                 r.x2(), r.x3());
  } else {
    std::fprintf(file_.get(), "# averaged over the entire volume\n");
  }
  std::fprintf(file_.get(), "# %s\n", to_string(out_par.td_dens_type));
  std::fprintf(file_.get(), "# time [fm/c], ");
  if (out_par_.td_rho_eckart) {
    std::fprintf(file_.get(), "%s [fm^-3], ",
                 to_string(ThermodynamicQuantity::EckartDensity));
  }
  if (out_par_.td_tmn) {
    if (out_par_.td_smearing) {
      std::fprintf(file_.get(), "%s [GeV/fm^3] 00 01 02 03 11 12 13 22 23 33, ",
                   to_string(ThermodynamicQuantity::Tmn));
    } else {
      std::fprintf(file_.get(), "%s [GeV] 00 01 02 03 11 12 13 22 23 33, ",
                   to_string(ThermodynamicQuantity::Tmn));
    }
  }
  if (out_par_.td_tmn_landau) {
    if (out_par_.td_smearing) {
      std::fprintf(file_.get(), "%s [GeV/fm^3] 00 01 02 03 11 12 13 22 23 33, ",
                   to_string(ThermodynamicQuantity::TmnLandau));
    } else {
      std::fprintf(file_.get(), "%s [GeV] 00 01 02 03 11 12 13 22 23 33, ",
                   to_string(ThermodynamicQuantity::TmnLandau));
    }
  }
  if (out_par_.td_v_landau) {
    std::fprintf(file_.get(), "%s x y z, ",
                 to_string(ThermodynamicQuantity::LandauVelocity));
  }
  if (out_par_.td_jQBS) {
    if (out_par_.td_smearing) {
      std::fprintf(file_.get(), "j_QBS [(Q,B,S)/fm^3] (0 1 2 3)x3");
    } else {
      std::fprintf(file_.get(), "j_QBS [(Q,B,S)] (0 1 2 3)x3");
    }
  }
  std::fprintf(file_.get(), "\n");
}

ThermodynamicLatticeOutput::~ThermodynamicLatticeOutput() {}

void ThermodynamicLatticeOutput::at_eventstart(
    const std::vector<Particles> & /*particles*/, const int event_number) {
  std::fprintf(file_.get(), "# event %i\n", event_number);
}

void ThermodynamicLatticeOutput::at_eventend(
    const std::vector<Particles> & /*particles*/, const int /*event_number*/) {
  std::fflush(file_.get());
}

void ThermodynamicLatticeOutput::at_intermediate_time(
    const std::vector<Particles> &ensembles,
    const std::unique_ptr<Clock> &clock, const DensityParameters &dens_param) {
  const double n_ensembles = ensembles.size();
  std::fprintf(file_.get(), "%6.2f ", clock->current_time());
  constexpr bool compute_gradient = false;
  if (out_par_.td_rho_eckart) {
    FourVector jmu = FourVector();
    for (const Particles &particles : ensembles) {
      jmu += std::get<1>(current_eckart(
          out_par_.td_position, particles, dens_param, out_par_.td_dens_type,
          compute_gradient, out_par_.td_smearing));
    }
    std::fprintf(file_.get(), "%7.4f ", jmu.abs() / n_ensembles);
  }
  if (out_par_.td_tmn || out_par_.td_tmn_landau || out_par_.td_v_landau) {
    EnergyMomentumTensor Tmn;
    for (const Particles &particles : ensembles) {
      for (const auto &p : particles) {
        const double dens_factor =
            density_factor(p.type(), out_par_.td_dens_type) / n_ensembles;
        if (std::abs(dens_factor) < really_small) {
          continue;
        }
        if (out_par_.td_smearing) {
          const auto sf =
              unnormalized_smearing_factor(
                  p.position().threevec() - out_par_.td_position, p.momentum(),
                  1.0 / p.momentum().abs(), dens_param, compute_gradient)
                  .first;
          if (sf < really_small) {
            continue;
          }
          Tmn.add_particle(p, dens_factor * sf * dens_param.norm_factor_sf());
        } else {
          Tmn.add_particle(p, dens_factor);
        }
      }
    }
    const FourVector u = Tmn.landau_frame_4velocity();
    const EnergyMomentumTensor Tmn_L = Tmn.boosted(u);
    if (out_par_.td_tmn) {
      for (int i = 0; i < 10; i++) {
        std::fprintf(file_.get(), "%15.12f ", Tmn[i]);
      }
    }
    if (out_par_.td_tmn_landau) {
      for (int i = 0; i < 10; i++) {
        std::fprintf(file_.get(), "%7.4f ", Tmn_L[i]);
      }
    }
    if (out_par_.td_v_landau) {
      std::fprintf(file_.get(), "%7.4f %7.4f %7.4f ", -u[1] / u[0],
                   -u[2] / u[0], -u[3] / u[0]);
    }
  }
  if (out_par_.td_jQBS) {
    FourVector jQ = FourVector(), jB = FourVector(), jS = FourVector();
    for (const Particles &particles : ensembles) {
      jQ += std::get<1>(current_eckart(out_par_.td_position, particles,
                                       dens_param, DensityType::Charge,
                                       compute_gradient, out_par_.td_smearing));
      jB += std::get<1>(current_eckart(out_par_.td_position, particles,
                                       dens_param, DensityType::Baryon,
                                       compute_gradient, out_par_.td_smearing));
      jS += std::get<1>(current_eckart(out_par_.td_position, particles,
                                       dens_param, DensityType::Strangeness,
                                       compute_gradient, out_par_.td_smearing));
    }
    jQ /= n_ensembles;
    jS /= n_ensembles;
    jB /= n_ensembles;
    std::fprintf(file_.get(), "%15.12f %15.12f %15.12f %15.12f ", jQ[0], jQ[1],
                 jQ[2], jQ[3]);
    std::fprintf(file_.get(), "%15.12f %15.12f %15.12f %15.12f ", jB[0], jB[1],
                 jB[2], jB[3]);
    std::fprintf(file_.get(), "%15.12f %15.12f %15.12f %15.12f ", jS[0], jS[1],
                 jS[2], jS[3]);
  }
  std::fprintf(file_.get(), "\n");
}

void ThermodynamicLatticeOutput::density_along_line(
    const char *file_name, const ParticleList &plist,
    const DensityParameters &param, DensityType dens_type,
    const ThreeVector &line_start, const ThreeVector &line_end, int n_points) {
  ThreeVector r;
  std::ofstream a_file;
  a_file.open(file_name, std::ios::out);
  const bool compute_gradient = false;
  const bool smearing = true;

  for (int i = 0; i <= n_points; i++) {
    r = line_start + (line_end - line_start) * (1.0 * i / n_points);
    double rho_eck = std::get<0>(
        current_eckart(r, plist, param, dens_type, compute_gradient, smearing));
    a_file << r.x1() << " " << r.x2() << " " << r.x3() << " " << rho_eck
           << "\n";
  }
}

}  // namespace smash
