/*
 *
 *    Copyright (c) 2014-2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/thermodynamicoutput.h"

#include <filesystem>
#include <fstream>
#include <memory>

#include "smash/clock.h"
#include "smash/config.h"
#include "smash/density.h"
#include "smash/energymomentumtensor.h"
#include "smash/experimentparameters.h"
#include "smash/vtkoutput.h"

namespace smash {

/*!\Userguide
 * \page thermodyn_output_user_guide_ ASCII Thermodynamics Output
 *
 * The thermodynamics output (thermodynamics.dat) is used to output either:
 * \li thermodynamic quantities at one point (\key Smearing: on)
 * \li thermodynamic quantities averaged over every particle
 *     in an event (\key Smearing: off - mostly useful for box setups).
 *
 * Smearing can be switched on or off in the configuration file, see
 * \ref input_output_content_specific_ "content-specific output options".
 *
 * The calculated quantities of the file can include: the Eckart density, the
 * energy-momentum tensor in the lab and/or Landau frame, the Landau
 * velocity and the electric/baryonic/strange currents (they will also be
 * outputted in this order); note that the
 * Eckart density is currently not affected by turning smearing off, and
 * will always return the density at a point.
 * \n
 * Which of these quantities are outputted needs to be specified in the config
 * file. They will be calculated using **one** of the following density types:
 * hadron, baryon, baryonic isospin, pion, or none.\n
 * The output file does contain any information about whether all hadrons are
 * included or \ref key_output_thermo_only_part_ "only participants".
 *
 * The format of the file is the following: \n
 *
 * \n
 * **Header in case of evaluation at one point**
 * \code
 * # **smash_version** thermodynamics output
 * # @ point ( **position** ) [fm]
 * # **density_type**
 * # time [fm], ** a list of all columns that will be printed, + units **
 * \endcode
 *
 * The header consists of 4 lines starting with a '#', containing the following
 * information:
 * -# SMASH-version and the information, that 'thermodynamic output' is provided
 * -# The point of evaluation
 * -# The density type, as specified in the config file
 * -# The header with all column names
 *
 * \n
 * **Header in case of average over the entire volume**
 * \code
 * # **smash_version** thermodynamics output
 * # averaged over the entire volume
 * # **density_type**
 * # time [fm], ** a list of all columns that will be printed, + units **
 * \endcode
 *
 * The header consists of 4 lines starting with a '#', containing the following
 * information:
 * -# SMASH-version and the information, that 'thermodynamic output' is provided
 * -# The info that the quantities are 'averaged over the entire volume'
 * -# The density type, as specified in the config file
 * -# The header with all column names
 *
 * \n
 * **Event Header** \n
 * Each event is indicated with an event starting line:
 * \code
 * # event number
 * \endcode
 * where
 * \li \key number: Event number
 *
 * Note, that 'event' is not a variable but a word that is printed. \n
 *
 * The event indication line is followed by the data lines formatted as:
 * <div class="fragment">
 * <div class="line"> <span class="preprocessor">
 * time [density] [10 cols Tmunu_Lab] [10 cols Tmunu_Landau] [3 cols v_Landau]
 * [4 cols el_current] [4 cols bar_current] [4 cols str_current] </span></div>
 * </div>
 * where
 * \li \key density: The density specified in the configuration file.
 * \li \key Tmunu_Lab: Energy-momentum tensor in the lab frame (10 columns).
 * \li \key Tmunu_Landau: Energy-momentum tensor in the Landau frame
 *          (10 columns).
 * \li \key v_Landau: The velocity in Landau frame (3 columns).
 * \li \key el_current: The electric current in the lab frame (4 columns).
 * \li \key bar_current: The baryonic current in the lab frame (4 columns).
 * \li \key str_current: The strange current in the lab frame (4 columns).
 *
 * Note that the number of columns depends on what was specified in the
 * configuration file,
 * i.e. all quantities in brackets will only be there if specifically asked for.
 */

ThermodynamicOutput::ThermodynamicOutput(const std::filesystem::path &path,
                                         const std::string &name,
                                         const OutputParameters &out_par)
    : OutputInterface(name),
      file_{path / "thermodynamics.dat", "w"},
      out_par_(out_par) {
  std::fprintf(file_.get(), "# %s thermodynamics output\n", SMASH_VERSION);
  const ThreeVector r = out_par.td_position;
  if (out_par_.td_only_participants) {
    std::fprintf(file_.get(), "# only participants are taken into account\n");
  }
  if (out_par_.td_smearing) {
    std::fprintf(file_.get(), "# @ point (%6.2f, %6.2f, %6.2f) [fm]\n", r.x1(),
                 r.x2(), r.x3());
  } else {
    std::fprintf(file_.get(), "# averaged over the entire volume\n");
  }
  std::fprintf(file_.get(), "# %s\n", to_string(out_par.td_dens_type));
  std::fprintf(file_.get(), "# time [fm], ");
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

ThermodynamicOutput::~ThermodynamicOutput() {}

void ThermodynamicOutput::at_eventstart(
    const std::vector<Particles> & /*particles*/, const int event_number) {
  std::fprintf(file_.get(), "# event %i\n", event_number);
}

void ThermodynamicOutput::at_eventend(
    const std::vector<Particles> & /*particles*/, const int /*event_number*/) {
  std::fflush(file_.get());
}

void ThermodynamicOutput::at_intermediate_time(
    const std::vector<Particles> &ensembles,
    const std::unique_ptr<Clock> &clock, const DensityParameters &dens_param) {
  std::fprintf(file_.get(), "%6.2f ", clock->current_time());
  constexpr bool compute_gradient = false;
  if (out_par_.td_rho_eckart) {
    FourVector jmu = FourVector();
    for (const Particles &particles : ensembles) {
      jmu += std::get<1>(current_eckart(
          out_par_.td_position, particles, dens_param, out_par_.td_dens_type,
          compute_gradient, out_par_.td_smearing));
    }
    std::fprintf(file_.get(), "%15.12f ", jmu.abs());
  }
  if (out_par_.td_tmn || out_par_.td_tmn_landau || out_par_.td_v_landau) {
    EnergyMomentumTensor Tmn;
    for (const Particles &particles : ensembles) {
      for (const auto &p : particles) {
        if (dens_param.only_participants()) {
          // if this condition holds, the hadron is a spectator and we skip it
          if (p.get_history().collisions_per_particle == 0) {
            continue;
          }
        }
        const double dens_factor =
            density_factor(p.type(), out_par_.td_dens_type);
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
        std::fprintf(file_.get(), "%15.12f ", Tmn_L[i]);
      }
    }
    if (out_par_.td_v_landau) {
      std::fprintf(file_.get(), "%15.12f %15.12f %15.12f ", -u[1] / u[0],
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
    std::fprintf(file_.get(), "%15.12f %15.12f %15.12f %15.12f ", jQ[0], jQ[1],
                 jQ[2], jQ[3]);
    std::fprintf(file_.get(), "%15.12f %15.12f %15.12f %15.12f ", jB[0], jB[1],
                 jB[2], jB[3]);
    std::fprintf(file_.get(), "%15.12f %15.12f %15.12f %15.12f ", jS[0], jS[1],
                 jS[2], jS[3]);
  }
  std::fprintf(file_.get(), "\n");
}

void ThermodynamicOutput::density_along_line(
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
