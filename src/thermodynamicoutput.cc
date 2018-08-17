/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/thermodynamicoutput.h"

#include <fstream>
#include <memory>

#include <boost/filesystem.hpp>

#include "smash/clock.h"
#include "smash/config.h"
#include "smash/density.h"
#include "smash/energymomentumtensor.h"
#include "smash/experimentparameters.h"
#include "smash/forwarddeclarations.h"
#include "smash/particles.h"
#include "smash/vtkoutput.h"

namespace smash {

ThermodynamicOutput::ThermodynamicOutput(const bf::path &path,
                                         const std::string &name,
                                         const OutputParameters &out_par)
    : OutputInterface(name),
      file_{path / "thermodynamics.dat", "w"},
      out_par_(out_par) {
  std::fprintf(file_.get(), "# %s thermodynamics output\n", VERSION_MAJOR);
  const ThreeVector r = out_par.td_position;
  std::fprintf(file_.get(), "# @ point (%6.2f, %6.2f, %6.2f) [fm]\n", r.x1(),
               r.x2(), r.x3());
  std::fprintf(file_.get(), "# %s\n", to_string(out_par.td_dens_type));
  std::fprintf(file_.get(), "# time [fm/c], ");
  if (out_par_.td_rho_eckart) {
    std::fprintf(file_.get(), "%s [fm^-3], ",
                 to_string(ThermodynamicQuantity::EckartDensity));
  }
  if (out_par_.td_tmn) {
    std::fprintf(file_.get(), "%s [GeV/fm^-3] 00 01 02 03 11 12 13 22 23 33, ",
                 to_string(ThermodynamicQuantity::Tmn));
  }
  if (out_par_.td_tmn_landau) {
    std::fprintf(file_.get(), "%s [GeV/fm^-3] 00 01 02 03 11 12 13 22 23 33, ",
                 to_string(ThermodynamicQuantity::TmnLandau));
  }
  if (out_par_.td_v_landau) {
    std::fprintf(file_.get(), "%s x y z ",
                 to_string(ThermodynamicQuantity::LandauVelocity));
  }
  std::fprintf(file_.get(), "\n");
}

ThermodynamicOutput::~ThermodynamicOutput() {}

void ThermodynamicOutput::at_eventstart(const Particles & /*particles*/,
                                        const int event_number) {
  std::fprintf(file_.get(), "# event %i\n", event_number);
}

void ThermodynamicOutput::at_eventend(const Particles & /*particles*/,
                                      const int /*event_number*/,
                                      double /*impact_parameter*/) {
  std::fflush(file_.get());
}

void ThermodynamicOutput::at_intermediate_time(
    const Particles &particles, const Clock &clock,
    const DensityParameters &dens_param) {
  std::fprintf(file_.get(), "%6.2f ", clock.current_time());
  constexpr bool compute_gradient = false;
  if (out_par_.td_rho_eckart) {
    const double rho = rho_eckart(out_par_.td_position, particles, dens_param,
                                  out_par_.td_dens_type, compute_gradient)
                           .first;
    std::fprintf(file_.get(), "%7.4f ", rho);
  }
  if (out_par_.td_tmn || out_par_.td_tmn_landau || out_par_.td_v_landau) {
    EnergyMomentumTensor Tmn;
    for (const auto &p : particles) {
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
      std::fprintf(file_.get(), "%7.4f %7.4f %7.4f", -u[1] / u[0], -u[2] / u[0],
                   -u[3] / u[0]);
    }
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

  for (int i = 0; i <= n_points; i++) {
    r = line_start + (line_end - line_start) * (1.0 * i / n_points);
    double rho_eck =
        rho_eckart(r, plist, param, dens_type, compute_gradient).first;
    a_file << r.x1() << " " << r.x2() << " " << r.x3() << " " << rho_eck
           << "\n";
  }
}

}  // namespace smash
