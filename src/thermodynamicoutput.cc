/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/thermodynamicoutput.h"

#include <fstream>
#include <memory>

#include <boost/filesystem.hpp>

#include "include/clock.h"
#include "include/config.h"
#include "include/density.h"
#include "include/energymomentumtensor.h"
#include "include/experimentparameters.h"
#include "include/filedeleter.h"
#include "include/forwarddeclarations.h"
#include "include/particles.h"
#include "include/vtkoutput.h"

namespace Smash {

ThermodynamicOutput::ThermodynamicOutput(const bf::path &path,
                                         Configuration &&config)
    : file_{std::fopen((path / ("thermodynamics.dat")).native().c_str(), "w")},
      td_set_(config.take({"Quantities"}).convert_for(td_set_)),
      dens_type_(config.take({"Type"})),
      smearing_(config.take({"Smearing"}, true)) {
  const std::array<double, 3> a = config.take({"Position"});
  r_ = ThreeVector(a[0], a[1], a[2]);
  std::fprintf(file_.get(), "# %s thermodynamics output\n", VERSION_MAJOR);
  std::fprintf(file_.get(), "# @ point (%6.2f, %6.2f, %6.2f) [fm]\n",
                                      r_.x1(), r_.x2(), r_.x3());
  std::fprintf(file_.get(), "# %s\n",  to_string(dens_type_));
  std::fprintf(file_.get(), "# time [fm/c], ");
  if (td_set_.count(ThermodynamicQuantity::EckartDensity) > 0) {
    std::fprintf(file_.get(), "%s [fm^-3], ",
       to_string(ThermodynamicQuantity::EckartDensity));
  }
  if (td_set_.count(ThermodynamicQuantity::Tmn) > 0) {
    std::fprintf(file_.get(), "%s [GeV/fm^-3] 00 01 02 03 11 12 13 22 23 33, ",
       to_string(ThermodynamicQuantity::Tmn));
  }
  if (td_set_.count(ThermodynamicQuantity::TmnLandau) > 0) {
    std::fprintf(file_.get(), "%s [GeV/fm^-3] 00 01 02 03 11 12 13 22 23 33, ",
       to_string(ThermodynamicQuantity::TmnLandau));
  }
  if (td_set_.count(ThermodynamicQuantity::LandauVelocity) > 0) {
    std::fprintf(file_.get(), "%s x y z ",
       to_string(ThermodynamicQuantity::LandauVelocity));
  }
  std::fprintf(file_.get(), "\n");
}

ThermodynamicOutput::~ThermodynamicOutput() {
}

void ThermodynamicOutput::at_eventstart(const Particles &/*particles*/,
                               const int event_number) {
  std::fprintf(file_.get(), "# event %i\n", event_number);
}

void ThermodynamicOutput::at_eventend(const Particles &/*particles*/,
                               const int /*event_number*/) {
  std::fflush(file_.get());
}

void ThermodynamicOutput::at_intermediate_time(const Particles &particles,
                                         const Clock &clock,
                                         const DensityParameters &dens_param) {
  std::fprintf(file_.get(), "%6.2f ", clock.current_time());
  constexpr bool compute_gradient = false;
  if (td_set_.count(ThermodynamicQuantity::EckartDensity) > 0) {
    const double rho = rho_eckart(r_, particles, dens_param, dens_type_,
                                compute_gradient).first;
    std::fprintf(file_.get(), "%7.4f ", rho);
  }
  if (td_set_.count(ThermodynamicQuantity::Tmn) > 0 ||
      td_set_.count(ThermodynamicQuantity::TmnLandau) > 0 ||
      td_set_.count(ThermodynamicQuantity::LandauVelocity) > 0) {
    EnergyMomentumTensor Tmn;
    for (const auto &p : particles) {
      const double dens_factor = density_factor(p.type(), dens_type_);
      if (std::abs(dens_factor) < really_small) {
        continue;
      }
      if (smearing_) {
        const auto sf = unnormalized_smearing_factor(
                              p.position().threevec() -r_,
                              p.momentum(), 1.0/p.momentum().abs(),
                              dens_param, compute_gradient).first;
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
    if (td_set_.count(ThermodynamicQuantity::Tmn) > 0) {
      for (int i = 0; i < 10; i++) {
        std::fprintf(file_.get(), "%15.12f ", Tmn[i]);
      }
    }
    if (td_set_.count(ThermodynamicQuantity::TmnLandau) > 0) {
      for (int i = 0; i < 10; i++) {
        std::fprintf(file_.get(), "%7.4f ", Tmn_L[i]);
      }
    }
    if (td_set_.count(ThermodynamicQuantity::LandauVelocity) > 0) {
      std::fprintf(file_.get(), "%7.4f %7.4f %7.4f",
                              -u[1]/u[0], -u[2]/u[0], -u[3]/u[0]);
    }
  }
  std::fprintf(file_.get(), "\n");
}

void ThermodynamicOutput::density_along_line(const char * file_name,
                        const ParticleList &plist,
                        const DensityParameters &param,
                         DensityType dens_type,
                        const ThreeVector &line_start,
                        const ThreeVector &line_end, int n_points) {
  ThreeVector r;
  std::ofstream a_file;
  a_file.open(file_name, std::ios::out);
  const bool compute_gradient = false;

  for (int i = 0; i <= n_points; i++) {
    r = line_start + (line_end - line_start) * (1.0 * i / n_points);
    double rho_eck = rho_eckart(r, plist, param,
                                dens_type, compute_gradient).first;
    a_file << r.x1() << " " <<
              r.x2() << " " <<
              r.x3() << " " << rho_eck << "\n";
  }
  a_file.close();
}

}  // namespace Smash
