/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/densityoutput.h"

#include <boost/filesystem.hpp>
#include <fstream>
#include <memory>
#include <include/config.h>

#include "include/clock.h"
#include "include/density.h"
#include "include/experimentparameters.h"
#include "include/filedeleter.h"
#include "include/forwarddeclarations.h"
#include "include/particles.h"

namespace Smash {

DensityOutput::DensityOutput(bf::path path, Configuration &&config)
    : file_{std::fopen((path / ("density_out.dat")).native().c_str(), "w")},
      r_(ThreeVector(config.take({"x"}),
                     config.take({"y"}),
                     config.take({"z"}))) {
  fprintf(file_.get(), "# %s density output\n", VERSION_MAJOR);
  fprintf(file_.get(), "# time[fm/c] density[fm^-3] @ (%6.2f, %6.2f, %6.2f)\n",
                               r_.x1(), r_.x2(), r_.x3());
}

DensityOutput::~DensityOutput() {
}

void DensityOutput::at_eventstart(const Particles &/*particles*/,
                               const int event_number) {
  fprintf(file_.get(), "# event %i\n", event_number);
}

void DensityOutput::at_eventend(const Particles &/*particles*/,
                               const int /*event_number*/) {
  std::fflush(file_.get());
}

void DensityOutput::at_intermediate_time(const Particles &/*particles*/,
                                      const int /*event_number*/,
                                      const Clock &/*t*/) {
}

void DensityOutput::thermodynamics_output(const Particles &particles,
                                          const ExperimentParameters &param) {
 const ParticleList plist = ParticleList(particles.data().begin(),
                                         particles.data().end());
 const double rho = four_current(r_, plist, param.gaussian_sigma, baryon_density,
                                 param.testparticles).abs();
 fprintf(file_.get(), "%g %g\n", param.labclock.current_time(), rho);
}

void DensityOutput::density_along_line(const char * file_name,
                        const ParticleList &plist,
                        double gs_sigma, Density_type dens_type, int ntest,
                        const ThreeVector &line_start,
                        const ThreeVector &line_end, int n_points) {
  ThreeVector r;
  double rho_eck;
  std::ofstream a_file;
  a_file.open(file_name, std::ios::out);

  for (int i = 0; i <= n_points; i++) {
    r = line_start + (line_end - line_start) * (1.0 * i / n_points);
    rho_eck = four_current(r, plist, gs_sigma, dens_type, ntest).abs();
    a_file << r.x1() << " " <<
              r.x2() << " " <<
              r.x3() << " " << rho_eck << "\n";
  }
  a_file.close();
}

}  // namespace Smash
