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
 Density_type dens_type = baryon;
 const ParticleList plist = ParticleList(particles.data().begin(),
                                         particles.data().end());
 const double rho = four_current(r_, plist, param.gaussian_sigma, dens_type,
                                 param.testparticles).abs();
 fprintf(file_.get(), "%g %g\n", param.labclock.current_time(), rho);
}
}  // namespace Smash
