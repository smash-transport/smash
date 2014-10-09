/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include <cinttypes>
#include <list>

#include "include/modusdefault.h"
#include "include/constants.h"
#include "include/experiment.h"
#include "include/logging.h"
#include "include/threevector.h"

namespace Smash {

// general propagation routine

void ModusDefault::propagate(Particles *particles,
                             const ExperimentParameters &parameters,
                             const OutputsList &) {
  const auto &log = logger<LogArea::ModusDefault>();
  FourVector distance, position;
  for (ParticleData &data : particles->data()) {
    /* propagation for this time step */
    distance = FourVector(0.0,
                          data.velocity() * parameters.timestep_duration());
    log.debug("Particle ", data, " motion: ", distance);
    position = data.position() + distance;
    position.set_x0(parameters.new_particle_time());
    data.set_4position(position);
  }
}

FourVector ModusDefault::baryon_jmu(ThreeVector r,
                                    const ParticleList &plist,
                                    double gs_sigma) {
  FourVector jmu(0.0, 0.0, 0.0, 0.0);
  ThreeVector ri, betai, dr_rest;
  double inv_gammai, hlp, norm;

  for (const auto &p : plist) {
    if (!p.is_baryon()) continue;
    ri = p.position().threevec();
    betai = p.velocity();
    // printf("betai = %12.4f %12.4f %12.4f\n",betai.x1(),betai.x2(),betai.x3());

    inv_gammai = p.inverse_gamma();
    // printf("gamma_inv = %12.4f\n", inv_gammai);

    // Get distance between particle and r in the particle rest frame
    hlp = ((ri - r) * betai) / (inv_gammai * (1. + inv_gammai));
    dr_rest = r - ri + betai * hlp;

    // Calculate the argument of exponential and check if it is too large
    hlp = 0.5 * dr_rest.sqr() / (gs_sigma * gs_sigma);
    if (hlp > 10.) continue;

    hlp = p.pdgcode().baryon_number() * std::exp(- hlp) / inv_gammai;
    jmu += FourVector(1., betai) * hlp;
  }

  norm = twopi * std::sqrt(twopi) * gs_sigma * gs_sigma * gs_sigma;
  jmu /= norm;

  // j^0 = jmu.x0() is computational frame density
  // jmu.abs() = sqrt(j^mu j_mu) is Eckart rest frame density
  return jmu;
}

}  // namespace Smash
