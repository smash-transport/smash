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
  double tmp;

  for (const auto &p : plist) {
    if (!p.is_baryon()) {
      continue;
    }
    const ThreeVector ri = p.position().threevec();
    // If particle is too far - reject it immediately: its input is too small
    if ((r - ri).sqr() > (6*gs_sigma) * (6*gs_sigma)) {
      continue;
    }

    const ThreeVector betai = p.velocity();
    const double inv_gammai = p.inverse_gamma();

    // Get distance between particle and r in the particle rest frame
    tmp = ((r - ri) * betai) / (inv_gammai * (1. + inv_gammai));
    const ThreeVector dr_rest = r - ri + betai * tmp;

    tmp = 0.5 * dr_rest.sqr() / (gs_sigma * gs_sigma);
    tmp = p.pdgcode().baryon_number() * std::exp(- tmp) / inv_gammai;
    jmu += FourVector(1., betai) * tmp;
  }

  const double norm = twopi * std::sqrt(twopi) * gs_sigma*gs_sigma*gs_sigma;
  jmu /= norm;

  // j^0 = jmu.x0() is computational frame density
  // jmu.abs() = sqrt(j^mu j_mu) is Eckart rest frame density
  return jmu;
}

double ModusDefault::potential(ThreeVector r, const ParticleList &plist,
                                                         double gs_sigma) {
  const double rho_eckart = baryon_jmu(r, plist, gs_sigma).abs();
  return skyrme_a * (rho_eckart/rho0) +
                     skyrme_b * std::pow(rho_eckart/rho0, skyrme_tau);
}

ThreeVector ModusDefault::potential_gradient(ThreeVector r,
                              const ParticleList &plist, double gs_sigma) {
  // baryon four-current in computational frame
  FourVector jbmu(0.0, 0.0, 0.0, 0.0);
  // derivatives of baryon four-current in computational frame
  FourVector djbmu_dx(0.0, 0.0, 0.0, 0.0);
  FourVector djbmu_dy(0.0, 0.0, 0.0, 0.0);
  FourVector djbmu_dz(0.0, 0.0, 0.0, 0.0);
  double tmp1, tmp2;

  for (const auto &p : plist) {
    if (!p.is_baryon()) {
      continue;
    }
    const ThreeVector ri = p.position().threevec();
    // If particle is too far - reject it immediately: its input is too small
    if ((r - ri).sqr() > (6*gs_sigma) * (6*gs_sigma)) {
      continue;
    }

    const ThreeVector betai = p.velocity();
    const double inv_gammai = p.inverse_gamma();

    // Get distance between particle and r in the particle rest frame
    tmp1 = inv_gammai * (1. + inv_gammai);
    const ThreeVector dr_rest = r - ri + betai * (((r - ri) * betai) / tmp1);

    tmp2 = 0.5 * dr_rest.sqr() / (gs_sigma * gs_sigma);
    tmp2 = p.pdgcode().baryon_number() * std::exp(- tmp2) / inv_gammai;
    jbmu += FourVector(1., betai) * tmp2;
    djbmu_dx += FourVector(1., betai) * (tmp2 * dr_rest.x1() * (1.0 + betai.x1() * betai.x1() / tmp1));
    djbmu_dy += FourVector(1., betai) * (tmp2 * dr_rest.x2() * (1.0 + betai.x2() * betai.x2() / tmp1));
    djbmu_dx += FourVector(1., betai) * (tmp2 * dr_rest.x3() * (1.0 + betai.x3() * betai.x3() / tmp1));
  }
  const double norm1 = twopi * std::sqrt(twopi) * gs_sigma*gs_sigma*gs_sigma;
  jbmu /= norm1;
  const double norm2 = - norm1 * gs_sigma*gs_sigma;
  djbmu_dx /= norm2;
  djbmu_dy /= norm2;
  djbmu_dz /= norm2;

  // Eckart rest frame baryon density
  const double rhob = jbmu.abs();

  // Gradient of Eckart rest frame baryon density
  const ThreeVector drho_dr = ThreeVector(jbmu.Dot(djbmu_dx),
                                          jbmu.Dot(djbmu_dy),
                                          jbmu.Dot(djbmu_dz)) / rhob;

  // Derivative of potential with respect to density
  tmp1 = skyrme_tau * std::pow(rhob/rho0, skyrme_tau);
  const double dpotential_drho = (skyrme_a  + skyrme_b * tmp1) / rho0;

  return drho_dr * dpotential_drho;

}


}  // namespace Smash
