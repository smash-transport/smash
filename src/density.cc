/*
 *
 *    Copyright (c) 2013-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/constants.h"
#include "include/density.h"

namespace Smash {

bool particle_in_denstype(const PdgCode pdg, Density_type dens_type) {
  switch (dens_type) {
  case baryon_density:
  case baryonic_isospin_density:
    return pdg.is_baryon();
  default:
    return false;
  }
}

FourVector four_current(const ThreeVector &r, const ParticleList &plist,
                        double gs_sigma, Density_type dens_type, int ntest) {
  FourVector jmu(0.0, 0.0, 0.0, 0.0);
  double tmp;

  for (const auto &p : plist) {
    if (!particle_in_denstype(p.pdgcode(), dens_type)) {
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

    tmp = std::exp(- 0.5 * dr_rest.sqr() / (gs_sigma * gs_sigma)) / inv_gammai;
    switch (dens_type) {
    case baryon_density:
      tmp *= p.pdgcode().baryon_number();
      break;
    case baryonic_isospin_density:
      tmp *= p.pdgcode().isospin3_rel();
      break;
    default:
      break;
    }
    jmu += FourVector(1., betai) * tmp;
  }

  const double norm = twopi * std::sqrt(twopi) * gs_sigma*gs_sigma*gs_sigma;
  jmu /= norm;

  // j^0 = jmu.x0() is computational frame density
  // jmu.abs() = sqrt(j^mu j_mu) is Eckart rest frame density
  return jmu / ntest;
}

std::pair<double, ThreeVector> rho_eckart_gradient(const ThreeVector &r,
                                const ParticleList &plist, double gs_sigma,
                                Density_type dens_type, int ntest) {
  // baryon four-current in computational frame
  FourVector jmu(0.0, 0.0, 0.0, 0.0);
  // derivatives of baryon four-current in computational frame
  FourVector djmu_dx(0.0, 0.0, 0.0, 0.0);
  FourVector djmu_dy(0.0, 0.0, 0.0, 0.0);
  FourVector djmu_dz(0.0, 0.0, 0.0, 0.0);
  double tmp1, tmp2;

  for (const auto &p : plist) {
    if (!particle_in_denstype(p.pdgcode(), dens_type)) {
      continue;
    }
    const ThreeVector ri = p.position().threevec();
    // If particle is too far - reject it immediately: its input is too small
    if ((r - ri).sqr() > (6*gs_sigma) * (6*gs_sigma)) {
      continue;
    }

    const ThreeVector betai = p.velocity();
    // std::cout << "Velocity: " << betai << std::endl;
    const double inv_gammai = p.inverse_gamma();
    // std::cout << "1/gamma: " << inv_gammai << std::endl;

    // Get distance between particle and r in the particle rest frame
    tmp1 = inv_gammai * (1. + inv_gammai);
    const ThreeVector dr_rest = r - ri + betai * (((r - ri) * betai) / tmp1);

    tmp2 = std::exp(- 0.5 * dr_rest.sqr() / (gs_sigma*gs_sigma)) / inv_gammai;
    switch (dens_type) {
    case baryon_density:
      tmp2 *= p.pdgcode().baryon_number();
      break;
    case baryonic_isospin_density:
      tmp2 *= p.pdgcode().isospin3_rel();
      break;
    default:
      break;
    }
    jmu += FourVector(1., betai) * tmp2;

    // Calculate the gradient: d(0.5 * r_rest^2)/d \vec{r}
    const ThreeVector drest2_dr = dr_rest + betai * ((dr_rest * betai) / tmp1);

    djmu_dx += FourVector(1., betai) * (tmp2 * drest2_dr.x1());
    djmu_dy += FourVector(1., betai) * (tmp2 * drest2_dr.x2());
    djmu_dz += FourVector(1., betai) * (tmp2 * drest2_dr.x3());
  }
  const double norm1 = twopi * std::sqrt(twopi) * gs_sigma*gs_sigma*gs_sigma;
  jmu /= norm1;
  const double norm2 = - norm1 * gs_sigma*gs_sigma;
  djmu_dx /= norm2;
  djmu_dy /= norm2;
  djmu_dz /= norm2;

  // Eckart rest frame density
  const double rho = jmu.abs();

  // std::cout << "Dens: " << jmu << " " << djmu_dz << std::endl;

  // Eckart rest frame density and its gradient
  if (fabs(rho) > really_small) {
    return std::make_pair(rho / ntest, ThreeVector(jmu.Dot(djmu_dx),
                                           jmu.Dot(djmu_dy),
                                           jmu.Dot(djmu_dz)) / (rho * ntest));
  } else {
    return std::make_pair(0.0, ThreeVector(0.0, 0.0, 0.0));
  }
}

}  // namespace Smash
