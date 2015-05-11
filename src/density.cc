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
#include "include/logging.h"
#include "include/particles.h"

namespace Smash {

bool particle_in_denstype(const PdgCode pdg, DensityType dens_type) {
  switch (dens_type) {
  case DensityType::baryon:
  case DensityType::baryonic_isospin:
    return pdg.is_baryon();
  default:
    return false;
  }
}

template <typename /*ParticlesContainer*/ T>
static FourVector four_current_impl(const ThreeVector &r, const T &plist,
                                    double gs_sigma, DensityType dens_type,
                                    int ntest) {
  FourVector jmu(0.0, 0.0, 0.0, 0.0);

  for (const auto &p : plist) {
    if (!particle_in_denstype(p.pdgcode(), dens_type)) {
      continue;
    }
    const ThreeVector ri = p.position().threevec();
    // If particle is too far - reject it immediately: its input is too small
    if ((r - ri).sqr() > (4*gs_sigma) * (4*gs_sigma)) {
      continue;
    }

    const ThreeVector betai = p.velocity();
    const double inv_gammai = p.inverse_gamma();


    // Get distance between particle and r in the particle rest frame
    double tmp = ((r - ri) * betai) / (inv_gammai * (1. + inv_gammai));
    const ThreeVector dr_rest = r - ri + betai * tmp;
    const double drr_sqr = dr_rest.sqr();
    if (drr_sqr > (4*gs_sigma) * (4*gs_sigma)) {
      continue;
    }

    tmp = std::exp(- 0.5 * drr_sqr / (gs_sigma * gs_sigma)) / inv_gammai;
    switch (dens_type) {
    case DensityType::baryon:
      tmp *= p.pdgcode().baryon_number();
      break;
    case DensityType::baryonic_isospin:
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
FourVector four_current(const ThreeVector &r, const ParticleList &plist,
                        double gs_sigma, DensityType dens_type, int ntest) {
  return four_current_impl(r, plist, gs_sigma, dens_type, ntest);
}
FourVector four_current(const ThreeVector &r, const Particles &plist,
                        double gs_sigma, DensityType dens_type, int ntest) {
  return four_current_impl(r, plist, gs_sigma, dens_type, ntest);
}

std::pair<double, ThreeVector> rho_eckart_gradient(const ThreeVector &r,
                                const ParticleList &plist, double gs_sigma,
                                DensityType dens_type, int ntest) {
  const auto &log = logger<LogArea::Density>();

  // baryon four-current in computational frame
  FourVector jmu(0.0, 0.0, 0.0, 0.0);
  // derivatives of baryon four-current in computational frame
  FourVector djmu_dx(0.0, 0.0, 0.0, 0.0);
  FourVector djmu_dy(0.0, 0.0, 0.0, 0.0);
  FourVector djmu_dz(0.0, 0.0, 0.0, 0.0);

  for (const auto &p : plist) {
    if (!particle_in_denstype(p.pdgcode(), dens_type)) {
      continue;
    }
    const ThreeVector ri = p.position().threevec();
    // If particle is too far - reject it immediately: its input is too small
    if ((r - ri).sqr() > (4*gs_sigma) * (4*gs_sigma)) {
      continue;
    }

    const ThreeVector betai = p.velocity();
    // std::cout << "Velocity: " << betai << std::endl;
    const double inv_gammai = p.inverse_gamma();
    // std::cout << "1/gamma: " << inv_gammai << std::endl;

    // Get distance between particle and r in the particle rest frame
    double tmp1 = inv_gammai * (1. + inv_gammai);
    const ThreeVector dr_rest = r - ri + betai * (((r - ri) * betai) / tmp1);
    const double drr_sqr = dr_rest.sqr();
    if (drr_sqr > (4*gs_sigma) * (4*gs_sigma)) {
      continue;
    }

    double tmp2 = std::exp(-0.5 * drr_sqr / (gs_sigma*gs_sigma)) / inv_gammai;
    switch (dens_type) {
    case DensityType::baryon:
      tmp2 *= p.pdgcode().baryon_number();
      break;
    case DensityType::baryonic_isospin:
      tmp2 *= p.pdgcode().isospin3_rel();
      break;
    default:
      break;
    }
    log.debug("Summand to jmu: ", FourVector(1., betai) * tmp2, "² = ",
              (FourVector(1., betai) * tmp2).sqr(), "\n      with jmu: ", jmu,
              "² = ", jmu.sqr());
    jmu += FourVector(1., betai) * tmp2;
    if (jmu.sqr() < 0.) {
      log.debug("negative jmu²: ", jmu.sqr(), ", jmu = ", jmu);
    }

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
  const double rho2 = jmu.sqr();
  // Due to numerical reasons it can happen that rho2 is (slighly) smaller than
  // zero, while analytically it should be positive. Negative values can be
  // reached only for numerical reasons and are thus small. Therefore I(oliiny)
  // set rho to zero if rho is negative. This is done under the assumption that
  // discarding small jmu^2 values doesn't introduce a bias.
  const double rho = (rho2 > 0.0) ? std::sqrt(rho2) : 0.0;

  // Eckart rest frame density and its gradient
  if (std::abs(rho) > really_small) {
    return std::make_pair(rho / ntest, ThreeVector(jmu.Dot(djmu_dx),
                                           jmu.Dot(djmu_dy),
                                           jmu.Dot(djmu_dz)) / (rho * ntest));
  } else {
    return std::make_pair(0.0, ThreeVector(0.0, 0.0, 0.0));
  }
}

std::ostream& operator<<(std::ostream& os, DensityType dens_type) {
  switch(dens_type) {
    case DensityType::baryon:
      os << "baryon density";
      break;
    case DensityType::baryonic_isospin:
      os << "baryonic isospin density";
      break;
    default:
      os.setstate(std::ios_base::failbit);
  }
  return os;
}

}  // namespace Smash
