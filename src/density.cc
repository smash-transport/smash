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

float density_factor(const PdgCode pdg, DensityType dens_type) {
  switch (dens_type) {
    case DensityType::particle:
      return 1.f;
    case DensityType::baryon:
      return static_cast<float>(pdg.baryon_number());
    case DensityType::baryonic_isospin:
      if (pdg.is_baryon()) {
        return pdg.isospin3_rel();
      } else {
        return 0.f;
      }
    case DensityType::pion:
      {
        const auto pdg_code = pdg.code();
        if (pdg_code == 0x111      // pi0
            || pdg_code == 0x211   // pi+
            || pdg_code == -0x211  // pi-
          ) {
          return 1.f;
        } else {
          return 0.f;
        }
      }
    default:
      return 0.f;
  }
}

std::pair<double, ThreeVector> unnormalized_smearing_factor(
                       const ThreeVector &r, const FourVector &p,
                       const double two_sigma_sqr, const double r_cut_sqr,
                       const bool compute_gradient) {
  const double r_sqr = r.sqr();
  if (r_sqr > r_cut_sqr) {
    // Distance from particle to point of interest > r_cut
    return std::make_pair(0.0, ThreeVector(0.0, 0.0, 0.0));
  }
  const double m_sqr = p.sqr();
  if (m_sqr < really_small) {
    const auto &log = logger<LogArea::Density>();
    log.warn("Gaussian smearing is undefined for momentum ", p);
    return std::make_pair(0.0, ThreeVector(0.0, 0.0, 0.0));
  }

  const FourVector u = p / std::sqrt(m_sqr);
  const double u_r_scalar = r * u.threevec();
  const double r_rest_sqr = r_sqr + u_r_scalar * u_r_scalar;

  if (r_rest_sqr > r_cut_sqr) {
  // Lorentz contracted distance from particle to point of interest > r_cut
    return std::make_pair(0.0, ThreeVector(0.0, 0.0, 0.0));
  }
  const double sf = std::exp(- r_rest_sqr / two_sigma_sqr) * u.x0();
  const ThreeVector sf_grad = compute_gradient ?
            sf * (r + u.threevec() * u_r_scalar) : ThreeVector(0.0, 0.0, 0.0);

  return std::make_pair(sf, sf_grad);
}

template <typename /*ParticlesContainer*/ T>
std::pair<double, ThreeVector> rho_eckart_impl(const ThreeVector &r,
                                   const T &plist, double gs_sigma,
                                   DensityType dens_type,
                                   int ntest, bool compute_gradient) {
  /* In the array first FourVector is jmu and next 3 are d jmu / dr.
   Division into positive and negative charges is necessary to avoid
   problems with the Eckart frame definition. Example of problem:
   get Eckart frame for two identical oppositely flying bunches of
   electrons and positrons. For this case jmu = (0, 0, 0, non-zero),
   so jmu.abs does not exist and Eckart frame is not defined.
   If one takes rho = jmu_pos.abs - jmu_neg.abs, it is still Lorentz-
   invariant and gives the right limit in non-relativistic case, but
   it gives no such problem.
  */
  std::array<FourVector, 4> jmu_pos, jmu_neg;

  for (const auto &p : plist) {
    const float dens_factor = density_factor(p.pdgcode(), dens_type);
    if (std::abs(dens_factor) < really_small) {
      continue;
    }
    const auto sf_and_grad = unnormalized_smearing_factor(
                            p.position().threevec() -r,
                            p.momentum(),
                            2*gs_sigma*gs_sigma, (4*gs_sigma)*(4*gs_sigma),
                            compute_gradient);
    if (sf_and_grad.first < really_small) {
      continue;
    }
    const FourVector tmp = p.momentum() * (dens_factor / p.momentum().x0());
    if (dens_factor > 0.f) {
      jmu_pos[0] += tmp * sf_and_grad.first;
      if (compute_gradient) {
        for (int k = 1; k <= 3; k++) {
          jmu_pos[k] += tmp * sf_and_grad.second[k-1];
        }
      }
    } else {
      jmu_neg[0] += tmp * sf_and_grad.first;
      if (compute_gradient) {
        for (int k = 1; k <= 3; k++) {
          jmu_neg[k] += tmp * sf_and_grad.second[k-1];
        }
      }
    }
  }

  const double rho_eck = (jmu_pos[0].abs() - jmu_neg[0].abs()) /
            smearing_factor_norm(2*gs_sigma*gs_sigma) / ntest;

  ThreeVector rho_eck_grad;
  if (compute_gradient) {
    for (int i = 1; i < 4; i++) {
      rho_eck_grad[i-1] = 0.;
      if (jmu_pos[0].x0() > really_small) {
        rho_eck_grad[i-1] += jmu_pos[i].Dot(jmu_pos[0]) / jmu_pos[0].abs();
      }
      if (jmu_pos[0].x0() < -really_small) {
        rho_eck_grad[i-1] -= jmu_neg[i].Dot(jmu_neg[0]) / jmu_neg[0].abs();
      }
    }
    rho_eck_grad /= smearing_factor_grad_norm(2*gs_sigma*gs_sigma) * ntest;
  }
  return std::make_pair(rho_eck, rho_eck_grad);
}

std::pair<double, ThreeVector> rho_eckart(const ThreeVector &r,
                  const ParticleList &plist,
                  double gs_sigma, DensityType dens_type, int ntest,
                  bool compute_gradient) {
  return rho_eckart_impl(r, plist, gs_sigma, dens_type, ntest,
                         compute_gradient);
}
std::pair<double, ThreeVector> rho_eckart(const ThreeVector &r,
                  const Particles &plist,
                  double gs_sigma, DensityType dens_type, int ntest,
                  bool compute_gradient) {
  return rho_eckart_impl(r, plist, gs_sigma, dens_type, ntest,
                         compute_gradient);
}

void update_density_lattice(DensityLattice* lat,
                            const LatticeUpdate update,
                            const DensityType dens_type,
                            const ExperimentParameters &par,
                            const Particles &particles) {
  // Do not proceed if lattice does not exists/update not required
  if (lat == nullptr || lat->when_update() != update) {
    return;
  }
  // Add particles to lattice jmus
  lat->reset();
  const int ntest = par.testparticles;
  const double sig = par.gaussian_sigma;
  const double two_sig_sqr = 2 * sig * sig;
  const double r_cut = 4 * sig;
  const double r_cut_sqr = r_cut * r_cut;
  for (const auto &part: particles) {
    const float dens_factor = density_factor(part.pdgcode(), dens_type);
    if (std::abs(dens_factor) < really_small) {
      continue;
    }
    const ThreeVector pos = part.position().threevec();
    lat->iterate_in_radius(pos, r_cut,
      [&](DensityOnLattice &node, int ix, int iy, int iz){
        const ThreeVector r = lat->cell_center(ix, iy, iz);
        const double sf = unnormalized_smearing_factor(pos - r,
                               part.momentum(), two_sig_sqr, r_cut_sqr).first;
        if (sf > really_small) {
          /*std::cout << "Adding particle " << part << " to lattice with"
                    << " smearing factor " << sf <<
                    " and density factor " << dens_factor << std::endl;*/
          node.add_particle(part, sf * dens_factor);
        }
      });
  }
  // Compute density from jmus, take care about smearing factor normalization
  for (auto &node: *lat) {
    node.compute_density(smearing_factor_norm(two_sig_sqr) * ntest);
  }
}


std::ostream& operator<<(std::ostream& os, DensityType dens_type) {
  switch (dens_type) {
    case DensityType::particle:
      os << "particle density";
      break;
    case DensityType::baryon:
      os << "baryon density";
      break;
    case DensityType::baryonic_isospin:
      os << "baryonic isospin density";
      break;
    case DensityType::pion:
      os << "pion density";
      break;
    default:
      os.setstate(std::ios_base::failbit);
  }
  return os;
}

}  // namespace Smash
