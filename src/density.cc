/*
 *
 *    Copyright (c) 2013-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/density.h"
#include "smash/constants.h"
#include "smash/logging.h"
#include "smash/particles.h"

namespace smash {

double density_factor(const ParticleType &type, DensityType dens_type) {
  switch (dens_type) {
    case DensityType::Hadron:
      return type.is_hadron() ? 1. : 0.;
    case DensityType::Baryon:
      return static_cast<double>(type.baryon_number());
    case DensityType::BaryonicIsospin:
      return type.is_baryon() ? type.isospin3_rel() : 0.;
    case DensityType::Pion:
      return type.pdgcode().is_pion() ? 1. : 0.;
    default:
      return 0.;
  }
}

std::pair<double, ThreeVector> unnormalized_smearing_factor(
    const ThreeVector &r, const FourVector &p, const double m_inv,
    const DensityParameters &dens_par, const bool compute_gradient) {
  const double r_sqr = r.sqr();
  // Distance from particle to point of interest > r_cut
  if (r_sqr > dens_par.r_cut_sqr()) {
    return std::make_pair(0.0, ThreeVector(0.0, 0.0, 0.0));
  }

  const FourVector u = p * m_inv;
  const double u_r_scalar = r * u.threevec();
  const double r_rest_sqr = r_sqr + u_r_scalar * u_r_scalar;

  // Lorentz contracted distance from particle to point of interest > r_cut
  if (r_rest_sqr > dens_par.r_cut_sqr()) {
    return std::make_pair(0.0, ThreeVector(0.0, 0.0, 0.0));
  }
  const double sf = std::exp(-r_rest_sqr * dens_par.two_sig_sqr_inv()) * u.x0();
  const ThreeVector sf_grad = compute_gradient
                                  ? sf * (r + u.threevec() * u_r_scalar) *
                                        dens_par.two_sig_sqr_inv() * 2.0
                                  : ThreeVector(0.0, 0.0, 0.0);

  return std::make_pair(sf, sf_grad);
}

/// \copydoc smash::rho_eckart
template <typename /*ParticlesContainer*/ T>
std::tuple<double, ThreeVector, ThreeVector, ThreeVector> rho_eckart_impl(
    const ThreeVector &r, const T &plist, const DensityParameters &par,
    DensityType dens_type, bool compute_gradient) {
  /* The current density of the positively and negatively charged particles.
   * Division into positive and negative charges is necessary to avoid
   * problems with the Eckart frame definition. Example of problem:
   * get Eckart frame for two identical oppositely flying bunches of
   * electrons and positrons. For this case jmu = (0, 0, 0, non-zero),
   * so jmu.abs does not exist and Eckart frame is not defined.
   * If one takes rho = jmu_pos.abs - jmu_neg.abs, it is still Lorentz-
   * invariant and gives the right limit in non-relativistic case, but
   * it gives no such problem. */
  FourVector jmu_pos, jmu_neg;
  /* The array of the derivatives of the current density.
   * The zeroth component is the time derivative,
   * while the next 3 ones are spacial derivatives. */
  std::array<FourVector, 4> djmu_dx;

  for (const auto &p : plist) {
    const double dens_factor = density_factor(p.type(), dens_type);
    if (std::abs(dens_factor) < really_small) {
      continue;
    }
    const FourVector mom = p.momentum();
    const double m = mom.abs();
    if (m < really_small) {
      continue;
    }
    const double m_inv = 1.0 / m;
    const auto sf_and_grad = unnormalized_smearing_factor(
        p.position().threevec() - r, mom, m_inv, par, compute_gradient);
    if (sf_and_grad.first < really_small) {
      continue;
    }
    const FourVector tmp = mom * (dens_factor / mom.x0());
    if (dens_factor > 0.) {
      jmu_pos += tmp * sf_and_grad.first;
    } else {
      jmu_neg += tmp * sf_and_grad.first;
    }
    if (compute_gradient) {
      for (int k = 1; k <= 3; k++) {
        djmu_dx[k] += tmp * sf_and_grad.second[k - 1];
        djmu_dx[0] -= tmp * sf_and_grad.second[k - 1] * tmp.threevec()[k - 1] /
                      dens_factor;
      }
    }
  }

  // Eckart density
  const double rho_eck = (jmu_pos.abs() - jmu_neg.abs()) * par.norm_factor_sf();

  // $\partial_t \vec j$
  const ThreeVector dj_dt = compute_gradient
                                ? djmu_dx[0].threevec() * par.norm_factor_sf()
                                : ThreeVector(0.0, 0.0, 0.0);

  // Gradient of density
  ThreeVector rho_grad;
  // Curl of current density
  ThreeVector j_rot;
  if (compute_gradient) {
    j_rot.set_x1(djmu_dx[2].x3() - djmu_dx[3].x2());
    j_rot.set_x2(djmu_dx[3].x1() - djmu_dx[1].x3());
    j_rot.set_x3(djmu_dx[1].x2() - djmu_dx[2].x1());
    j_rot *= par.norm_factor_sf();
    for (int i = 1; i < 4; i++) {
      rho_grad[i - 1] += djmu_dx[i].x0() * par.norm_factor_sf();
    }
  }
  return std::make_tuple(rho_eck, rho_grad, dj_dt, j_rot);
}

std::tuple<double, ThreeVector, ThreeVector, ThreeVector> rho_eckart(
    const ThreeVector &r, const ParticleList &plist,
    const DensityParameters &par, DensityType dens_type,
    bool compute_gradient) {
  return rho_eckart_impl(r, plist, par, dens_type, compute_gradient);
}
std::tuple<double, ThreeVector, ThreeVector, ThreeVector> rho_eckart(
    const ThreeVector &r, const Particles &plist, const DensityParameters &par,
    DensityType dens_type, bool compute_gradient) {
  return rho_eckart_impl(r, plist, par, dens_type, compute_gradient);
}

void update_density_lattice(RectangularLattice<DensityOnLattice> *lat,
                            const LatticeUpdate update,
                            const DensityType dens_type,
                            const DensityParameters &par,
                            const Particles &particles,
                            const bool compute_gradient) {
  update_general_lattice(lat, update, dens_type, par, particles,
                         compute_gradient);
}

void update_Tmn_lattice(RectangularLattice<EnergyMomentumTensor> *lat,
                        const LatticeUpdate update, const DensityType dens_type,
                        const DensityParameters &par,
                        const Particles &particles) {
  update_general_lattice(lat, update, dens_type, par, particles);
}

std::ostream &operator<<(std::ostream &os, DensityType dens_type) {
  switch (dens_type) {
    case DensityType::Hadron:
      os << "hadron density";
      break;
    case DensityType::Baryon:
      os << "baryon density";
      break;
    case DensityType::BaryonicIsospin:
      os << "baryonic isospin density";
      break;
    case DensityType::Pion:
      os << "pion density";
      break;
    case DensityType::None:
      os << "none";
      break;
    default:
      os.setstate(std::ios_base::failbit);
  }
  return os;
}

}  // namespace smash
