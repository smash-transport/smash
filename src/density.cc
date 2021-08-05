/*
 *
 *    Copyright (c) 2013-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/density.h"
#include "smash/constants.h"
#include "smash/logging.h"

namespace smash {

double density_factor(const ParticleType &type, DensityType dens_type) {
  switch (dens_type) {
    case DensityType::Hadron:
      return type.is_hadron() ? 1. : 0.;
    case DensityType::Baryon:
      return static_cast<double>(type.baryon_number());
    case DensityType::BaryonicIsospin:
      return type.is_baryon() || type.is_nucleus() ? type.isospin3_rel() : 0.;
    case DensityType::Pion:
      return type.pdgcode().is_pion() ? 1. : 0.;
    case DensityType::Isospin3_tot:
      return type.is_hadron() ? type.isospin3() : 0.;
    case DensityType::Charge:
      return static_cast<double>(type.charge());
    case DensityType::Strangeness:
      return static_cast<double>(type.strangeness());
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

/// \copydoc smash::current_eckart
template <typename /*ParticlesContainer*/ T>
std::tuple<double, FourVector, ThreeVector, ThreeVector, FourVector, FourVector,
           FourVector, FourVector>
current_eckart_impl(const ThreeVector &r, const T &plist,
                    const DensityParameters &par, DensityType dens_type,
                    bool compute_gradient, bool smearing) {
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
  std::array<FourVector, 4> djmu_dxnu;

  for (const auto &p : plist) {
    const double dens_factor = density_factor(p.type(), dens_type);
    if (std::fabs(dens_factor) < really_small) {
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
    const FourVector tmp = mom * (dens_factor / mom.x0());
    if (smearing) {
      if (dens_factor > 0.) {
        jmu_pos += tmp * sf_and_grad.first;
      } else {
        jmu_neg += tmp * sf_and_grad.first;
      }
    } else {
      if (dens_factor > 0.) {
        jmu_pos += tmp;
      } else {
        jmu_neg += tmp;
      }
    }
    if (compute_gradient) {
      for (int k = 1; k <= 3; k++) {
        djmu_dxnu[k] += tmp * sf_and_grad.second[k - 1];
        djmu_dxnu[0] -= tmp * sf_and_grad.second[k - 1] *
                        tmp.threevec()[k - 1] / dens_factor;
      }
    }
  }

  // Eckart density
  const double rho_eck = (jmu_pos.abs() - jmu_neg.abs()) * par.norm_factor_sf();

  // $\partial_t j^{\mu}$
  const FourVector djmu_dt = compute_gradient
                                 ? djmu_dxnu[0] * par.norm_factor_sf()
                                 : FourVector(0.0, 0.0, 0.0, 0.0);
  // $\partial_x j^{\mu}$
  const FourVector djmu_dx = compute_gradient
                                 ? djmu_dxnu[1] * par.norm_factor_sf()
                                 : FourVector(0.0, 0.0, 0.0, 0.0);
  // $\partial_y j^{\mu}$
  const FourVector djmu_dy = compute_gradient
                                 ? djmu_dxnu[2] * par.norm_factor_sf()
                                 : FourVector(0.0, 0.0, 0.0, 0.0);
  // $\partial_z j^{\mu}$
  const FourVector djmu_dz = compute_gradient
                                 ? djmu_dxnu[3] * par.norm_factor_sf()
                                 : FourVector(0.0, 0.0, 0.0, 0.0);

  // Gradient of density
  ThreeVector grad_j0 = ThreeVector(0.0, 0.0, 0.0);
  // Curl of current density
  ThreeVector curl_vecj = ThreeVector(0.0, 0.0, 0.0);
  if (compute_gradient) {
    curl_vecj.set_x1(djmu_dxnu[2].x3() - djmu_dxnu[3].x2());
    curl_vecj.set_x2(djmu_dxnu[3].x1() - djmu_dxnu[1].x3());
    curl_vecj.set_x3(djmu_dxnu[1].x2() - djmu_dxnu[2].x1());
    curl_vecj *= par.norm_factor_sf();
    for (int i = 1; i < 4; i++) {
      grad_j0[i - 1] += djmu_dxnu[i].x0() * par.norm_factor_sf();
    }
  }
  return std::make_tuple(rho_eck, jmu_pos + jmu_neg, grad_j0, curl_vecj,
                         djmu_dt, djmu_dx, djmu_dy, djmu_dz);
}

std::tuple<double, FourVector, ThreeVector, ThreeVector, FourVector, FourVector,
           FourVector, FourVector>
current_eckart(const ThreeVector &r, const ParticleList &plist,
               const DensityParameters &par, DensityType dens_type,
               bool compute_gradient, bool smearing) {
  return current_eckart_impl(r, plist, par, dens_type, compute_gradient,
                             smearing);
}
std::tuple<double, FourVector, ThreeVector, ThreeVector, FourVector, FourVector,
           FourVector, FourVector>
current_eckart(const ThreeVector &r, const Particles &plist,
               const DensityParameters &par, DensityType dens_type,
               bool compute_gradient, bool smearing) {
  return current_eckart_impl(r, plist, par, dens_type, compute_gradient,
                             smearing);
}

void update_lattice(
    RectangularLattice<DensityOnLattice> *lat,
    RectangularLattice<FourVector> *old_jmu,
    RectangularLattice<FourVector> *new_jmu,
    RectangularLattice<std::array<FourVector, 4>> *four_grad_lattice,
    const LatticeUpdate update, const DensityType dens_type,
    const DensityParameters &par, const std::vector<Particles> &ensembles,
    const double time_step, const bool compute_gradient) {
  // Do not proceed if lattice does not exists/update not required
  if (lat == nullptr || lat->when_update() != update) {
    return;
  }
  const std::array<int, 3> lattice_n_cells = lat->n_cells();
  const int number_of_nodes =
      lattice_n_cells[0] * lattice_n_cells[1] * lattice_n_cells[2];

  /*
   * Take the provided DensityOnLattice lattice and use the information about
   * the current to create a lattice of current FourVectors. Because the lattice
   * hasn't been updated at this point yet, it provides the t_0 time step
   * information on the currents.
   */
  // copy values of jmu at t_0 onto old_jmu;
  // proceed only if finite difference gradients are calculated
  if (par.derivatives() == DerivativesMode::FiniteDifference) {
    for (int i = 0; i < number_of_nodes; i++) {
      old_jmu->assign_value(i, ((*lat)[i]).jmu_net());
    }
  }

  update_lattice(lat, update, dens_type, par, ensembles, compute_gradient);

  // calculate the gradients for finite difference derivatives
  if (par.derivatives() == DerivativesMode::FiniteDifference) {
    // copy values of jmu FourVectors at t_0 + time_step onto new_jmu
    for (int i = 0; i < number_of_nodes; i++) {
      new_jmu->assign_value(i, ((*lat)[i]).jmu_net());
    }

    // compute time derivatives and gradients of all components of jmu
    new_jmu->compute_four_gradient_lattice(*old_jmu, time_step,
                                           *four_grad_lattice);

    // substitute new derivatives
    int node_number = 0;
    for (auto &node : *lat) {
      auto tmp = (*four_grad_lattice)[node_number];
      node.overwrite_djmu_dxnu(tmp[0], tmp[1], tmp[2], tmp[3]);
      node_number++;
    }
  }  // if (par.derivatives() == DerivativesMode::FiniteDifference)

  // calculate gradients of rest frame density
  if (par.rho_derivatives() == RestFrameDensityDerivativesMode::On) {
    for (auto &node : *lat) {
      // the rest frame density
      double rho = node.rho();
      const int sgn = rho > 0 ? 1 : -1;
      if (std::abs(rho) < very_small_double) {
        rho = sgn * very_small_double;
      }

      // the computational frame j^mu
      const FourVector jmu = node.jmu_net();
      // computational frame array of derivatives of j^mu
      const std::array<FourVector, 4> djmu_dxnu = node.djmu_dxnu();

      const double drho_dt =
          (1 / rho) *
          (jmu.x0() * djmu_dxnu[0].x0() - jmu.x1() * djmu_dxnu[0].x1() -
           jmu.x2() * djmu_dxnu[0].x2() - jmu.x3() * djmu_dxnu[0].x3());

      const double drho_dx =
          (1 / rho) *
          (jmu.x0() * djmu_dxnu[1].x0() - jmu.x1() * djmu_dxnu[1].x1() -
           jmu.x2() * djmu_dxnu[1].x2() - jmu.x3() * djmu_dxnu[1].x3());

      const double drho_dy =
          (1 / rho) *
          (jmu.x0() * djmu_dxnu[2].x0() - jmu.x1() * djmu_dxnu[2].x1() -
           jmu.x2() * djmu_dxnu[2].x2() - jmu.x3() * djmu_dxnu[2].x3());

      const double drho_dz =
          (1 / rho) *
          (jmu.x0() * djmu_dxnu[3].x0() - jmu.x1() * djmu_dxnu[3].x1() -
           jmu.x2() * djmu_dxnu[3].x2() - jmu.x3() * djmu_dxnu[3].x3());

      const FourVector drho_dxnu = {drho_dt, drho_dx, drho_dy, drho_dz};

      node.overwrite_drho_dxnu(drho_dxnu);
    }

  }  // if (par.rho_derivatives() == RestFrameDensityDerivatives::On){

}  // void update_lattice()

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
    case DensityType::Isospin3_tot:
      os << "total isospin3 density";
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
