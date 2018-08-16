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
                                  ? sf * (r + u.threevec() * u_r_scalar)
                                    * dens_par.two_sig_sqr_inv() * 2.0
                                  : ThreeVector(0.0, 0.0, 0.0);

  return std::make_pair(sf, sf_grad);
}

/**
 * Calculates Eckart rest frame density and optionally the gradient of the
 * density in an arbitary frame, the curl of the 3-currrent and the time
 * derivative of the 3-current.
 * \f[j^{\mu} = (\sqrt{2\pi} \sigma )^{-3} \sum_{i=1}^N C_i u^{\mu}_i
 * exp \left(- \frac{(\vec r -\vec r_i + \frac{\gamma_i^2}{1 + \gamma_i}
 * \vec \beta_i (\vec \beta_i, \vec r - \vec r_i))^2}{2\sigma^2} \right)\f]
 * \f[ \rho^{Eckart} = \sqrt{j^{\mu} j_{\mu}} \f]
 * Here \f$ C_i \f$ is a corresponding value of "charge". If baryon
 * current option is selected then \f$ C_i \f$ is 1 for baryons,
 * -1 for antibaryons and 0 otherwise. For proton/neutron
 * current \f$ C_i = 1\f$ for proton/neutron and 0 otherwise.
 *
 * To avoid the problems with Eckart frame definition, densities for
 * positive and negative charges, \f$\rho_+ \f$ and \f$ \rho_-\f$,
 * are computed separately and result is \f$\rho_+ - \rho_-\f$.
 *
 * \param[in] r Arbitrary space point where 4-current is calculated [fm]
 * \param[in] plist List of all particles to be used in \f$j^{\mu}\f$
 *            calculation. If the distance between particle and calculation
 *            point r, \f$ |r-r_i| > r_{cut} \f$ then particle input
 *            to density will be ignored.
 *
 * Next three values are taken from ExperimentalParameters structure:
 *
 * \param[in] par Set of parameters packed in one structure.
 *            From them the cutting radius r_cut \f$ r_{cut} / \sigma \f$,
 *            number of test-particles ntest and the gaussian width
 *            gs_sigma are needed.
 * \param[in] dens_type type of four-currect to be calculated:
 *            baryon, proton or neutron options are currently available
 * \param[in] compute_gradient true - compute gradient, false - no
 * \fpPrecision Density gradient is returned as double, because it is
 *   ThreeVector and ThreeVector currently comes only as double.
 *   Density itself is double for uniformity: if gradient is double,
 *   density should also be.
 * \tparam T ParticlesContainer
 * \return (density in the local Eckart frame [fm\$f^{-3}\$f],
 *          \f$ \nabla\cdots\rho \f$ or a 0 3-vector,
 *          \f$ \partial_t \vec j\f$ or a 0 3-vector,
 *          \f$ \nabla \times \vec j \f$ or a 0 3-vector).
 */
template <typename /*ParticlesContainer*/ T>
std::tuple<double, ThreeVector, ThreeVector, ThreeVector> rho_eckart_impl(
                                               const ThreeVector &r,
                                               const T &plist,
                                               const DensityParameters &par,
                                               DensityType dens_type,
                                               bool compute_gradient) {
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
        djmu_dx[0] -= tmp * sf_and_grad.second[k - 1]
                      * tmp.threevec()[k-1] / dens_factor;
      }
    }
  }

  // Eckart density
  const double rho_eck =
      (jmu_pos.abs() - jmu_neg.abs()) * par.norm_factor_sf();

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
                                          const ThreeVector &r,
                                          const ParticleList &plist,
                                          const DensityParameters &par,
                                          DensityType dens_type,
                                          bool compute_gradient) {
  return rho_eckart_impl(r, plist, par, dens_type, compute_gradient);
}
std::tuple<double, ThreeVector, ThreeVector, ThreeVector> rho_eckart(
                                          const ThreeVector &r,
                                          const Particles &plist,
                                          const DensityParameters &par,
                                          DensityType dens_type,
                                          bool compute_gradient) {
  return rho_eckart_impl(r, plist, par, dens_type, compute_gradient);
}

void update_density_lattice(RectangularLattice<DensityOnLattice> *lat,
                            const LatticeUpdate update,
                            const DensityType dens_type,
                            const DensityParameters &par,
                            const Particles &particles,
                            bool compute_gradient) {
  // Do not proceed if lattice does not exists/update not required
  if (lat == nullptr || lat->when_update() != update) {
    return;
  }
  lat->reset();
  const double norm_factor = par.norm_factor_sf();
  for (const auto &part : particles) {
    const double dens_factor = density_factor(part.type(), dens_type);
    if (std::abs(dens_factor) < really_small) {
      continue;
    }
    const FourVector p = part.momentum();
    const double m = p.abs();
    if (unlikely(m < really_small)) {
      const auto &log = logger<LogArea::Density>();
      log.warn("Gaussian smearing is undefined for momentum ", p);
      continue;
    }
    const double m_inv = 1.0 / m;

    const ThreeVector pos = part.position().threevec();
    lat->iterate_in_radius(
        pos, par.r_cut(), [&](DensityOnLattice &node, int ix, int iy, int iz) {
          const ThreeVector r = lat->cell_center(ix, iy, iz);
          const double sf =
              norm_factor *
              unnormalized_smearing_factor(pos - r, p, m_inv,
                                           par, compute_gradient).first;
          const ThreeVector sf_grad =
              norm_factor *
              unnormalized_smearing_factor(pos - r, p, m_inv,
                                           par, compute_gradient).second;
          if (sf > really_small) {
            node.add_particle(part, dens_factor, sf, compute_gradient, sf_grad);
          }
        });
  }
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
