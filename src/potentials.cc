/*
 *
 *    Copyright (c) 2014-2015,2017-2023
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/potentials.h"

#include "smash/constants.h"
#include "smash/density.h"
#include "smash/input_keys.h"

namespace smash {

Potentials::Potentials(Configuration conf, const DensityParameters &param)
    : param_(param),
      use_skyrme_(conf.has_value({"Skyrme"})),
      use_symmetry_(conf.has_value({"Symmetry"})),
      use_coulomb_(conf.has_value({"Coulomb"})),
      use_vdf_(conf.has_value({"VDF"})),
      use_momentum_dependence_(conf.has_value({"Momentum_Dependence"})),
      use_potentials_outside_lattice_(
          conf.take({"Use_Potentials_Outside_Lattice"},
                    InputKeys::potentials_use_potentials_outside_lattice
                        .default_value())) {
  if (use_skyrme_) {
    skyrme_a_ = conf.take({"Skyrme", "Skyrme_A"});
    skyrme_b_ = conf.take({"Skyrme", "Skyrme_B"});
    skyrme_tau_ = conf.take({"Skyrme", "Skyrme_Tau"});
  }
  if (use_momentum_dependence_) {
    mom_dependence_Lambda_ = conf.take({"Momentum_Dependence", "Lambda"});
    mom_dependence_C_ = conf.take({"Momentum_Dependence", "C"});
  }

  if (use_symmetry_) {
    symmetry_S_Pot_ = conf.take({"Symmetry", "S_Pot"});
    if (conf.has_value({"Symmetry", "gamma"})) {
      symmetry_gamma_ = conf.take({"Symmetry", "gamma"});
      symmetry_is_rhoB_dependent_ = true;
    }
  }
  if (use_coulomb_) {
    coulomb_r_cut_ = conf.take({"Coulomb", "R_Cut"});
  }
  if (use_vdf_) {
    saturation_density_ = conf.take({"VDF", "Sat_rhoB"});
    std::vector<double> aux_coeffs = conf.take({"VDF", "Coeffs"});
    std::vector<double> aux_powers = conf.take({"VDF", "Powers"});
    if (aux_coeffs.size() != aux_powers.size()) {
      throw std::invalid_argument(
          "The number of coefficients should equal the number of powers.");
    }
    const int n_terms = aux_powers.size();
    for (int i = 0; i < n_terms; i++) {
      if (aux_powers[i] < 0.0) {
        throw std::invalid_argument("Powers need to be positive real numbers.");
      }
      // coefficients are provided in MeV, but the code uses GeV
      coeffs_.push_back(aux_coeffs[i] * mev_to_gev);
      powers_.push_back(aux_powers[i]);
    }
  }
}

Potentials::~Potentials() {}

double Potentials::skyrme_pot(const double baryon_density, const double A,
                              const double B, const double tau) {
  const double tmp = baryon_density / nuclear_density;
  /* U = U(|rho|) * sgn , because the sign of the potential changes
   * under a charge reversal transformation. */
  const int sgn = tmp > 0 ? 1 : -1;
  // Return in GeV
  return mev_to_gev * sgn *
         (A * std::abs(tmp) + B * std::pow(std::abs(tmp), tau));
}

double Potentials::symmetry_S(const double baryon_density) const {
  if (symmetry_is_rhoB_dependent_) {
    return 12.3 * std::pow(baryon_density / nuclear_density, 2. / 3.) +
           20.0 * std::pow(baryon_density / nuclear_density, symmetry_gamma_);
  } else {
    return 0.;
  }
}
double Potentials::symmetry_pot(const double baryon_isospin_density,
                                const double baryon_density) const {
  double pot = mev_to_gev * 2. * symmetry_S_Pot_ * baryon_isospin_density /
               nuclear_density;
  if (symmetry_is_rhoB_dependent_) {
    pot += mev_to_gev * symmetry_S(baryon_density) * baryon_isospin_density *
           baryon_isospin_density / (baryon_density * baryon_density);
  }
  return pot;
}

FourVector Potentials::vdf_pot(double rhoB, const FourVector jmuB_net) const {
  // this needs to be used in order to prevent trying to calculate something
  // like
  // (-rho_B)^{3.4}
  const int sgn = rhoB > 0 ? 1 : -1;
  double abs_rhoB = std::abs(rhoB);
  // to prevent NAN expressions
  if (abs_rhoB < very_small_double) {
    abs_rhoB = very_small_double;
  }
  // F_2 is a multiplicative factor in front of the baryon current
  // in the VDF potential
  double F_2 = 0.0;
  for (int i = 0; i < number_of_terms(); i++) {
    F_2 += coeffs_[i] * std::pow(abs_rhoB, powers_[i] - 2.0) /
           std::pow(saturation_density_, powers_[i] - 1.0);
  }
  F_2 = F_2 * sgn;
  // Return in GeV
  return F_2 * jmuB_net;
}

double Potentials::potential(const ThreeVector &r, const ParticleList &plist,
                             const ParticleType &acts_on) const {
  double total_potential = 0.0;
  const bool compute_gradient = false;
  const bool smearing = true;
  const auto scale = force_scale(acts_on);

  if (!(acts_on.is_baryon() || acts_on.is_nucleus())) {
    return total_potential;
  }
  const auto baryon_density_and_gradient = current_eckart(
      r, plist, param_, DensityType::Baryon, compute_gradient, smearing);
  const double rhoB = std::get<0>(baryon_density_and_gradient);
  if (use_skyrme_) {
    total_potential += scale.first * skyrme_pot(rhoB);
  }
  if (use_symmetry_) {
    const double rho_iso = std::get<0>(
        current_eckart(r, plist, param_, DensityType::BaryonicIsospin,
                       compute_gradient, smearing));
    const double sym_pot = symmetry_pot(rho_iso, rhoB) * acts_on.isospin3_rel();
    total_potential += scale.second * sym_pot;
  }

  if (use_vdf_) {
    const FourVector jmuB = std::get<1>(baryon_density_and_gradient);
    const FourVector VDF_potential = vdf_pot(rhoB, jmuB);
    total_potential += scale.first * VDF_potential.x0();
  }

  return total_potential;
}

std::pair<double, int> Potentials::force_scale(const ParticleType &data) {
  const auto &pdg = data.pdgcode();
  const double skyrme_or_VDF_scale =
      (3 - std::abs(pdg.strangeness())) / 3. * pdg.baryon_number();
  const int symmetry_scale = pdg.baryon_number();
  return std::make_pair(skyrme_or_VDF_scale, symmetry_scale);
}

std::pair<ThreeVector, ThreeVector> Potentials::skyrme_force(
    const double rhoB, const ThreeVector grad_j0B, const ThreeVector dvecjB_dt,
    const ThreeVector curl_vecjB) const {
  ThreeVector E_component(0.0, 0.0, 0.0), B_component(0.0, 0.0, 0.0);
  if (use_skyrme_) {
    const int sgn = rhoB > 0 ? 1 : -1;
    const double abs_rhoB = std::abs(rhoB);
    const double dV_drho = sgn *
                           (skyrme_a_ + skyrme_b_ * skyrme_tau_ *
                                            std::pow(abs_rhoB / nuclear_density,
                                                     skyrme_tau_ - 1)) *
                           mev_to_gev / nuclear_density;
    E_component -= dV_drho * (grad_j0B + dvecjB_dt);
    B_component += dV_drho * curl_vecjB;
  }
  return std::make_pair(E_component, B_component);
}

std::pair<ThreeVector, ThreeVector> Potentials::symmetry_force(
    const double rhoI3, const ThreeVector grad_j0I3,
    const ThreeVector dvecjI3_dt, const ThreeVector curl_vecjI3,
    const double rhoB, const ThreeVector grad_j0B, const ThreeVector dvecjB_dt,
    const ThreeVector curl_vecjB) const {
  ThreeVector E_component(0.0, 0.0, 0.0), B_component(0.0, 0.0, 0.0);
  if (use_symmetry_) {
    E_component -= dVsym_drhoI3(rhoB, rhoI3) * (grad_j0I3 + dvecjI3_dt) +
                   dVsym_drhoB(rhoB, rhoI3) * (grad_j0B + dvecjB_dt);
    B_component += dVsym_drhoI3(rhoB, rhoI3) * curl_vecjI3 +
                   dVsym_drhoB(rhoB, rhoI3) * curl_vecjB;
  }
  return std::make_pair(E_component, B_component);
}

std::pair<ThreeVector, ThreeVector> Potentials::vdf_force(
    double rhoB, const double drhoB_dt, const ThreeVector grad_rhoB,
    const ThreeVector gradrhoB_cross_vecjB, const double j0B,
    const ThreeVector grad_j0B, const ThreeVector vecjB,
    const ThreeVector dvecjB_dt, const ThreeVector curl_vecjB) const {
  // this needs to be used to prevent trying to calculate something like
  // (-rhoB)^{3.4}
  const int sgn = rhoB > 0 ? 1 : -1;
  ThreeVector E_component(0.0, 0.0, 0.0), B_component(0.0, 0.0, 0.0);
  // to prevent NAN expressions
  double abs_rhoB = std::abs(rhoB);
  if (abs_rhoB < very_small_double) {
    abs_rhoB = very_small_double;
  }
  if (use_vdf_) {
    // F_1 and F_2 are multiplicative factors in front of the baryon current
    // in the VDF potential
    double F_1 = 0.0;
    for (int i = 0; i < number_of_terms(); i++) {
      F_1 += coeffs_[i] * (powers_[i] - 2.0) *
             std::pow(abs_rhoB, powers_[i] - 3.0) /
             std::pow(saturation_density_, powers_[i] - 1.0);
    }
    F_1 = F_1 * sgn;

    double F_2 = 0.0;
    for (int i = 0; i < number_of_terms(); i++) {
      F_2 += coeffs_[i] * std::pow(abs_rhoB, powers_[i] - 2.0) /
             std::pow(saturation_density_, powers_[i] - 1.0);
    }
    F_2 = F_2 * sgn;

    E_component -= (F_1 * (grad_rhoB * j0B + drhoB_dt * vecjB) +
                    F_2 * (grad_j0B + dvecjB_dt));
    B_component += F_1 * gradrhoB_cross_vecjB + F_2 * curl_vecjB;
  }
  return std::make_pair(E_component, B_component);
}

// overload of the above
std::pair<ThreeVector, ThreeVector> Potentials::vdf_force(
    const ThreeVector grad_A_0, const ThreeVector dA_dt,
    const ThreeVector curl_A) const {
  ThreeVector E_component(0.0, 0.0, 0.0), B_component(0.0, 0.0, 0.0);
  if (use_vdf_) {
    E_component -= (grad_A_0 + dA_dt);
    B_component += curl_A;
  }
  return std::make_pair(E_component, B_component);
}

double Potentials::dVsym_drhoI3(const double rhoB, const double rhoI3) const {
  double term1 = 2. * symmetry_S_Pot_ / nuclear_density;
  if (symmetry_is_rhoB_dependent_) {
    double term2 = 2. * rhoI3 * symmetry_S(rhoB) / (rhoB * rhoB);
    return mev_to_gev * (term1 + term2);
  } else {
    return mev_to_gev * term1;
  }
}

double Potentials::dVsym_drhoB(const double rhoB, const double rhoI3) const {
  if (symmetry_is_rhoB_dependent_) {
    double rhoB_over_rho0 = rhoB / nuclear_density;
    double term1 = 8.2 * std::pow(rhoB_over_rho0, -1. / 3.) / nuclear_density +
                   20. * symmetry_gamma_ *
                       std::pow(rhoB_over_rho0, symmetry_gamma_) / rhoB;
    double term2 = -2. * symmetry_S(rhoB) / rhoB;
    return mev_to_gev * (term1 + term2) * rhoI3 * rhoI3 / (rhoB * rhoB);
  } else {
    return 0.;
  }
}

std::tuple<ThreeVector, ThreeVector, ThreeVector, ThreeVector>
Potentials::all_forces(const ThreeVector &r, const ParticleList &plist) const {
  const bool compute_gradient = true;
  const bool smearing = true;
  auto F_skyrme_or_VDF =
      std::make_pair(ThreeVector(0., 0., 0.), ThreeVector(0., 0., 0.));
  auto F_symmetry =
      std::make_pair(ThreeVector(0., 0., 0.), ThreeVector(0., 0., 0.));

  const auto baryon_density_and_gradient = current_eckart(
      r, plist, param_, DensityType::Baryon, compute_gradient, smearing);
  double rhoB = std::get<0>(baryon_density_and_gradient);
  const ThreeVector grad_j0B = std::get<2>(baryon_density_and_gradient);
  const ThreeVector curl_vecjB = std::get<3>(baryon_density_and_gradient);
  const FourVector djmuB_dt = std::get<4>(baryon_density_and_gradient);
  if (use_skyrme_) {
    F_skyrme_or_VDF =
        skyrme_force(rhoB, grad_j0B, djmuB_dt.threevec(), curl_vecjB);
  }

  if (use_symmetry_) {
    const auto density_and_gradient =
        current_eckart(r, plist, param_, DensityType::BaryonicIsospin,
                       compute_gradient, smearing);
    const double rhoI3 = std::get<0>(density_and_gradient);
    const ThreeVector grad_j0I3 = std::get<2>(density_and_gradient);
    const ThreeVector curl_vecjI3 = std::get<3>(density_and_gradient);
    const FourVector dvecjI3_dt = std::get<4>(density_and_gradient);
    F_symmetry =
        symmetry_force(rhoI3, grad_j0I3, dvecjI3_dt.threevec(), curl_vecjI3,
                       rhoB, grad_j0B, djmuB_dt.threevec(), curl_vecjB);
  }

  if (use_vdf_) {
    const FourVector jmuB = std::get<1>(baryon_density_and_gradient);
    const FourVector djmuB_dx = std::get<5>(baryon_density_and_gradient);
    const FourVector djmuB_dy = std::get<6>(baryon_density_and_gradient);
    const FourVector djmuB_dz = std::get<7>(baryon_density_and_gradient);

    // safety check to not divide by zero
    const int sgn = rhoB > 0 ? 1 : -1;
    if (std::abs(rhoB) < very_small_double) {
      rhoB = sgn * very_small_double;
    }

    const double drhoB_dt =
        (1 / rhoB) * (jmuB.x0() * djmuB_dt.x0() - jmuB.x1() * djmuB_dt.x1() -
                      jmuB.x2() * djmuB_dt.x2() - jmuB.x3() * djmuB_dt.x3());

    const double drhoB_dx =
        (1 / rhoB) * (jmuB.x0() * djmuB_dx.x0() - jmuB.x1() * djmuB_dx.x1() -
                      jmuB.x2() * djmuB_dx.x2() - jmuB.x3() * djmuB_dx.x3());

    const double drhoB_dy =
        (1 / rhoB) * (jmuB.x0() * djmuB_dy.x0() - jmuB.x1() * djmuB_dy.x1() -
                      jmuB.x2() * djmuB_dy.x2() - jmuB.x3() * djmuB_dy.x3());

    const double drhoB_dz =
        (1 / rhoB) * (jmuB.x0() * djmuB_dz.x0() - jmuB.x1() * djmuB_dz.x1() -
                      jmuB.x2() * djmuB_dz.x2() - jmuB.x3() * djmuB_dz.x3());

    const FourVector drhoB_dxnu = {drhoB_dt, drhoB_dx, drhoB_dy, drhoB_dz};

    const ThreeVector grad_rhoB = drhoB_dxnu.threevec();
    const ThreeVector vecjB = jmuB.threevec();
    const ThreeVector Drho_cross_vecj = grad_rhoB.cross_product(vecjB);

    F_skyrme_or_VDF = vdf_force(
        rhoB, drhoB_dt, drhoB_dxnu.threevec(), Drho_cross_vecj, jmuB.x0(),
        grad_j0B, jmuB.threevec(), djmuB_dt.threevec(), curl_vecjB);
  }

  return std::make_tuple(F_skyrme_or_VDF.first, F_skyrme_or_VDF.second,
                         F_symmetry.first, F_symmetry.second);
}

}  // namespace smash
