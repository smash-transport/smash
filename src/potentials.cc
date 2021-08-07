/*
 *
 *    Copyright (c) 2014-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/potentials.h"

#include "smash/constants.h"
#include "smash/density.h"

namespace smash {

Potentials::Potentials(Configuration conf, const DensityParameters &param)
    : param_(param),
      use_skyrme_(conf.has_value({"Skyrme"})),
      use_symmetry_(conf.has_value({"Symmetry"})),
      use_vdf_(conf.has_value({"VDF"})) {
  /*!\Userguide
   * \page input_potentials_ Potentials
   * SMASH simulation supports two sets of potentials:\n
   * 1) Skyrme and/or Symmetry potentials\n
   * 2) VDF (vector density functional) model potentials,
   * https://arxiv.org/pdf/2011.06635.pdf \n
   *
   * Skyrme and VDF potentials both describe the behavior of symmetric nuclear
   * matter. The symmetry potential can adjust the Skyrme potential (but not the
   * VDF potential) to include effects due to isospin. The Skyrme and Symmetry
   * potentials are semi-relativistic, while the VDF potential is fully
   * relativistic.
   *
   * \li \subpage potentials_skyrme_
   * \li \subpage potentials_sym_
   * \li \subpage potentials_VDF_
   *
   * \page potentials_skyrme_ Skyrme
   *
   * The Skyrme potential has the form
   * \f[ U_{Sk} = A(\rho/\rho_0) + B (\rho/\rho_0)^{\tau} \,, \f]
   * where \f$ \rho \f$ is baryon density in the local Eckart rest frame.
   *
   * \key Skyrme_A (double, required, no default): \n
   *      Parameter A of Skyrme potential in MeV
   *
   * \key Skyrme_B (double, required, no default): \n
   *      Parameter B of Skyrme potential in MeV
   *
   * \key Skyrme_Tau (double, required, no default): \n
   *      Parameter \f$\tau\f$ of Skyrme potential.
   *
   * \page input_potentials_ Potentials
   * \n
   * **Example: Configuring Skyrme Potentials**\n
   *
   * The following extract from the configuration file configures SMASH such
   * that the Skyrme as well as the Symmetry potential are activated for the
   * simulation. There is however no requirement to include both simultaneously.
   They
   * can be switched on and off individually.
   * \n
   *\verbatim
   Potentials:
       Skyrme:
           Skyrme_A: -209.2
           Skyrme_B: 156.4
           Skyrme_Tau: 1.35
       Symmetry:
           S_Pot: 18.0
   \endverbatim
   */
  if (use_skyrme_) {
    skyrme_a_ = conf.take({"Skyrme", "Skyrme_A"});
    skyrme_b_ = conf.take({"Skyrme", "Skyrme_B"});
    skyrme_tau_ = conf.take({"Skyrme", "Skyrme_Tau"});
  }

  /*!\Userguide
   * \page potentials_sym_ Symmetry
   *
   * The symmetry potential has the form
   * \f[ U_{Sym} = \pm 2 S_{pot} \frac{I_3}{I} \frac{\rho_{I_3}}{\rho_0}
   * + S(\rho_B)\left(\frac{\rho_{I_3}}{\rho_B}\right)^2 \,, \f]
   * where \f$ \rho_{I_3}\f$ is the density of the relative isospin \f$ I_3/I
   * \f$ and \f$ \rho_B \f$ is the net baryon density.
   *
   * \key S_pot (double, required, no default): \n
   *      Parameter \f$S_{pot}\f$ of symmetry potential in MeV
   *
   * \key gamma (double, no default): \n
   *      Power \f$ \gamma \f$ in formula for \f$ S(\rho_B) \f$:
   * \f[ S(\rho_B)=12.3\,\mathrm{MeV}\times
   * \left(\frac{\rho_B}{\rho_0}\right)^{2/3}+
   * 20\,\mathrm{MeV}\times\left(\frac{\rho_B}{\rho_0}\right)^\gamma \f]
   * If gamma is specified, the baryon density dependence is included in the
   * potential. Otherwise only the first term of the potential will be taken
   * into account.
   */
  if (use_symmetry_) {
    symmetry_S_Pot_ = conf.take({"Symmetry", "S_Pot"});
    if (conf.has_value({"Symmetry", "gamma"})) {
      symmetry_gamma_ = conf.take({"Symmetry", "gamma"});
      symmetry_is_rhoB_dependent_ = true;
    }
  }
  /*!\Userguide
   * \page potentials_VDF_ VDF
   *
   * The VDF potential is a four-vector of the form
   * \f[A^{\mu} = \sum_{i=1}^N C_i \left(\frac{\rho}{\rho_0}\right)^{b_i - 2}
   * \frac{j^{\mu}}{\rho_0} \,,\f] where \f$j^{\mu}\f$ is baryon 4-current,
   * \f$\rho\f$ is baryon density in the local Eckart rest frame, and
   * \f$\rho_0\f$ is the saturation density. The parameters of the potential,
   * the coefficients \f$C_i\f$ and the powers \f$b_i\f$, are fitted to
   * reproduce a chosen set of properties of dense nuclear matter, and in
   * particular these may include describing two first order phase transitions:
   * the well-known phase transition in ordinary nuclear matter, and a transiton
   * at high baryon densities meant to model a possible QCD phase transition (a
   * "QGP-like" phase transition); see https://arxiv.org/pdf/2011.06635.pdf for
   * details and example parameter sets for the case \f$N=4\f$.
   * \n
   * The user can decide how many terms \f$N\f$ should enter the potential by
   * populating the coefficients and powers vectors in the config file with a
   * chosen number of entries. The number of coefficients must match the number
   * of powers.
   *
   * \key Sat_rhoB (double, required, default 0.160): \n
   *      The saturation density of nuclear matter, in fm\f$^{-3}\f$
   *
   * \key Coeffs (array<double, N>, required, no default): \n
   *      Parameters \f$C_i\f$ of the VDF potential in MeV (note: the code
   *      automatically converts these to GeV)
   *
   * \key Powers (array<double, N>, required, no default): \n
   *      Parameters \f$b_i\f$ of the VDF potential
   *
   * \page input_potentials_ Potentials
   * \n
   * **Example: Configuring VDF Potentials**\n
   *
   * The following extracts from the configuration file configure SMASH such
   * that the VDF potential is activated for the simulation. In the first
   * example, VDF potentials are configured to reproduce the default SMASH
   * Skyrme potentials (without the symmetry potential, as it is not described
   * within the VDF model):
   * \n
   *\verbatim
   Potentials:
       VDF:
           Sat_rhoB: 0.168
           Powers: [2.0, 2.35]
           Coeffs: [-209.2, 156.5]
   \endverbatim
   *
   * In the second example, VDF potentials are configured to describe nuclear
   * matter with saturation density of \f$\rho_0 =\f$ 0.160 fm\f$^{-3}\f$,
   * binding energy of \f$B_0 = -16.3\f$ MeV, the critical point of the ordinary
   * nuclear liquid-gas phase transition at \f$T_c^{(N)} = 18\f$ MeV and
   * \f$\rho_c^{(N)} = 0.375 \rho_0\f$, the critical point of the conjectured
   * "QGP-like" phase transition at \f$T_c^{(Q)} = 100\f$ MeV and
   * \f$\rho_c^{(Q)} = 3.0\rho_0\f$, and the boundaries of the spinodal region
   * of the "QGP-like" phase transition at \f$\eta_L = 2.50 \rho_0\f$ and
   * \f$\eta_R = 3.315 \rho_0\f$:
   * \n
   *\verbatim
   Potentials:
       VDF:
           Sat_rhoB: 0.160
           Powers: [1.7681391, 3.5293515, 5.4352788, 6.3809822]
           Coeffs: [-8.450948e+01, 3.843139e+01, -7.958557e+00, 1.552594e+00]
   \endverbatim
   */
  if (use_vdf_) {
    saturation_density_ = conf.take({"VDF", "Sat_rhoB"});
    std::vector<double> aux_coeffs = conf.take({"VDF", "Coeffs"});
    std::vector<double> aux_powers = conf.take({"VDF", "Powers"});
    if (aux_coeffs.size() != aux_powers.size()){
      throw std::invalid_argument(
        "The number of coefficients should equal the number of powers.");
    }
    number_of_terms_ = aux_powers.size();
    for (int i = 0; i < number_of_terms_; i++){
      if (aux_powers[i] < 0.0){
	throw std::invalid_argument(
          "Powers need to be positive real numbers.");
      }
      // coefficients are provided in MeV, but the code uses GeV
      coeffs_.push_back(aux_coeffs[i] * mev_to_gev);
      powers_.push_back(aux_powers[i]);
    }
  }
}

Potentials::~Potentials() {}

double Potentials::skyrme_pot(const double baryon_density) const {
  const double tmp = baryon_density / nuclear_density;
  /* U = U(|rho|) * sgn , because the sign of the potential changes
   * under a charge reversal transformation. */
  const int sgn = tmp > 0 ? 1 : -1;
  // Return in GeV
  return 1.0e-3 * sgn *
         (skyrme_a_ * std::abs(tmp) +
          skyrme_b_ * std::pow(std::abs(tmp), skyrme_tau_));
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
  double pot =
      1.0e-3 * 2. * symmetry_S_Pot_ * baryon_isospin_density / nuclear_density;
  if (symmetry_is_rhoB_dependent_) {
    pot += 1.0e-3 * symmetry_S(baryon_density) * baryon_isospin_density *
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
  for (int i = 0; i < number_of_terms_; i++){
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
    for (int i = 0; i < number_of_terms_; i++){
      F_1 += coeffs_[i] * (powers_[i] - 2.0) *
	std::pow(abs_rhoB, powers_[i] - 3.0) /
	std::pow(saturation_density_, powers_[i] - 1.0);
    }
    F_1 = F_1 * sgn;

    double F_2 = 0.0;
    for (int i = 0; i < number_of_terms_; i++){
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
