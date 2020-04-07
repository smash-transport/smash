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
      use_symmetry_(conf.has_value({"Symmetry"})) {
  /*!\Userguide
   * \page input_potentials_ Potentials
   * Skyrme and/or Symmetry potentials can be accounted for within a SMASH
   * simulation.
   *
   * Currently, potentials are just added without re-adjusting the energy
   * and momenta of the colliding nucleons. This can be done, because the
   * binding energy of nucleons is between 0 and 8 MeV per nucleon while the
   * kinetic energies, at which SMASH operates, are at least 400 MeV per
   * nucleon. The binding energy can thus be neglected compared to the
   * kinetic energy.
   *
   * \li \subpage potentials_skyrme_
   * \li \subpage potentials_sym_
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
   * **Example: Configuring Potentials**\n
   *
   * The following extract from the configuration file configures SMASH such
   * that the Skyrme as well as the Symmetry potential are activated for the
   * simulation. There is however no necessity to include both simultaneously.
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

double Potentials::potential(const ThreeVector &r, const ParticleList &plist,
                             const ParticleType &acts_on) const {
  double total_potential = 0.0;
  const bool compute_gradient = false;
  const bool smearing = true;
  const auto scale = force_scale(acts_on);

  if (!(acts_on.is_baryon() || acts_on.is_nucleus())) {
    return total_potential;
  }
  const double rho_B = std::get<0>(current_eckart(
      r, plist, param_, DensityType::Baryon, compute_gradient, smearing));
  if (use_skyrme_) {
    total_potential += scale.first * skyrme_pot(rho_B);
  }
  if (use_symmetry_) {
    const double rho_iso = std::get<0>(
        current_eckart(r, plist, param_, DensityType::BaryonicIsospin,
                       compute_gradient, smearing));
    const double sym_pot =
        symmetry_pot(rho_iso, rho_B) * acts_on.isospin3_rel();
    total_potential += scale.second * sym_pot;
  }
  return total_potential;
}

std::pair<double, int> Potentials::force_scale(const ParticleType &data) {
  const auto &pdg = data.pdgcode();
  const double skyrme_scale =
      (3 - std::abs(pdg.strangeness())) / 3. * pdg.baryon_number();
  const int symmetry_scale = pdg.baryon_number();
  return std::make_pair(skyrme_scale, symmetry_scale);
}

std::pair<ThreeVector, ThreeVector> Potentials::skyrme_force(
    const double density, const ThreeVector grad_rho, const ThreeVector dj_dt,
    const ThreeVector rot_j) const {
  ThreeVector E_component(0.0, 0.0, 0.0), B_component(0.0, 0.0, 0.0);
  if (use_skyrme_) {
    const double dV_drho =
        (skyrme_a_ + skyrme_b_ * skyrme_tau_ *
                         std::pow(density / nuclear_density, skyrme_tau_ - 1)) *
        mev_to_gev / nuclear_density;
    E_component -= dV_drho * (grad_rho + dj_dt);
    B_component += dV_drho * rot_j;
  }
  return std::make_pair(E_component, B_component);
}

std::pair<ThreeVector, ThreeVector> Potentials::symmetry_force(
    const double rhoI3, const ThreeVector grad_rhoI3, const ThreeVector djI3_dt,
    const ThreeVector rot_jI3, const double rhoB, const ThreeVector grad_rhoB,
    const ThreeVector djB_dt, const ThreeVector rot_jB) const {
  ThreeVector E_component(0.0, 0.0, 0.0), B_component(0.0, 0.0, 0.0);
  if (use_symmetry_) {
    E_component -= dVsym_drhoI3(rhoB, rhoI3) * (grad_rhoI3 + djI3_dt) +
                   dVsym_drhoB(rhoB, rhoI3) * (grad_rhoB + djB_dt);
    B_component +=
        dVsym_drhoI3(rhoB, rhoI3) * rot_jI3 + dVsym_drhoB(rhoB, rhoI3) * rot_jB;
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
  auto F_skyrme =
      std::make_pair(ThreeVector(0., 0., 0.), ThreeVector(0., 0., 0.));
  auto F_symmetry =
      std::make_pair(ThreeVector(0., 0., 0.), ThreeVector(0., 0., 0.));

  const auto baryon_density_and_gradient = current_eckart(
      r, plist, param_, DensityType::Baryon, compute_gradient, smearing);
  const double rhoB = std::get<0>(baryon_density_and_gradient);
  const ThreeVector grad_rhoB = std::get<2>(baryon_density_and_gradient);
  const ThreeVector djB_dt = std::get<3>(baryon_density_and_gradient);
  const ThreeVector rot_jB = std::get<4>(baryon_density_and_gradient);
  if (use_skyrme_) {
    F_skyrme = skyrme_force(rhoB, grad_rhoB, djB_dt, rot_jB);
  }

  if (use_symmetry_) {
    const auto density_and_gradient =
        current_eckart(r, plist, param_, DensityType::BaryonicIsospin,
                       compute_gradient, smearing);
    const double rhoI3 = std::get<0>(density_and_gradient);
    const ThreeVector grad_rhoI3 = std::get<2>(density_and_gradient);
    const ThreeVector djI3_dt = std::get<3>(density_and_gradient);
    const ThreeVector rot_jI3 = std::get<4>(density_and_gradient);
    F_symmetry = symmetry_force(rhoI3, grad_rhoI3, djI3_dt, rot_jI3, rhoB,
                                grad_rhoB, djB_dt, rot_jB);
  }

  return std::make_tuple(F_skyrme.first, F_skyrme.second, F_symmetry.first,
                         F_symmetry.second);
}

}  // namespace smash
