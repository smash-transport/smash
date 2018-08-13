/*
 *
 *    Copyright (c) 2014-2018
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
   * Example: Configuring Potentials
   * --------------
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
   * \f[ U_{Sym} = \pm 2 S_{pot} \frac{\rho_n - \rho_p}{\rho_0} \,, \f]
   * where \f$ \rho_n\f$ is neutron density and \f$ \rho_p\f$ is proton
   * density. Definition and implementation are still to be worked out.
   *
   * \key S_pot (double, required, no default): \n
   *      Parameter \f$S_{pot}\f$ of symmetry potential in MeV
   */
  if (use_symmetry_) {
    symmetry_s_ = conf.take({"Symmetry", "S_Pot"});
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

double Potentials::symmetry_pot(const double baryon_isospin_density) const {
  return 1.0e-3 * 2. * symmetry_s_ * baryon_isospin_density / nuclear_density;
}

double Potentials::potential(const ThreeVector &r, const ParticleList &plist,
                             const ParticleType &acts_on) const {
  double total_potential = 0.0;
  const bool compute_gradient = false;
  const auto scale = force_scale(acts_on);

  if (!acts_on.is_baryon()) {
    return total_potential;
  }

  if (use_skyrme_) {
    const double rho_eck =
        rho_eckart(r, plist, param_, DensityType::Baryon, compute_gradient)
            .first;
    total_potential += scale.first * skyrme_pot(rho_eck);
  }
  if (use_symmetry_) {
    const double rho_iso =
        rho_eckart(r, plist, param_, DensityType::BaryonicIsospin,
                   compute_gradient)
            .first;
    const double sym_pot = symmetry_pot(rho_iso) * acts_on.isospin3_rel();
    total_potential += scale.second * sym_pot;
  }
  return total_potential;
}

std::pair<double, int> Potentials::force_scale(const ParticleType &data) const {
  double skyrme_scale = 1.0;
  if (data.pdgcode().is_hyperon()) {
    if (data.pdgcode().is_Xi()) {
      skyrme_scale = 1. / 3.;
    } else if (data.pdgcode().is_Omega()) {
      skyrme_scale = 0.;
    } else {
      skyrme_scale = 2. / 3.;
    }
  }
  skyrme_scale = skyrme_scale * data.pdgcode().baryon_number();
  // Symmetry force acts only on the neutron and proton.
  const int symmetry_scale =
      data.pdgcode().is_nucleon() ? data.pdgcode().baryon_number() : 0;
  return std::make_pair(skyrme_scale, symmetry_scale);
}

std::pair<ThreeVector, ThreeVector> Potentials::potential_gradient(
    const ThreeVector &r, const ParticleList &plist) const {
  const bool compute_gradient = true;
  const double MeV_to_GeV = 1.0e-3;
  ThreeVector dUB_dr(0.0, 0.0, 0.0);
  if (use_skyrme_) {
    const auto density_and_gradient =
        rho_eckart(r, plist, param_, DensityType::Baryon, compute_gradient);
    const double rho = density_and_gradient.first;
    const ThreeVector drho_dr = density_and_gradient.second;

    // Derivative of potential with respect to density
    double tmp = skyrme_tau_ * std::pow(rho / nuclear_density, skyrme_tau_ - 1);
    dUB_dr =
        MeV_to_GeV * drho_dr * (skyrme_a_ + skyrme_b_ * tmp) / nuclear_density;
  }

  ThreeVector dUsym_dr(0.0, 0.0, 0.0);
  if (use_symmetry_) {
    const ThreeVector p_iso =
        rho_eckart(r, plist, param_, DensityType::BaryonicIsospin,
                   compute_gradient)
            .second;
    dUsym_dr = MeV_to_GeV * 2. * symmetry_s_ * p_iso / nuclear_density;
  }
  return std::make_pair(dUB_dr, dUsym_dr);
}

}  // namespace smash
