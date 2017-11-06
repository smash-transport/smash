/*
 *
 *    Copyright (c) 2014-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
*/

#include "include/potentials.h"

#include "include/constants.h"
#include "include/density.h"

namespace Smash {

/**
 * Potentials constructor. Gets parameters of potentials from configuration.
 */
Potentials::Potentials(Configuration conf, const DensityParameters &param)
    : param_(param),
      use_skyrme_(conf.has_value({"Skyrme"})),
      use_symmetry_(conf.has_value({"Symmetry"})) {
  /*!\Userguide
   * \page input_potentials_ Potentials
   * To switch potentials on / off one can just uncomment/comment the
   * section, to switch on only Skyrme or Symmetry potentials uncomment
   * only the part you want to switch on.
   *
   *
   * Currently potentials are just added without re-adjusting the energy
   * and momenta of colliding nucleons. This can be done, because the
   * binding energy of nucleons is from 0 to 8 MeV per nucleon and
   * kinetic energies, at which SMASH operates are at least 400 MeV per
   * nucleon, so the binding energy can be neglected compared to the
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
   * \key Skyrme_A (double, required): \n
   *      Parameter A of Skyrme potential in MeV
   *
   * \key Skyrme_B (double, required): \n
   *      Parameter B of Skyrme potential in MeV
   *
   * \key Skyrme_Tau (double, required): \n
   *      Parameter \f$\tau\f$ of Skyrme potential.
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
   * \key S_pot (double, required): \n
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
     under a charge reversal transformation. */
  const int sgn = tmp > 0 ? 1 : -1;
  // Return in GeV
  return 1.0e-3 * sgn * (skyrme_a_ * std::abs(tmp) +
                         skyrme_b_ * std::pow(std::abs(tmp), skyrme_tau_));
}

double Potentials::symmetry_pot(const double baryon_isospin_density) const {
  // Return in GeV -> 10^-3 coefficient
  return 1.0e-3 * 2. * symmetry_s_ * baryon_isospin_density / nuclear_density;
}

double Potentials::potential(const ThreeVector &r, const ParticleList &plist,
                             const ParticleType &acts_on) const {
  double total_potential = 0.0;
  const bool compute_gradient = false;

  if (!acts_on.is_baryon()) {
    return total_potential;
  }

  if (use_skyrme_) {
    const double rho_eck =
        rho_eckart(r, plist, param_, DensityType::Baryon, compute_gradient)
            .first;
    total_potential += skyrme_pot(rho_eck);
  }
  if (use_symmetry_) {
    // use isospin density
    const double rho_iso =
        rho_eckart(r, plist, param_, DensityType::BaryonicIsospin,
                   compute_gradient)
            .first;
    const double sym_pot = symmetry_pot(rho_iso) * acts_on.isospin3_rel();
    total_potential += sym_pot;
  }
  return total_potential;
}

std::pair<double, int> Potentials::force_scale(const ParticleType &data) const {
  /* For Lambda and Sigma, since they carry 2 light (u or d) quarks, they
  are affected by 2/3 of the Skyrme force. Xi carries 1 light quark, it
  is affected by 1/3 of the Skyrme force. Omega carries no light quark, so
  it's not affected by the Skyrme force. Anti-baryons are affected by the
  force as large as the force acting on baryons but with an opposite
  direction.*/
  double skyrme_scale = 1.0;
  if (data.pdgcode().is_hyperon()) {
    if (data.pdgcode().is_xi1321()) {
      skyrme_scale = 1. / 3.;
    } else if (data.pdgcode().is_Omega1672()) {
      skyrme_scale = 0.;
    } else {
      skyrme_scale = 2. / 3.;
    }
  }
  skyrme_scale = skyrme_scale * data.pdgcode().baryon_number();
  /* Hyperons are not affected by the symmetry force.*/
  const int symmetry_scale =
      data.pdgcode().is_hyperon() ? 0 : data.pdgcode().baryon_number();
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
    // use isospin density
    const ThreeVector p_iso =
        rho_eckart(r, plist, param_, DensityType::BaryonicIsospin,
                   compute_gradient)
            .second;
    dUsym_dr = MeV_to_GeV * 2. * symmetry_s_ * p_iso / nuclear_density;
  }
  return std::make_pair(dUB_dr, dUsym_dr);
}

}  // namespace Smash
