/*
 *
 *    Copyright (c) 2014-2015
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
Potentials::Potentials(Configuration conf, const ExperimentParameters &param)
    : param_(param),
      use_skyrme_(conf.has_value({"Skyrme"})),
      use_symmetry_(conf.has_value({"Symmetry"})) {
  /*!\Userguide
   * \page input_potentials_ Potentials
   * To switch potentials on / off one can just uncomment/comment the
   * section, to switch on only Skyrme or Symmetry potentials uncomment
   * only the part you want to switch on.
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
   * \key Skyrme_A (float, required): \n
   *      Parameter A of Skyrme potential in MeV
   *
   * \key Skyrme_B (float, required): \n
   *      Parameter B of Skyrme potential in MeV
   *
   * \key Skyrme_Tau (float, required): \n
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
   * \key S_pot (float, required): \n
   *      Parameter \f$S_{pot}\f$ of symmetry potential in MeV
   */
  if (use_symmetry_) {
    symmetry_s_ = conf.take({"Symmetry", "S_Pot"});
  }
}

Potentials::~Potentials() {
}

double Potentials::skyrme_pot(const double baryon_density) const {
  const double tmp = baryon_density/nuclear_density;
  // Return in GeV
  return 1.0e-3 * (skyrme_a_*tmp + skyrme_b_*std::pow(tmp, skyrme_tau_));
}

double Potentials::symmetry_pot(const double baryon_isospin_density) const {
  // Return in GeV -> 10^-3 coefficient
  return 1.0e-3 * 2.*symmetry_s_ * baryon_isospin_density/nuclear_density;
}

double Potentials::potential(const ThreeVector &r,
                             const ParticleList &plist,
                             const PdgCode acts_on) const {
  double total_potential = 0.0;
  const bool compute_gradient = false;

  if (!acts_on.is_baryon()) {
    return total_potential;
  }

  if (use_skyrme_) {
    const double rho_eck = rho_eckart(r, plist, param_, DensityType::Baryon,
                                      compute_gradient).first;
    total_potential += skyrme_pot(rho_eck);
  }
  if (use_symmetry_) {
    // use isospin density
    const double rho_iso = rho_eckart(r, plist, param_,
                                      DensityType::BaryonicIsospin,
                                      compute_gradient).first;
    const double sym_pot = symmetry_pot(rho_iso) * acts_on.isospin3_rel();
    total_potential += sym_pot;
  }
  return total_potential;
}

ThreeVector Potentials::potential_gradient(const ThreeVector &r,
                                           const ParticleList &plist,
                                           const PdgCode acts_on) const {
  ThreeVector total_gradient(0.0, 0.0, 0.0);

  if (!acts_on.is_baryon()) {
    return total_gradient;
  }

  const bool compute_gradient = true;
  if (use_skyrme_) {
    const auto density_and_gradient = rho_eckart(r, plist, param_,
                           DensityType::Baryon, compute_gradient);
    const double rho = density_and_gradient.first;
    const ThreeVector drho_dr = density_and_gradient.second;

    // Derivative of potential with respect to density
    double tmp = skyrme_tau_ * std::pow(rho/nuclear_density, skyrme_tau_ - 1);
    total_gradient += drho_dr * (skyrme_a_ + skyrme_b_*tmp) / nuclear_density;
  }

  if (use_symmetry_) {
    // use isospin density
    const ThreeVector p_iso = rho_eckart(r, plist, param_,
                                         DensityType::BaryonicIsospin,
                                         compute_gradient).second;
    const ThreeVector dUsym_dr = 2.*symmetry_s_ * p_iso/nuclear_density
                                 * acts_on.isospin3_rel();
    total_gradient += dUsym_dr;
  }
  // Return in GeV
  return total_gradient * 1.0e-3;
}



}  // namespace Smash
