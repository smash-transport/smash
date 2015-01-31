/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
*/

#include "include/potentials.h"
#include "include/density.h"

namespace Smash {

/**
 * Potentials constructor. Gets parameters of potentials from configuration.
 */
Potentials::Potentials(Configuration conf, const ExperimentParameters &param)
    : use_skyrme_(conf.has_value({"Skyrme"})),
      use_symmetry_(conf.has_value({"Symmetry"})) {
  ntest_ = param.testparticles;
  sigma_ = param.gaussian_sigma;

  /*!\Userguide
   * \page potentials Potentials
   * Skyrme potential:
   * -----------------
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
   *      Parameter \f$\tau\f$ of Skyrme potent.
   */
  if (use_skyrme_) {
    skyrme_a_ = conf.take({"Skyrme", "Skyrme_A"});
    skyrme_b_ = conf.take({"Skyrme", "Skyrme_B"});
    skyrme_tau_ = conf.take({"Skyrme", "Skyrme_Tau"});
  }

  /*!\Userguide
   * \page potentials Potentials
   * Symmetry potential:
   * -------------------
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

double Potentials::potential(const ThreeVector &r,
                             const ParticleList &plist,
                             const PdgCode acts_on) const {
  double total_potential = 0.0;

  if (!acts_on.is_baryon()) {
    return total_potential;
  }

  if (use_skyrme_) {
    const double rho_eckart = four_current(r, plist, sigma_,
                                           baryon_density, ntest_).abs();
    total_potential += skyrme_a_ * (rho_eckart/rho0) +
                       skyrme_b_ * std::pow(rho_eckart/rho0, skyrme_tau_);
  }
  if (use_symmetry_) {
    // use isospin density
    const double rho_iso = four_current(r, plist, sigma_,
                                        isospin_density, ntest_).abs();
    const double sym_pot = 2.*symmetry_s_ * rho_iso/rho0
                           * float(acts_on.isospin3())/acts_on.isospin_total();
    total_potential += sym_pot;
  }
  // Return in GeV
  return total_potential * 1.0e-3;
}

ThreeVector Potentials::potential_gradient(const ThreeVector &r,
                                           const ParticleList &plist,
                                           const PdgCode acts_on) const {
  ThreeVector total_gradient(0.0, 0.0, 0.0);

  if (!acts_on.is_baryon()) {
    return total_gradient;
  }

  if (use_skyrme_) {
    const auto density_and_gradient = rho_eckart_gradient(r, plist,
                                                 sigma_, baryon_density, ntest_);
    const double rho = density_and_gradient.first;
    const ThreeVector drho_dr = density_and_gradient.second;

    // Derivative of potential with respect to density
    double tmp = skyrme_tau_ * std::pow(rho/rho0, skyrme_tau_ - 1);
    const double dpotential_drho = (skyrme_a_  + skyrme_b_ * tmp) / rho0;
    total_gradient += drho_dr * dpotential_drho;
  }

  if (use_symmetry_) {
    // use isospin density
    const ThreeVector p_iso = rho_eckart_gradient(r, plist, sigma_,
                                                  isospin_density, ntest_).second;
    const ThreeVector dUsym_dr = 2.*symmetry_s_ * p_iso/rho0
                                 * float(acts_on.isospin3())/acts_on.isospin_total();
    total_gradient += dUsym_dr;
  }
  // Return in GeV
  return total_gradient * 1.0e-3;
}



}  // namespace Smash
