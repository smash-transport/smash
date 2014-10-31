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
      use_symmetry_(conf.has_value({"Symmetry", "Enable"})) {
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
    symmetry_s_ = conf.take({"Symmetry", "S_pot"});
  }
}

Potentials::~Potentials() {
}

double Potentials::potential(const ThreeVector &r,
                             const ParticleList &plist) const {
  double total_potential = 0.0;

  if (use_skyrme_) {
    const Density_type dens_type = baryon;
    const double rho_eckart = four_current(r, plist, sigma_,
                                           dens_type, ntest_).abs();
    total_potential += skyrme_a * (rho_eckart/rho0) +
                       skyrme_b * std::pow(rho_eckart/rho0, skyrme_tau);
  }
  if (use_symmetry_) {
    // TODO(oliiny): use neutron-proton density here or isospin density?
    // total_potential +=
  }
  return total_potential;
}

ThreeVector Potentials::potential_gradient(const ThreeVector &r,
                                           const ParticleList &plist) const {
  ThreeVector total_gradient(0.0, 0.0, 0.0);
  double tmp;

  if (use_skyrme_) {
    const Density_type dens_type = baryon;
    const auto density_and_gradient = rho_eckart_gradient(r, plist,
                                                 sigma_, dens_type, ntest_);
    const double rho = density_and_gradient.first;
    const ThreeVector drho_dr = density_and_gradient.second;

    // Derivative of potential with respect to density
    tmp = skyrme_tau * std::pow(rho/rho0, skyrme_tau - 1);
    const double dpotential_drho = (skyrme_a  + skyrme_b * tmp) / rho0;
    total_gradient += drho_dr * dpotential_drho;
  }

  if (use_symmetry_) {
    // TODO(oliiny): use neutron-proton density here or isospin density?
    // total_gradient +=
  }
  return total_gradient;
}



}  // namespace Smash
