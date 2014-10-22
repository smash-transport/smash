/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/forwarddeclarations.h"
#include "include/potentials.h"

namespace Smash {

/**
 * Potentials constructor. Gets parameters of potentials from configuration.
 */
Potentials::Potentials(Configuration conf)
    : use_skyrme_(conf.has_value({"Skyrme"})),
      use_symmetry_(conf.has_value({"Symmetry", "Enable"})) {
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
    skyrme_a_ = conf.take({"Skyrme","Skyrme_A"});
    skyrme_b_ = conf.take({"Skyrme","Skyrme_B"});
    skyrme_tau_ = conf.take({"Skyrme","Skyrme_Tau"});
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
    symmetry_s_ = conf.take({"Symmetry","S_pot"});
  }
}

Potentials::~Potentials() {
}

}  // namespace Smash
