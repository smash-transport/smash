/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include <cmath>
#include <cstdio>

#include "include/constants.h"
#include "include/crosssections.h"
#include "include/fourvector.h"
#include "include/outputroutines.h"
#include "include/parametrizations.h"
#include "include/particledata.h"
#include "include/particles.h"
#include "include/particletype.h"

namespace Smash {

void CrossSections::compute_kinematics(const ParticleData &data_a,
                                       const ParticleData &data_b) {
  const FourVector momentum_a = data_a.momentum();
  const FourVector momentum_b = data_b.momentum();
  /* Mandelstam s = (p_a + p_b)^2 = square of CMS energy */
  mandelstam_s_ = (momentum_a + momentum_b).sqr();
  /* In case we have resonances, let's use the relation m^2 = p^2 for masses */
  squared_mass_a_ = momentum_a.sqr();
  squared_mass_b_ = momentum_b.sqr();
  /* Beam momentum (assuming particle A as "beam") */
  p_lab_ = std::sqrt((mandelstam_s_ - squared_mass_a_ - squared_mass_b_)
                     * (mandelstam_s_ - squared_mass_a_ - squared_mass_b_)
                     - 4 * squared_mass_a_ * squared_mass_b_)
           / (2 * std::sqrt(squared_mass_b_));
}


float CrossSections::total(const PdgCode &pdg_a, const PdgCode &pdg_b) const {
  /* For now, the meson-meson and meson-baryon
   *  "total" cross section is just zero. */
  if (pdg_a.baryon_number() == 0 || pdg_b.baryon_number() == 0) {
    return 0.f;
  }

  /* For baryon-baryon, we have to check the parametrized cross sections */
  /* pp-scattering */
  if (pdg_a == pdg_b) {
    return pp_total(p_lab_);
  /* ppbar-scattering */
  } else if (pdg_a.is_antiparticle_of(pdg_b)) {
    return ppbar_total(p_lab_);
  /* np-scattering */
  } else {
    return np_total(p_lab_);
  }
}

}  // namespace Smash
