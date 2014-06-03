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

void CrossSections::compute_kinematics(const Particles &particles, int id_a,
                                       int id_b) {
  const FourVector momentum_a = particles.data(id_a).momentum();
  const FourVector momentum_b = particles.data(id_b).momentum();
  /* Mandelstam s = (p_a + p_b)^2 = square of CMS energy */
  mandelstam_s_ = (momentum_a + momentum_b).Dot(momentum_a + momentum_b);
  /* In case we have resonances, let's use the relation m^2 = p^2 for masses */
  squared_mass_a_ = momentum_a.Dot(momentum_a);
  squared_mass_b_ = momentum_b.Dot(momentum_b);
  /* Beam momentum (assuming particle A as "beam") */
  p_lab_ = sqrt((mandelstam_s_ - squared_mass_a_ - squared_mass_b_)
                * (mandelstam_s_ - squared_mass_a_ - squared_mass_b_)
                - 4 * squared_mass_a_ * squared_mass_b_)
           / (2 * sqrt(squared_mass_b_));
}

float CrossSections::elastic(const Particles &particles, int id_a,
                             int id_b) const {
  const ParticleType &type_a = particles.type(id_a);
  const ParticleType &type_b = particles.type(id_b);

  const int spin_a = type_a.spin();
  const int spin_b = type_b.spin();
  /* For now, the meson-meson and meson-baryon elastic cross sections
   * are simply given by the cross section parameter
   */
  if (spin_a % 2 == 0 || spin_b % 2 == 0) {
    return elastic_parameter_;
  }

  /* For baryon-baryon, we have to check the parametrized cross sections */
  /* pp-scattering */
  if (type_a.pdgcode() == type_b.pdgcode()) {
    return pp_elastic(p_lab_, mandelstam_s_, sqrt(squared_mass_a_));
  /* ppbar-scattering */
  } else if (type_a.pdgcode().is_antiparticle_of(type_b.pdgcode())) {
    return ppbar_elastic(p_lab_);
  /* np-scattering */
  } else {
    return np_elastic(p_lab_, mandelstam_s_, sqrt(squared_mass_a_));
  }
}

float CrossSections::total(const Particles &particles, int id_a,
                           int id_b) const {
  const ParticleType &type_a = particles.type(id_a);
  const ParticleType &type_b = particles.type(id_b);

  const int spin_a = type_a.spin();
  const int spin_b = type_b.spin();
  /* For now, the meson-meson and meson-baryon
   *  "total" cross section is just zero
   */
  if (spin_a % 2 == 0 || spin_b % 2 == 0) {
    return 0.f;
  }

  /* For baryon-baryon, we have to check the parametrized cross sections */
  /* pp-scattering */
  if (type_a.pdgcode() == type_b.pdgcode()) {
    return pp_total(p_lab_);
  /* ppbar-scattering */
  } else if (type_a.pdgcode().is_antiparticle_of(type_b.pdgcode())) {
    return ppbar_total(p_lab_);
  /* np-scattering */
  } else {
    return np_total(p_lab_);
  }
}

}  // namespace Smash
