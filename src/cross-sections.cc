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
  p_lab_ = sqrt((mandelstam_s_ - squared_mass_a_ - squared_mass_b_)
                * (mandelstam_s_ - squared_mass_a_ - squared_mass_b_)
                - 4 * squared_mass_a_ * squared_mass_b_)
           / (2 * sqrt(squared_mass_b_));
}


float CrossSections::elastic(const ParticleData &data_a,
                             const ParticleData &data_b) {

  compute_kinematics(data_a, data_b);

  const PdgCode &pdg_a = data_a.type().pdgcode();
  const PdgCode &pdg_b = data_b.type().pdgcode();

  /* For now, the meson-meson and meson-baryon elastic cross sections
   * are simply given by the cross section parameter. */
  if (pdg_a.baryon_number() == 0 || pdg_b.baryon_number() == 0) {
    return elastic_parameter_;
  }

  /* For baryon-baryon, we have to check the parametrized cross sections */
  float sig;
  /* pp-scattering */
  if (pdg_a == pdg_b) {
    sig = pp_elastic(mandelstam_s_);
  /* ppbar-scattering */
  } else if (pdg_a.is_antiparticle_of(pdg_b)) {
    sig = ppbar_elastic(p_lab_);
  /* np-scattering */
  } else {
    sig = np_elastic(mandelstam_s_);
  }

  if (sig>0.) {
    return sig;
  } else {
    std::stringstream ss;
    ss << "problem in CrossSections::elastic: " << pdg_a.string().c_str()
       << " " << pdg_b.string().c_str() << " " << pdg_a.spin() << " "
       << pdg_b.spin() << " " << sig << " " << p_lab_ << " " << mandelstam_s_
       << " " << sqrt(squared_mass_a_);
    throw std::runtime_error(ss.str());
  }
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
