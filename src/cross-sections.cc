/*
 *
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "include/CrossSections.h"
#include "include/constants.h"
#include "include/FourVector.h"
#include "include/outputroutines.h"
#include "include/parametrizations.h"
#include "include/Particles.h"
#include "include/ParticleData.h"
#include "include/ParticleType.h"

void CrossSections::compute_kinematics(Particles *particles,
  int id_a, int id_b) {
  /* Mandelstam s = (p_a + p_b)^2 = square of CMS energy */
  mandelstam_s_ =
    (particles->data(id_a).momentum() + particles->data(id_b).momentum()).Dot(
     particles->data(id_a).momentum() + particles->data(id_b).momentum());
  /* In case we have resonances, let's use the relation m^2 = p^2 for masses */
  squared_mass_a_
    = particles->data(id_a).momentum().Dot(particles->data(id_a).momentum());
  squared_mass_b_
    = particles->data(id_b).momentum().Dot(particles->data(id_b).momentum());
  /* Beam momentum (assuming particle A as "beam") */
  p_lab_ = sqrt((mandelstam_s_ - squared_mass_a_ - squared_mass_b_)
                / (2 * sqrt(squared_mass_b_)));
}

float CrossSections::parametrization_(std::vector<float> parameters) const {
  switch (parameters.size() - 1) {
  case 1:
    /* just a constant */
    return parameters[1];
    break;
  case 2:
    /* A + B*mass_a/(mandelstam s - 2*mass_a^2 - 2*mass_b^2) */
    return parameters[1] + parameters[2] * sqrt(squared_mass_a_)
      / (mandelstam_s_ - 2 * squared_mass_a_ - 2 * squared_mass_b_);
    break;
  case 3:
    /* A * p^n * exp[-B*(log(p))^2] */
    return parameters[1] * pow(p_lab_, parameters[2])
      * exp(-parameters[3] * log(p_lab_) * log(p_lab_));
    break;
  case 4:
    /* A + B*(p-b)*exp(-C*p) */
    return parameters[1] + parameters[2] * (p_lab_ - parameters[3])
      * exp(-parameters[4] * p_lab_);
    break;
  case 5:
    /* A + B*p^n + C*(log(p))^2 + D*log(p) */
    return parameters[1] + parameters[2] * pow(p_lab_, parameters[3])
      + parameters[4] * log(p_lab_) * log(p_lab_) + parameters[5] * log(p_lab_);
    break;
  case 7:
    /* A + B*(abs(p+b))^k + C*(p+c)^m */
    return parameters[1]
      + parameters[2] * pow(fabs(p_lab_ + parameters[3]), parameters[4])
      + parameters[5] * pow(p_lab_ + parameters[6], parameters[7]);
  case 8:
    /* A + B*{C + c*exp[-D*(p+d)^m]}^n */
    return parameters[1] + parameters[2]
      * pow(parameters[3] + parameters[4] * exp(-parameters[5]
            * pow(p_lab_ + parameters[6], parameters[7])), parameters[8]);
    break;
  case 10:
    /* A + B*(abs(p+b))^k + C*(p+c)^m + D*(p+d)^n */
    return parameters[1]
      + parameters[2] * pow(fabs(p_lab_ + parameters[3]), parameters[4])
      + parameters[5] * pow(p_lab_ + parameters[6], parameters[7])
      + parameters[8] * pow(p_lab_ + parameters[9], parameters[10]);
    break;
  default:
    return 0.0;
  }
}

float CrossSections::elastic(Particles *particles, int id_a, int id_b)
  const {
  const int spin_a = particles->type(id_a).spin(),
    spin_b = particles->type(id_b).spin();
  /* For now, the meson-meson and meson-baryon elastic cross sections
   * are simply given by the cross section parameter
   */
  if (spin_a % 2 == 0 || spin_b % 2 == 0)
    return elastic_parameter_;

  /* For baryon-baryon, we have to check the parametrized cross sections */
  /* pp-scattering */
  if (particles->type(id_a).pdgcode() == particles->type(id_b).pdgcode()) {
    return pp_elastic(p_lab_, mandelstam_s_, sqrt(squared_mass_a_));
  /* ppbar-scattering */
  } else if (abs(particles->type(id_a).pdgcode())
          == abs(particles->type(id_b).pdgcode())) {
    return ppbar_elastic(p_lab_);
  /* np-scattering */
  } else {
    return np_elastic(p_lab_, mandelstam_s_, sqrt(squared_mass_a_));
  }
}

float CrossSections::annihilation(Particles *particles,
                                  int id_a, int id_b) const {
  const int spin_a = particles->type(id_a).spin(),
    spin_b = particles->type(id_b).spin();
  /* For now, the meson-meson and meson-baryon
   * annihilation cross sections are zero
   */
  if (spin_a % 2 == 0 || spin_b % 2 == 0)
    return 0.0;

  /* For baryon-baryon, we have to check the parametrized cross sections */
  if (particles->type(id_a).pdgcode() == -(particles->type(id_b).pdgcode())) {
    size_t i = 0;
    while (p_lab_ < (ppbar_annihilation_[i])[0])
      i++;
    return parametrization_(ppbar_annihilation_[i]);
  } else {
    return 0.0;
  }
}

float CrossSections::total(Particles *particles, int id_a, int id_b)
  const {
  const int spin_a = particles->type(id_a).spin(),
    spin_b = particles->type(id_b).spin();
  /* For now, the meson-meson and meson-baryon
   *  "total" cross section is just zero
   */
  if (spin_a % 2 == 0 || spin_b % 2 == 0)
    return 0.0;

  /* For baryon-baryon, we have to check the parametrized cross sections */
  /* pp-scattering */
  if (particles->type(id_a).pdgcode() == particles->type(id_b).pdgcode()) {
    return pp_total(p_lab_);
  /* ppbar-scattering */
  } else if (abs(particles->type(id_a).pdgcode())
          == abs(particles->type(id_b).pdgcode())) {
    return ppbar_total(p_lab_);
  /* np-scattering */
  } else {
    return np_total(p_lab_);
  }
}
