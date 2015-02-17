/*
 *
 *    Copyright (c) 2012-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "include/width.h"

#include <cstdio>
#include <stdexcept>
#include <istream>
#include <gsl/gsl_integration.h>

#include "include/constants.h"
#include "include/resonances.h"

namespace Smash {


const float interaction_radius = 1. / hbarc;

/**
 * Returns the squared Blatt-Weisskopf functions,
 * which govern the mass dependence of the width of a resonance
 * decaying into two particles A and B. See e.g. Effenberger's thesis, page 28.
 *
 * \param x = p_ab*R  with
 *        p_ab = relative momentum of outgoing particles AB and
 *        R = interaction radius
 * \param L Angular momentum of outgoing particles AB.
 */
static float BlattWeisskopf(const float x, const int L)
#ifdef NDEBUG
    noexcept
#endif
{
  const auto x2 = x * x;
  switch (L) {
    case 0:
      return 1.f;
    case 1:
      return x2 / (1.f + x2);
    /* The following lines should be correct. But since nothing in SMASH uses
     * L > 1 this code is untested and dead. Therefore we only keep it as a
     * reference for later. (x4 = x2 * x2)
     * See also input sanitization in load_decaymodes in decaymodes.cc.
    case 2:
      return x4 / (9.f + 3.f * x2 + x4);
    case 3:
      return x4 * x2 / (225.f + 45.f * x2 + 6.f * x4 + x4 * x2);
    case 4:
      return x4 * x4 /
             (11025.f + 1575.f * x2 + 135.f * x4 + 10.f * x2 * x4 + x4 * x4);
    */
#ifndef NDEBUG
    default:
      throw std::invalid_argument(
          std::string("Wrong angular momentum in BlattWeisskopf: ") +
          std::to_string(L));
#endif
  }
  return 0.f;
}


float width_Manley_stable(const float mass, const float poleMass,
                          const float mass_a, const float mass_b,
                          const int L, const float partial_width_at_pole) {

  if (mass <= mass_a + mass_b) {
    return 0.;
  }

  // Determine momentum of outgoing particles in Restframe of Resonance
  const float p_ab_mass = pCM(mass, mass_a, mass_b);
  const float p_ab_pole = pCM(poleMass, mass_a, mass_b);

  // Evaluate rho_ab according to equ. (2.76) in Effenberger's thesis
  // rho_ab(mu)=p_ab/mu * BlattWeisskopf(pab*interaction_radius,L)
  const float rho_ab_mass = p_ab_mass / mass *
                            BlattWeisskopf(p_ab_mass*interaction_radius, L);

  const float rho_ab_pole = p_ab_pole / poleMass *
                            BlattWeisskopf(p_ab_pole*interaction_radius, L);

  return partial_width_at_pole * rho_ab_mass / rho_ab_pole;
}


/** Parameters for GSL integration. */
struct IntegParam {
  const ParticleType &type;  // type of daughter resonance
  double m2;                 // mass of stable particle
  double srts;               // sqrt(s) = mass of decaying resonance
  int L;                     // angular momentum
};


static double integrand_rho_Manley(double mass, void *parameters) {
  IntegParam *ip = reinterpret_cast<IntegParam*>(parameters);

  double stable_mass = ip->m2;
  double srts = ip->srts;

  if (srts < mass + stable_mass) {
    return 0.;
  }

  double resonance_width = ip->type.total_width(srts);
  /* center-of-mass momentum of final state particles */
  double p_f = pCM(srts, stable_mass, mass);

  return p_f/srts * BlattWeisskopf(p_f*interaction_radius, ip->L) * 2.*srts
          * spectral_function(mass, ip->type.mass(), resonance_width);
}


static double rho_Manley(IntegParam *ip) {
  gsl_integration_workspace *workspace
    = gsl_integration_workspace_alloc(1000);
  gsl_function integrand;
  integrand.function = integrand_rho_Manley;
  integrand.params = ip;
  size_t subintervals_max = 100;
  int gauss_points = 2;
  double accuracy_absolute = 1.0e-6;
  double accuracy_relative = 1.0e-4;
  double integral_value, integral_error;

  gsl_integration_qag(&integrand, ip->type.minimum_mass(), ip->srts - ip->m2,
                      accuracy_absolute, accuracy_relative,
                      subintervals_max, gauss_points, workspace,
                      &integral_value, &integral_error);

  gsl_integration_workspace_free(workspace);

  return integral_value;
}


static double Post_FF_sqr (double m, double M0, double s0, double L) {
  double FF = (L*L*L*L + (s0-M0*M0)*(s0-M0*M0)/4.) /
              (L*L*L*L + (m*m-(s0+M0*M0)/2.) * (m*m-(s0+M0*M0)/2.));
  return FF*FF;
}


float width_Manley_semistable(const float mass, const float poleMass,
                              const float mass_stable,
                              const ParticleType &type_unstable, const int L,
                              const float partial_width_at_pole) {
  IntegParam ip = {type_unstable, mass_stable, mass, L};

  double Lambda = (type_unstable.pdgcode().baryon_number()!=0) ? 2.0 : 1.6;

  // calculate rho function at off-shell mass
  double rho_off = rho_Manley(&ip);

  // calculate rho function at pole mass
  ip.srts = poleMass;
  double rho_pole = rho_Manley(&ip);

  return partial_width_at_pole * rho_off / rho_pole
         * Post_FF_sqr (mass, poleMass,
                        mass_stable+type_unstable.minimum_mass(), Lambda);
}


float in_width_Manley_semistable(const float mass, const float poleMass,
                                 const float mass_stable,
                                 const float mass_unstable,
                                 const ParticleType &type_unstable,
                                 const int L,
                                 const float partial_width_at_pole) {
  double p_f = pCM(mass, mass_stable, mass_unstable);
  IntegParam ip = {type_unstable, mass_stable, poleMass, L};

  double Lambda = (type_unstable.pdgcode().baryon_number()!=0) ? 2.0 : 1.6;

  return partial_width_at_pole * p_f
         * BlattWeisskopf(p_f*interaction_radius, L)
         * Post_FF_sqr (mass, poleMass,
                        mass_stable+type_unstable.minimum_mass(), Lambda)
         / ( mass * rho_Manley(&ip) );
}


}  // namespace Smash
