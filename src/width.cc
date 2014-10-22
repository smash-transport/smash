/*
 *
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "include/width.h"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <istream>
#include <gsl/gsl_integration.h>

#include "include/constants.h"
#include "include/resonances.h"

namespace Smash {


/**
 * Return the center-of-mass momentum of two particles,
 * given sqrt(s) and their masses.
 *
 * \param srts sqrt(s) of the process [GeV].
 * \param mass1 Mass of first particle [GeV].
 * \param mass2 Mass of second particle [GeV].
 */
static float pCM(const float srts, const float mass1, const float mass2) {
  float s, mass12, x;
  s = srts*srts;
  mass12 = mass1*mass1;
  x = s + mass12 - mass2*mass2;
  return std::sqrt(x*x / (4. * s) - mass12);
}


const float interactionRadius = 1. / hbarc;

/**
 * Returns the value of the Blatt-Weisskopf functions,
 *  which govern the mass dependence of the width of a resonance
 * decaying into two particles A and B. See e.g. Effenberger's thesis, page 28.
 *
 * \param x = p_ab*R  with
 *        p_ab = relative momentum of outgoing particles AB and
 *        R = interaction radius
 * \param L Angular momentum of outgoing particles AB.
 */
static float BlattWeisskopf(const float x, const int L) {
  switch (L) {
    case 0:
      return 1.;
    case 1:
      return x / std::sqrt(1. + x*x);
    case 2:
      return x*x / std::sqrt(9. + 3. * x*x + x*x*x*x);
    case 3:
      return x*x*x / std::sqrt(225. + 45. * x*x + 6.*x*x*x*x + x*x*x*x*x*x);
    case 4:
      return x*x*x*x / std::sqrt(11025. + 1575. * x*x + 135. * x*x*x*x
                            + 10. * x*x*x*x*x*x + x*x*x*x*x*x*x*x);
    default:
      throw std::invalid_argument(
        std::string("Wrong angular momentum in BlattWeisskopf: ")
        + std::to_string(L));
  }
}


float width_Manley_stable(const float mass, const float poleMass,
                          const float mass1, const float mass2,
                          const int L, const float partialWidth_pole) {

  if (mass <= mass1 + mass2) {
    return 0.;
  }

  // Determine momentum of outgoing particles in Restframe of Resonance
  const float p_ab_mass = pCM(mass, mass1, mass2);
  const float p_ab_pole = pCM(poleMass, mass1, mass2);

  // Evaluate rho_ab according to equ. (2.76) in Effenberger's thesis
  // rho_ab(mu)=p_ab/mu * BlattWeisskopf(pab*interactionRadius,L)
  float bw = BlattWeisskopf(p_ab_mass * interactionRadius, L);
  const float rho_ab_mass = p_ab_mass / mass * bw*bw;

  bw = BlattWeisskopf(p_ab_pole*interactionRadius, L);
  const float rho_ab_pole = p_ab_pole / poleMass * bw*bw;

  return partialWidth_pole * rho_ab_mass / rho_ab_pole;
}


/** Parameters for GSL integration. */
struct IntegParam {
  const ParticleType *type;  // type of daughter resonance
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

  double mandelstam_s = srts*srts;
  double resonance_width = ip->type->total_width(srts);
  /* center-of-mass momentum of final state particles */
  double p_f
    = std::sqrt((mandelstam_s - (stable_mass + mass)
            * (stable_mass + mass))
            * (mandelstam_s - (stable_mass - mass)
              * (stable_mass - mass))
            / (4 * mandelstam_s));

  return p_f/srts * BlattWeisskopf(p_f*interactionRadius, ip->L) * 2.*srts
          * spectral_function(mass, ip->type->mass(), resonance_width);
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

  gsl_integration_qag(&integrand, ip->type->minimum_mass(), ip->srts - ip->m2,
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
                              const ParticleType *type_unstable,
                              const int L, const float partialWidth_pole) {
  IntegParam ip;
  ip.type = type_unstable;
  ip.m2 = mass_stable;
  ip.L = L;

  double Lambda = (type_unstable->pdgcode().baryon_number()!=0) ? 2.0 : 1.6;

  // calculate rho function at off-shell mass
  ip.srts = mass;
  double rho_off = rho_Manley(&ip);

  // calculate rho function at pole mass
  ip.srts = poleMass;
  double rho_pole = rho_Manley(&ip);

  return partialWidth_pole * rho_off / rho_pole
         * Post_FF_sqr (mass, poleMass,
                        mass_stable+type_unstable->minimum_mass(), Lambda);
}

}  // namespace Smash
