/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "include/decaytype.h"

#include <algorithm>
#include <stdexcept>
#include <gsl/gsl_integration.h>

#include "include/constants.h"
#include "include/width.h"
#include "include/resonances.h"
#include "include/logging.h"

namespace Smash {

// TwoBodyDecayStable

TwoBodyDecay::TwoBodyDecay(ParticleTypePtrList part_types, int l)
	                  : DecayType(part_types, l) {
  if (part_types.size() != 2 ) {
    throw std::runtime_error(
      std::string("Wrong number of particles in TwoBodyDecay constructor: ") +
      std::to_string(part_types.size()));
  }
}

int TwoBodyDecay::particle_number() const {
  return 2;
}

bool TwoBodyDecay::has_particles (const ParticleType &t_a, const ParticleType &t_b) const {
  return (*particle_types_[0] == t_a && *particle_types_[1] == t_b) ||
         (*particle_types_[0] == t_b && *particle_types_[1] == t_a);
}


// TwoBodyDecayStable

TwoBodyDecayStable::TwoBodyDecayStable(ParticleTypePtrList part_types, int l)
                                      : TwoBodyDecay(part_types, l) {
  if (!(part_types[0]->is_stable() && part_types[1]->is_stable())) {
    throw std::runtime_error(
      std::string("Error: Unstable particle in TwoBodyDecayStable constructor: ") +
      part_types[0]->pdgcode().string() + " " +
      part_types[1]->pdgcode().string());
  }
}

float TwoBodyDecayStable::rho (float m) const {
  // Determine momentum of outgoing particles in rest frame of resonance
  const float p_ab = pCM(m, particle_types_[0]->mass(),
                            particle_types_[1]->mass());
  // determine rho(m)
  return p_ab / m * BlattWeisskopf(p_ab*interaction_radius, L_);
}

float TwoBodyDecayStable::width (float m0, float G0, float m) const {
  if (m > particle_types_[0]->mass() + particle_types_[1]->mass())
    return G0 * rho(m) / rho(m0);
  else
    return 0;
}

float TwoBodyDecayStable::in_width (float m0, float G0, float m, float, float) const {
  // in-width = out-width
  return width(m0, G0, m);
}


// TwoBodyDecaySemistable

TwoBodyDecaySemistable::TwoBodyDecaySemistable(ParticleTypePtrList part_types,
                                               int l)
                                              : TwoBodyDecay(part_types, l) {
  // mark sure that the first particle is the stable one
  if (particle_types_[1]->is_stable()) {
    std::swap(particle_types_[0], particle_types_[1]);
  }
  // make sure that this is really a "semi-stable" decay
  if (!particle_types_[0]->is_stable() || particle_types_[1]->is_stable()) {
    throw std::runtime_error(
      std::string("Error in TwoBodyDecaySemistable constructor: ") +
      particle_types_[0]->pdgcode().string() + " " +
      particle_types_[1]->pdgcode().string());
  }
  Lambda_ = (particle_types_[1]->baryon_number()!=0) ? 2.0 : 1.6;
  // initialize tabulation
  init_tabulation(1., 50);
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

float TwoBodyDecaySemistable::calc_rho (float m) const {
  gsl_integration_workspace *workspace
    = gsl_integration_workspace_alloc(1000);
  IntegParam ip = {*particle_types_[1], particle_types_[0]->mass(), m, L_};
  gsl_function integrand;
  integrand.function = integrand_rho_Manley;
  integrand.params = &ip;
  size_t subintervals_max = 100;
  int gauss_points = 2;
  double accuracy_absolute = 1.0e-6;
  double accuracy_relative = 1.0e-4;
  double integral_value, integral_error;

  gsl_integration_qag(&integrand, ip.type.minimum_mass(), ip.srts - ip.m2,
                      accuracy_absolute, accuracy_relative,
                      subintervals_max, gauss_points, workspace,
                      &integral_value, &integral_error);

  gsl_integration_workspace_free(workspace);

  return integral_value;
}

void TwoBodyDecaySemistable::init_tabulation(float range, unsigned int N) {
  const auto &log = logger<LogArea::DecayType>();
  M_min_ = particle_types_[0]->mass()+particle_types_[1]->minimum_mass();
  log.info("Tabulating width for decay into ",
           particle_types_[0]->name(), " ",
           particle_types_[1]->name(), " (L=", L_,", M_min=", M_min_, ")");
  dM_ = range/N;
  tabulation_.resize(N);
  for (unsigned int i=0; i<N; i++) {
    tabulation_[i] = calc_rho(M_min_ + i*dM_);
  }
}

float TwoBodyDecaySemistable::rho (float m) const {
  if (m<M_min_) {
    return 0.;
  }
  // lookup tabulated values
  unsigned int n = static_cast<unsigned int>(round((m-M_min_)/dM_));
  if (n>=tabulation_.size()) {
    return tabulation_.back();
  } else {
    return tabulation_[n];
  }
}

float TwoBodyDecaySemistable::width (float m0, float G0, float m) const {
  return G0 * rho(m) / rho(m0)
         * Post_FF_sqr (m, m0, particle_types_[0]->mass()
                               + particle_types_[1]->minimum_mass(), Lambda_);
}

float TwoBodyDecaySemistable::in_width (float m0, float G0, float m, float m1, float m2) const {
  double p_f = pCM(m, m1, m2);

  return G0 * p_f * BlattWeisskopf(p_f*interaction_radius, L_)
         * Post_FF_sqr (m, m0, particle_types_[0]->mass()
                               + particle_types_[1]->minimum_mass(), Lambda_)
         / ( m * rho(m0) );
}


// TwoBodyDecayUnstable

TwoBodyDecayUnstable::TwoBodyDecayUnstable(ParticleTypePtrList part_types,
                                           int l)
                                          : TwoBodyDecay(part_types, l) {
  if (part_types[0]->is_stable() || part_types[1]->is_stable()) {
    throw std::runtime_error(
      std::string("Error: Stable particle in TwoBodyDecayUnstable constructor: ") +
      part_types[0]->pdgcode().string() + " " +
      part_types[1]->pdgcode().string());
  }
}

float TwoBodyDecayUnstable::rho (float) const {
  return 1.;
}

float TwoBodyDecayUnstable::width (float, float G0, float) const {
  return G0;  // use on-shell width
}

float TwoBodyDecayUnstable::in_width (float, float G0, float, float, float) const {
  return G0;  // use on-shell width
}


// ThreeBodyDecay

ThreeBodyDecay::ThreeBodyDecay(ParticleTypePtrList part_types, int l)
                               : DecayType(part_types, l) {
  if (part_types.size() != 3) {
    throw std::runtime_error(
      std::string("Wrong number of particles in ThreeBodyDecay constructor: ") +
      std::to_string(part_types.size()));
  }
}

int ThreeBodyDecay::particle_number() const {
  return 3;
}

bool ThreeBodyDecay::has_particles (const ParticleType &, const ParticleType &) const {
  return false;
}

float ThreeBodyDecay::rho (float) const {
  return 1.;
}

float ThreeBodyDecay::width (float, float G0, float) const {
  return G0;  // use on-shell width
}

float ThreeBodyDecay::in_width (float, float G0, float, float, float) const {
  return G0;  // use on-shell width
}


}
