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

#include "include/constants.h"
#include "include/processbranch.h"
#include "include/decaymodes.h"

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
  float s,mass12,x;
  s = srts*srts;
  mass12 = mass1*mass1;
  x = s+mass12-mass2*mass2;
  return sqrt(x*x/(4.*s)-mass12);
}


/**
 * Returns the value of the Blatt-Weisskopf functions,
 *  which govern the mass dependence of the width of a resonance
 * decaying into two particles A and B. See e.g. Effenberger's thesis, page 28.
 *
 * \param x = p_ab*R  with
 *        p_ab = relative momentum of outgoing particles AB and
 *        R = interaction radius
 * \param l Angular momentum of outgoing particles AB.
 */
static float BlattWeisskopf(const float x, const int L) {
  switch (L)
  {
    case 0:
      return 1.;
    case 1:
      return x / sqrt(1. + x*x);
    case 2:
      return x*x / sqrt(9. + 3.*x*x + x*x*x*x);
    case 3:
      return x*x*x / sqrt(225. + 45.*x*x + 6.*x*x*x*x + x*x*x*x*x*x);
    case 4:
      return x*x*x*x / sqrt(11025. + 1575.*x*x + 135.*x*x*x*x + 10.*x*x*x*x*x*x + x*x*x*x*x*x*x*x);
    default:
      throw std::invalid_argument(std::string("Wrong angular momentum in BlattWeisskopf: ")
                                  + std::to_string(L));
  }
}


/**
 * Get the mass-dependent width of a two-body decay into stable particles
 * according to Manley/Saleski, Phys. Rev. D 45 (1992) 4002.
 * 
 * \param mass Actual mass of the decaying particle [GeV].
 * \param poleMass Pole mass of the decaying particle [GeV].
 * \param mass1 Mass of the first daughter particle [GeV].
 * \param mass2 Mass of the second daughter particle [GeV].
 * \param L Angular momentum of the decay.
 * \param partialWidth_pole Partial width at the pole mass [GeV].
 */
static float width_Manley_stable(const float mass, const float poleMass,
                                 const float mass1, const float mass2,
                                 const int L, const float partialWidth_pole) {

  float bw, p_ab_mass, p_ab_pole, rho_ab_mass, rho_ab_pole;
  float interactionRadius = 1./hbarc;

  if (mass <= mass1 + mass2) {
    return 0.;
  }

  // Determine momentum of outgoing particles in Restframe of Resonance
  p_ab_mass = pCM(mass,mass1,mass2);
  p_ab_pole = pCM(poleMass,mass1,mass2);

  // Evaluate rho_ab according to equ. (2.76) in Effenberger's thesis
  // rho_ab(mu)=p_ab/mu * BlattWeisskopf(pab*interactionRadius,L)
  bw = BlattWeisskopf(p_ab_mass*interactionRadius,L);
  rho_ab_mass = p_ab_mass/mass * bw*bw;

  bw = BlattWeisskopf(p_ab_pole*interactionRadius,L);
  rho_ab_pole = p_ab_pole/poleMass * bw*bw;

  return partialWidth_pole * rho_ab_mass/rho_ab_pole;
}


float width_total(const ParticleType *t, const float m) {
  float w = 0., partial_width_at_pole;
  if (t->is_stable()) return w;
  const std::vector<DecayBranch> decaymodes
        = DecayModes::find(t->pdgcode()).decay_mode_list();
  // loop over decay modes and sum up all partial widths
  for (const auto &mode : decaymodes) {
    partial_width_at_pole = t->width_at_pole()*mode.weight();
    const ParticleType &t1 = ParticleType::find(mode.pdg_list()[0]);
    const ParticleType &t2 = ParticleType::find(mode.pdg_list()[1]);
    if (mode.pdg_list().size()==2 && t1.is_stable() && t2.is_stable()) {
      // mass-dependent width for 2-body decays
      w = w + width_Manley_stable(m, t->mass(), t1.mass(), t2.mass(),
                                  mode.angular_momentum(),
                                  partial_width_at_pole);
    }
    else {
      // constant width for three-body decays
      w = w + partial_width_at_pole;
    }
  }

  return w;
}


}  // namespace Smash
