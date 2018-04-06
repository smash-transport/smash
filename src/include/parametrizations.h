/*
 *
 *    Copyright (c) 2013-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#ifndef SRC_INCLUDE_PARAMETRIZATIONS_H_
#define SRC_INCLUDE_PARAMETRIZATIONS_H_

#include <unordered_map>
#include <utility>

#include "particletype.h"

/* All quantities in this file use they same units as the rest of SMASH.
 * That is: GeV for energies and momenta, fm for distances and time, and mb for
 * cross-sections. */

namespace smash {

/**
 * total hadronic cross sections at high energies parametrized in the 2016 PDG
 * book(http://pdg.lbl.gov/2016/reviews/rpp2016-rev-cross-section-plots.pdf)
 *
 * This function is a utility function called from specific parametrizations.
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \param[in] is_opposite_charge whether the particles being collided have
 *                               opposite charges
 * \param[in] ma mass of first particle
 * \param[in] mb mass of second particle
 * \param[in] P Pomeranchuk's constant term
 * \param[in] R1 intensity of the first Regge pole contribution
 * \param[in] R2 intensity of the second Regge pole contribution
 * \return the parametrized cross-section
 */
double xs_high_energy(double mandelstam_s, bool is_opposite_charge, double ma,
                      double mb, double P, double R1, double R2);

/**
 * pp total cross section at high energies
 * \see xs_high_energy
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double pp_high_energy(double mandelstam_s);

/**
 * ppbar total cross section at high energies
 * \see xs_high_energy
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double ppbar_high_energy(double mandelstam_s);

/**
 * np total cross section at high energies
 * \see xs_high_energy
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double np_high_energy(double mandelstam_s);

/**
 * npbar total cross section at high energies
 * \see xs_high_energy
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double npbar_high_energy(double mandelstam_s);

/**
 * pi+p total cross section at high energies
 * \see xs_high_energy
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double piplusp_high_energy(double mandelstam_s);

/**
 * pi-p total cross section at high energies
 * \see xs_high_energy
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double piminusp_high_energy(double mandelstam_s);

/**
 * Utility function called by specific other parametrizations
 * Parametrized hard scattering cross section (with partonic scattering)
 * This parametrization is a direct fit to cross sections in PYTHIA
 * See \iref{Sjostrand:1987su}
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \param[in] xs_0 a fit parameter
 * \param[in] e_0 a fit parameter
 * \param[in] lambda_pow a fit parameter
 * \return the parametrized cross-section
 */
double xs_string_hard(double mandelstam_s, double xs_0, double e_0,
                      double lambda_pow);

/**
 * nucleon-nucleon hard scattering cross section (with partonic scattering)
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double NN_string_hard(double mandelstam_s);

/**
 * nucleon-pion hard scattering cross section (with partonic scattering)
 *
 * \param[in] mandelstam_s the rest frame total energy squared 
 * \return the parametrized cross-section
 */
double Npi_string_hard(double mandelstam_s);

/**
 * pion-pion hard scattering cross section (with partonic scattering)
 *
 * \param[in] mandelstam_s the rest frame total energy squared 
 * \return the parametrized cross-section
 */
double pipi_string_hard(double mandelstam_s);

/**
 * pi+p elastic cross section parametrization, PDG data.
 * Source: GiBUU:parametrizationBarMes_HighEnergy.f90
 *
 * The parametrizations of the elastic pion+nucleon cross sections
 * are still under tuning. The parametrizaton is employed to give a
 * non-zero cross section at high energies. To make sure it
 * doesn't affect the cross section at the low energies, I truncate
 * the parametrization at p_lab = 8 GeV, which correspons to square
 * root of s equal to 4 GeV.
 *
 * \param[in] mandelstam_s the rest frame total energy squared 
 * \return the parametrized cross-section
 */
double piplusp_elastic(double mandelstam_s);

/**
 * pi-p elastic cross section parametrization
 * Source: GiBUU:parametrizationBarMes_HighEnergy.f90
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double piminusp_elastic(double mandelstam_s);

/**
 * pp elastic cross section parametrization
 * Source: \iref{Weil:2013mya}, eq. (44)
 *
 * \param[in] mandelstam_s the rest frame total energy squared 
 * \return the parametrized cross-section
 */
double pp_elastic(double mandelstam_s);

/**
 * pp total cross section parametrization
 * Sources:
 * low-p: \iref{Cugnon:1996kh}
 * highest-p: \iref{Buss:2011mx}
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double pp_total(double mandelstam_s);

/**
 * np elastic cross section parametrization
 * Source: \iref{Weil:2013mya}, eq. (45)
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double np_elastic(double mandelstam_s);

/**
 * np total cross section parametrization
 * Sources:
 * low-p: \iref{Cugnon:1996kh}
 * highest-p: \iref{Buss:2011mx}
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double np_total(double mandelstam_s);

/**
 * ppbar elastic cross section parametrization
 * Source: \iref{Bass:1998ca}
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double ppbar_elastic(double mandelstam_s);

/**
 * ppbar total cross section parametrization
 * Source: \iref{Bass:1998ca}
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double ppbar_total(double mandelstam_s);

/**
 * ppbar total cross section parametrization;
 * Used for reverse cross-section from detailed balance (JB: no it's not ಠ_ಠ )
 * \todo unused (also not even defined)
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \param[in] m_proj mass of projectile?
 * \param[in] m_target mass of target?
 * \return the parametrized cross-section
 */
double ppbar_total(double mandelstam_s, double m_proj, double m_target);

/**
 * K+ p elastic background cross section parametrization.
 * sigma(K+n->K+n) = sigma(K+n->K0p) = 0.5 * sigma(K+p->K+p)
 * Source: \iref{Buss:2011mx}, B.3.8
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double kplusp_elastic_background(double mandelstam_s);

/**
 * K+ n elastic background cross section parametrization
 * sigma(K+n->K+n) = sigma(K+n->K0p) = 0.5 * sigma(K+p->K+p)
 * Source: \iref{Buss:2011mx}, B.3.8
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double kplusn_elastic_background(double mandelstam_s);

/**
 * K+ n charge exchange cross section parametrization.
 * sigma(K+n->K+n) = sigma(K+n->K0p) = 0.5 * sigma(K+p->K+p)
 * Source: \iref{Buss:2011mx}, B.3.8
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double kplusn_k0p(double mandelstam_s);

/**
 * K- p elastic background cross section parametrization
 * Source: \iref{Buss:2011mx}, B.3.9
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double kminusp_elastic_background(double mandelstam_s);

/**
 * K- n elastic background cross section parametrization
 * Source: \iref{Buss:2011mx}, B.3.9
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double kminusn_elastic_background(double mandelstam_s);

/**
 * K0 p elastic background cross section parametrization
 * Source: \iref{Buss:2011mx}, B.3.9
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double k0p_elastic_background(double mandelstam_s);

/**
 * K0 n elastic background cross section parametrization
 * Source: \iref{Buss:2011mx}, B.3.9
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double k0n_elastic_background(double mandelstam_s);

/**
 * Kbar0 p elastic background cross section parametrization
 * Source: \iref{Buss:2011mx}, B.3.9
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double kbar0p_elastic_background(double mandelstam_s);

/**
 * Kbar0 n elastic background cross section parametrization
 * Source: \iref{Buss:2011mx}, B.3.9
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double kbar0n_elastic_background(double mandelstam_s);

/**
 * K+ p inelastic background cross section parametrization
 * Source: \iref{Buss:2011mx}, B.3.8
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double kplusp_inelastic_background(double mandelstam_s);

/**
 * K+ n inelastic background cross section parametrization
 * Source: \iref{Buss:2011mx}, B.3.8
 *
 * This interpolates the experimental data of the total cross section and
 * subtracts the elastic and charge exchange cross section.
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double kplusn_inelastic_background(double mandelstam_s);

/**
 * Hash a pair of integers.
 *
 * Note that symmetric pairs and permutations yield identical hashes with this
 * implementation.
 */
struct pair_hash {
  std::size_t operator()(const std::pair<uint64_t, uint64_t>& p) const {
    auto h1 = std::hash<uint64_t>{}(p.first);
    auto h2 = std::hash<uint64_t>{}(p.second);

    // In our case the integers are PDG codes. We know they are different
    // and their order is defined, so we can simply combine the hashes
    // using XOR. Note that this yields 0 for h1 == h2. Also,
    // std::swap(h1, h2) does not not change the final hash.
    assert(h1 != h2);
    return h1 ^ h2;
  }
};

/**
 * Isospin weights for inelastic K+ N channels.
 */
class KplusNRatios {
 private:
  /// Internal representation of isospin weights once calculated
  mutable std::unordered_map<std::pair<uint64_t, uint64_t>, double, pair_hash>
      ratios_;

 public:
  /// Create an empty K+ N isospin ratio storage.
  KplusNRatios() : ratios_({}) {}

  /**
   * Return the isospin ratio of the given K+ N reaction's cross section.
   *
   * On the first call all ratios are calculated.
   */
  double get_ratio(const ParticleType& a, const ParticleType& b,
                   const ParticleType& c, const ParticleType& d) const;
};

extern /*thread_local (see #3075)*/ KplusNRatios kplusn_ratios;

/**
 * K- p <-> Kbar0 n cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.9
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \return the parametrized cross-section
 */
double kminusp_kbar0n(double mandelstam_s);

/**
 * K- p <-> pi- Sigma+ cross section parametrization
 * Taken from UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts the rest frame total energy
 * \return the parametrized cross-section
 */
double kminusp_piminussigmaplus(double sqrts);

/**
 * K- p <-> pi+ Sigma- cross section parametrization
 * Taken from UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts the rest frame total energy
 * \return the parametrized cross-section
 */
double kminusp_piplussigmaminus(double sqrts);

/**
 * K- p <-> pi0 Sigma0 cross section parametrization
 * Fit to Landolt-Börnstein instead of UrQMD values
 *
 * \param[in] sqrts the rest frame total energy
 * \return the parametrized cross-section
 */
double kminusp_pi0sigma0(double sqrts);

/**
 * K- p <-> pi0 Lambda cross section parametrization
 * Fit to Landolt-Börnstein instead of UrQMD values
 * \todo clarify this
 *
 * \param[in] sqrts the rest frame total energy
 * \return the parametrized cross-section
 */
double kminusp_pi0lambda(double sqrts);

/**
 * K- n <-> pi- Sigma0 cross section parametrization
 * Follow from the parametrization with the same strange
 * product via isospin symmetry.
 *
 * \param[in] sqrts the rest frame total energy
 * \return the parametrized cross-section
 */
double kminusn_piminussigma0(double sqrts);

/**
 * K- n <-> pi0 Sigma- cross section parametrization
 * Follow from the parametrization with the same strange
 * product via isospin symmetry.
 *
 * \param[in] sqrts the rest frame total energy
 * \return the parametrized cross-section
 */
double kminusn_pi0sigmaminus(double sqrts);

/**
 * K- n <-> pi- Lambda cross section parametrization
 * Follow from the parametrization with the same strange
 * product via isospin symmetry.
 *
 * \param[in] sqrts the rest frame total energy
 * \return the parametrized cross-section
 */
double kminusn_piminuslambda(double sqrts);

/**
 * Lambda Lambda <-> Xi- p cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy minus threshold energy
 * \param[in] p_N momentum of outgoing nucleon in center of mass frame
 * \param[in] p_lambda momentum of incoming lambda in center of mass frame
 * \return the parametrized cross-section
 */
double lambdalambda_ximinusp(double sqrts_sqrts0, double p_N, double p_lambda);

/**
 * Lambda Lambda <-> Xi0 n cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy minus threshold energy
 * \param[in] p_N momentum of outgoing nucleon in center of mass frame
 * \param[in] p_lambda momentum of incoming lambda in center of mass frame
 * \return the parametrized cross-section
 */
double lambdalambda_xi0n(double sqrts_sqrts0, double p_N, double p_lambda);

/**
 * Lambda Sigma+ <-> Xi0 p cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy minus threshold energy
 * \return the parametrized cross-section
 */
double lambdasigmaplus_xi0p(double sqrts_sqrts0);

/**
 * Lambda Sigma- <-> Xi- n cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy minus threshold energy
 * \return the parametrized cross-section
 */
double lambdasigmaminus_ximinusn(double sqrts_sqrts0);

/**
 * Lambda Sigma0 <-> Xi- p cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy minus threshold energy
 * \return the parametrized cross-section
 */
double lambdasigma0_ximinusp(double sqrts_sqrts0);

/**
 * Lambda Sigma0 <-> Xi0 n cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy minus threshold energy
 * \return the parametrized cross-section
 */
double lambdasigma0_xi0n(double sqrts_sqrts0);

/**
 * Sigma0 Sigma0 <-> Xi- p cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy minus threshold energy
 * \return the parametrized cross-section
 */
double sigma0sigma0_ximinusp(double sqrts_sqrts0);

/**
 * Sigma0 Sigma0 <-> Xi0 n cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * Note that there is a typo in the paper in equation (6):
 * "Lambda Sigma0 -> Xi0 n" should be "Sigma0 Sigma0 -> Xi0 n".
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy minus threshold energy
 * \return the parametrized cross-section
 */
double sigma0sigma0_xi0n(double sqrts_sqrts0);

/**
 * Sigma+ Sigma- <-> Xi0 p cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy minus threshold energy
 * \return the parametrized cross-section
 */
double sigmaplussigmaminus_xi0p(double sqrts_sqrts0);

/**
 * Sigma0 Sigma- <-> Xi- n cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy minus threshold energy
 * \return the parametrized cross-section
 */
double sigma0sigmaminus_ximinusn(double sqrts_sqrts0);

/**
 * Sigma+ Sigma- <-> Xi- p cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy minus threshold energy
 * \return the parametrized cross-section
 */
double sigmaplussigmaminus_ximinusp(double sqrts_sqrts0);

/**
 * Sigma+ Sigma- <-> Xi0 n cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy minus threshold energy
 * \return the parametrized cross-section
 */
double sigmaplussigmaminus_xi0n(double sqrts_sqrts0);

}  // namespace smash

#endif  // SRC_INCLUDE_PARAMETRIZATIONS_H_
