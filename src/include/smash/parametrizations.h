/*
 *
 *    Copyright (c) 2013-2018
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
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \param[in] is_opposite_charge whether the particles being collided have
 *                               opposite charges
 * \param[in] ma mass of first particle [GeV]
 * \param[in] mb mass of second particle [GeV]
 * \param[in] P Pomeranchuk's constant term [mb]
 * \param[in] R1 intensity of the first Regge pole contribution [mb]
 * \param[in] R2 intensity of the second Regge pole contribution [mb]
 * \return the parametrized cross-section [mb]
 */
double xs_high_energy(double mandelstam_s, bool is_opposite_charge, double ma,
                      double mb, double P, double R1, double R2);

/**
 * pp total cross section at high energies
 * \see xs_high_energy
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double pp_high_energy(double mandelstam_s);

/**
 * ppbar total cross section at high energies
 * \see xs_high_energy
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double ppbar_high_energy(double mandelstam_s);

/**
 * np total cross section at high energies
 * \see xs_high_energy
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double np_high_energy(double mandelstam_s);

/**
 * npbar total cross section at high energies
 * \see xs_high_energy
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double npbar_high_energy(double mandelstam_s);

/**
 * pi+p total cross section at high energies
 * \see xs_high_energy
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double piplusp_high_energy(double mandelstam_s);

/**
 * pi-p total cross section at high energies
 * \see xs_high_energy
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double piminusp_high_energy(double mandelstam_s);

/**
 * parametrized cross-section for proton-antiproton annihilation
 * used in the UrQMD model
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double xs_ppbar_annihilation(double mandelstam_s);

/**
 * Utility function called by specific other parametrizations
 * Parametrized hard scattering cross section (with partonic scattering)
 * This parametrization is a direct fit to cross sections in PYTHIA
 * See \iref{Sjostrand:1987su}
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \param[in] xs_0 a fit parameter [mb]
 * \param[in] e_0 a fit parameter [GeV]
 * \param[in] lambda_pow a fit parameter
 * \return the parametrized cross-section [mb]
 */
double xs_string_hard(double mandelstam_s, double xs_0, double e_0,
                      double lambda_pow);

/**
 * nucleon-nucleon hard scattering cross section (with partonic scattering)
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 *
 * \see xs_string_hard
 */
double NN_string_hard(double mandelstam_s);

/**
 * nucleon-pion hard scattering cross section (with partonic scattering)
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 *
 * \see xs_string_hard
 */
double Npi_string_hard(double mandelstam_s);

/**
 * pion-pion hard scattering cross section (with partonic scattering)
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 *
 * \see xs_string_hard
 */
double pipi_string_hard(double mandelstam_s);

/**
 * pi+p elactic cross section parametrization.
 * Source: GiBUU:parametrizationBarMes_HighEnergy.f90
 * Elastic contributions from decays are not subtracted, high energy
 * parametrization used at all energies (useful for AQM)
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \param[in] m1 the mass of the first particle [GeV]
 * \param[in] m2 the mass of the second particle [GeV]
 * \return the parametrized cross-section [mb]
 */
double piplusp_elastic_high_energy(double mandelstam_s, double m1, double m2);

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
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double piplusp_elastic(double mandelstam_s);

/**
 *  pi+ p to Sigma+ K+ cross section parametrization, PDG data.
 *
 * The PDG data is smoothed using the LOWESS algorithm. If more than one
 * cross section was given for one p_lab value, the corresponding cross sections
 * are averaged.
 */
double piplusp_sigmapluskplus_pdg(double mandelstam_s);

/**
 * pi-p elastic cross section parametrization
 * Source: GiBUU:parametrizationBarMes_HighEnergy.f90
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double piminusp_elastic(double mandelstam_s);

/**
 * pi- p -> Lambda K0 cross section parametrization, PDG data.
 *
 * The PDG data is smoothed using the LOWESS algorithm. If more than one
 * cross section was given for one p_lab value, the corresponding cross sections
 * are averaged.
 */
double piminusp_lambdak0_pdg(double mandelstam_s);

/**
 * pi- p -> Sigma- K+ cross section parametrization, PDG data.
 *
 * The PDG data is smoothed using the LOWESS algorithm. If more than one
 * cross section was given for one p_lab value, the corresponding cross sections
 * are averaged.
 */
double piminusp_sigmaminuskplus_pdg(double mandelstam_s);

/**
 * pi- p -> Sigma0 K0 cross section parametrization, resonance contribution.
 *
 * The data is smoothed using the LOWESS algorithm. If more than one
 * cross section was given for one sqrts value, the corresponding cross sections
 * are averaged.
 */
double piminusp_sigma0k0_res(double mandelstam_s);

/**
 * pp elastic cross section parametrization
 * Source: \iref{Weil:2013mya}, eq. (44)
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double pp_elastic(double mandelstam_s);

/**
 * pp elastic cross section parametrization, with only the high
 * energy part generalized to all energy regimes (used for AQM)
 * Source: \iref{Weil:2013mya}, eq. (44)
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \param[in] m1 the mass of the first particle [GeV]
 * \param[in] m2 the mass of the second particle [GeV]
 * \return the parametrized cross-section [mb]
 */
double pp_elastic_high_energy(double mandelstam_s, double m1, double m2);

/**
 * pp total cross section parametrization
 * Sources:
 * low-p: \iref{Cugnon:1996kh}
 * highest-p: \iref{Buss:2011mx}
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double pp_total(double mandelstam_s);

/**
 * np elastic cross section parametrization
 * Source: \iref{Weil:2013mya}, eq. (45)
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double np_elastic(double mandelstam_s);

/**
 * np total cross section parametrization
 * Sources:
 * low-p: \iref{Cugnon:1996kh}
 * highest-p: \iref{Buss:2011mx}
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double np_total(double mandelstam_s);

/**
 * ppbar elastic cross section parametrization
 * Source: \iref{Bass:1998ca}
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double ppbar_elastic(double mandelstam_s);

/**
 * ppbar total cross section parametrization
 * Source: \iref{Bass:1998ca}
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
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
 * Deuteron pion elastic cross-section [mb] parametrized
 * to fit pi-d elastic scattering data (the data collection
 * was be obtained from SAID data base, gwdac.phys.gwu.edu)
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double deuteron_pion_elastic(double mandelstam_s);

/**
 * Deuteron nucleon elastic cross-section [mb] parametrized
 * by \iref{Oh:2009gx}.
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double deuteron_nucleon_elastic(double mandelstam_s);

/**
 * K+ p elastic background cross section parametrization.
 * sigma(K+n->K+n) = sigma(K+n->K0p) = 0.5 * sigma(K+p->K+p)
 * Source: \iref{Buss:2011mx}, B.3.8
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double kplusp_elastic_background(double mandelstam_s);

/**
 * K+ n elastic background cross section parametrization
 * sigma(K+n->K+n) = sigma(K+n->K0p) = 0.5 * sigma(K+p->K+p)
 * Source: \iref{Buss:2011mx}, B.3.8
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double kplusn_elastic_background(double mandelstam_s);

/**
 * K+ n charge exchange cross section parametrization.
 * sigma(K+n->K+n) = sigma(K+n->K0p) = 0.5 * sigma(K+p->K+p)
 * Source: \iref{Buss:2011mx}, B.3.8
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double kplusn_k0p(double mandelstam_s);

/**
 * K- p elastic background cross section parametrization
 * Source: \iref{Buss:2011mx}, B.3.9
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double kminusp_elastic_background(double mandelstam_s);

/**
 * K- n elastic background cross section parametrization
 * Source: \iref{Buss:2011mx}, B.3.9
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double kminusn_elastic_background(double mandelstam_s);

/**
 * K0 p elastic background cross section parametrization
 * Source: \iref{Buss:2011mx}, B.3.9
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double k0p_elastic_background(double mandelstam_s);

/**
 * K0 n elastic background cross section parametrization
 * Source: \iref{Buss:2011mx}, B.3.9
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double k0n_elastic_background(double mandelstam_s);

/**
 * Kbar0 p elastic background cross section parametrization
 * Source: \iref{Buss:2011mx}, B.3.9
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double kbar0p_elastic_background(double mandelstam_s);

/**
 * Kbar0 n elastic background cross section parametrization
 * Source: \iref{Buss:2011mx}, B.3.9
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double kbar0n_elastic_background(double mandelstam_s);

/**
 * K+ p inelastic background cross section parametrization
 * Source: \iref{Buss:2011mx}, B.3.8
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double kplusp_inelastic_background(double mandelstam_s);

/**
 * K+ n inelastic background cross section parametrization
 * Source: \iref{Buss:2011mx}, B.3.8
 *
 * This interpolates the experimental data of the total cross section and
 * subtracts the elastic and charge exchange cross section.
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double kplusn_inelastic_background(double mandelstam_s);

/**
 * Hash a pair of integers.
 *
 * Note that symmetric pairs and permutations yield identical hashes with this
 * implementation.
 */
struct pair_hash {
  /// Hashing is done by this operator
  std::size_t operator()(const std::pair<uint64_t, uint64_t>& p) const {
    auto h1 = std::hash<uint64_t>{}(p.first);
    auto h2 = std::hash<uint64_t>{}(p.second);

    /* In our case the integers are PDG codes. We know they are different
     * and their order is defined, so we can simply combine the hashes
     * using XOR. Note that this yields 0 for h1 == h2. Also,
     * std::swap(h1, h2) does not not change the final hash. */
    assert(h1 != h2);
    return h1 ^ h2;
  }
};

/**
 * Calculate and store isospin ratios for K N -> K Delta reactions.
 *
 * The ratios are given by the squared Clebsch-Gordan coefficient for the
 * respective reaction, divided by the sum of the squared coefficients of all
 * possible isospin-symmetric reactions. They are used when calculating the
 * corresponding cross sections from the parametrizations of experimental data.
 */
class KaonNucleonRatios {
 private:
  /// Internal representation of isospin weights once calculated
  mutable std::unordered_map<std::pair<uint64_t, uint64_t>, double, pair_hash>
      ratios_;

 public:
  /// Create an empty K N -> K Delta isospin ratio storage.
  KaonNucleonRatios() : ratios_({}) {}

  /**
   * Return the isospin ratio of the given K N -> K Delta cross section.
   *
   * On the first call all ratios are calculated.
   */
  double get_ratio(const ParticleType& a, const ParticleType& b,
                   const ParticleType& c, const ParticleType& d) const;
};

extern /*thread_local (see #3075)*/ KaonNucleonRatios kaon_nucleon_ratios;

/**
 * K- p <-> Kbar0 n cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.9
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 */
double kminusp_kbar0n(double mandelstam_s);

/**
 * K- p <-> pi- Sigma+ cross section parametrization
 * Taken from UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double kminusp_piminussigmaplus(double sqrts);

/**
 * K- p <-> pi+ Sigma- cross section parametrization
 * Taken from UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double kminusp_piplussigmaminus(double sqrts);

/**
 * K- p <-> pi0 Sigma0 cross section parametrization
 * Fit to Landolt-Börnstein instead of UrQMD values
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double kminusp_pi0sigma0(double sqrts);

/**
 * K- p <-> pi0 Lambda cross section parametrization
 * Fit to Landolt-Börnstein instead of UrQMD values
 * \todo clarify this
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double kminusp_pi0lambda(double sqrts);

/**
 * K- n <-> pi- Sigma0 cross section parametrization
 * Follow from the parametrization with the same strange
 * product via isospin symmetry.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double kminusn_piminussigma0(double sqrts);

/**
 * K- n <-> pi0 Sigma- cross section parametrization
 * Follow from the parametrization with the same strange
 * product via isospin symmetry.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double kminusn_pi0sigmaminus(double sqrts);

/**
 * K- n <-> pi- Lambda cross section parametrization
 * Follow from the parametrization with the same strange
 * product via isospin symmetry.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double kminusn_piminuslambda(double sqrts);

/**
 * Lambda Lambda <-> Xi- p cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy
 *            minus threshold energy [GeV]
 * \param[in] p_N momentum of outgoing nucleon in center of mass frame [GeV]
 * \param[in] p_lambda momentum of incoming lambda in center of mass frame [GeV]
 * \return the parametrized cross-section [mb]
 */
double lambdalambda_ximinusp(double sqrts_sqrts0, double p_N, double p_lambda);

/**
 * Lambda Lambda <-> Xi0 n cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy
 *            minus threshold energy [GeV]
 * \param[in] p_N momentum of outgoing nucleon in center of mass frame [GeV]
 * \param[in] p_lambda momentum of incoming lambda in center of mass frame [GeV]
 * \return the parametrized cross-section [mb]
 */
double lambdalambda_xi0n(double sqrts_sqrts0, double p_N, double p_lambda);

/**
 * Lambda Sigma+ <-> Xi0 p cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy
 *            minus threshold energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double lambdasigmaplus_xi0p(double sqrts_sqrts0);

/**
 * Lambda Sigma- <-> Xi- n cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy
 *            minus threshold energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double lambdasigmaminus_ximinusn(double sqrts_sqrts0);

/**
 * Lambda Sigma0 <-> Xi- p cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy
 *            minus threshold energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double lambdasigma0_ximinusp(double sqrts_sqrts0);

/**
 * Lambda Sigma0 <-> Xi0 n cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy
 *            minus threshold energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double lambdasigma0_xi0n(double sqrts_sqrts0);

/**
 * Sigma0 Sigma0 <-> Xi- p cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy
 *            minus threshold energy [GeV]
 * \return the parametrized cross-section [mb]
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
 * \param[in] sqrts_sqrts0 the rest frame total energy
 *            minus threshold energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double sigma0sigma0_xi0n(double sqrts_sqrts0);

/**
 * Sigma+ Sigma- <-> Xi0 p cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy
 *            minus threshold energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double sigmaplussigmaminus_xi0p(double sqrts_sqrts0);

/**
 * Sigma0 Sigma- <-> Xi- n cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy
 *            minus threshold energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double sigma0sigmaminus_ximinusn(double sqrts_sqrts0);

/**
 * Sigma+ Sigma- <-> Xi- p cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy
 *            minus threshold energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double sigmaplussigmaminus_ximinusp(double sqrts_sqrts0);

/**
 * Sigma+ Sigma- <-> Xi0 n cross section parametrization
 * Two hyperon exchange, based on effective model by Feng Li,
 * as in UrQMD (\iref{Graef:2014mra}).
 *
 * \param[in] sqrts_sqrts0 the rest frame total energy
 *            minus threshold energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double sigmaplussigmaminus_xi0n(double sqrts_sqrts0);

}  // namespace smash

#endif  // SRC_INCLUDE_PARAMETRIZATIONS_H_
