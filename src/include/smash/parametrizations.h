/*
 *
 *    Copyright (c) 2013-2018,2020,2023-2024,2026
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#ifndef SRC_INCLUDE_SMASH_PARAMETRIZATIONS_H_
#define SRC_INCLUDE_SMASH_PARAMETRIZATIONS_H_

#include <optional>
#include <unordered_map>
#include <utility>

#include "particletype.h"

/* All quantities in this file use the same units as the rest of SMASH.
 * That is: GeV for energies and momenta, fm for distances and time, and mb for
 * cross-sections. */

namespace smash {

/**
 * Checks if supplied codes have existing parametrizations of total cross
 * sections
 *
 * \param[in] pdg_a PDG code of first incoming particle
 * \param[in] pdg_b PDG code of second incoming particle
 * \return Whether the parametrization exists
 */
bool parametrization_exists(const PdgCode& pdg_a, const PdgCode& pdg_b);

/**
 * total hadronic cross sections at high energies parametrized in the 2016 PDG
 * book (http://pdg.lbl.gov/2016/reviews/rpp2016-rev-cross-section-plots.pdf)
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
 *  pi+ pi- total cross section parametrized from PDG2018, smoothed using the
 * LOWESS algorithm. If the requested energy is out of the interpolation bounds,
 * the hard string value is returned.
 *
 *  \param[in] sqrts the rest frame total energy [GeV]
 *  \return the parametrized cross-section [mb]
 */
double pipluspiminus_total(double sqrts);

/**
 *  pi0 pi0 total cross section parametrized from PDG2018, smoothed using the
 * LOWESS algorithm. If the requested energy is out of the interpolation bounds,
 * the hard string value is returned.
 *
 *  \param[in] sqrts the rest frame total energy [GeV]
 *  \return the parametrized cross-section [mb]
 */
double pizeropizero_total(double sqrts);

/**
 *  pi+ p total cross section parametrized from PDG2018, smoothed using the
 * LOWESS algorithm. If the requested energy is out of the interpolation bounds,
 * the high energy cross section is returned.
 *
 *  \param[in] sqrts the rest frame total energy [GeV]
 *  \return the parametrized cross-section [mb]
 */
double piplusp_total(double sqrts);

/**
 * pi+p elactic cross section parametrization.
 *
 * Source: GiBUU:parametrizationBarMes_HighEnergy.f90
 * Elastic contributions from decays are not subtracted, high energy
 * parametrization used at all energies (useful for AQM).
 *
 * The very low part is replaced by a flat 7.5 mb cross section; used for
 * meson-meson interactions.
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \param[in] m1 the mass of the first particle [GeV]
 * \param[in] m2 the mass of the second particle [GeV]
 * \return the parametrized cross-section [mb]
 */
double piplusp_elastic_AQM(double mandelstam_s, double m1, double m2);

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
 *  pi- p total cross section parametrized from PDG2018, smoothed using the
 * LOWESS algorithm. If the requested energy is out of the interpolation bounds,
 * the high energy cross section is returned.
 *
 *  \param[in] sqrts the rest frame total energy [GeV]
 *  \return the parametrized cross-section [mb]
 */
double piminusp_total(double sqrts);

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
 * Parametrization of deuteron-pion inelastic cross section
 *
 * \param[in] pion_kinetic_energy pion kinetic energy [GeV]
 *             in the deuteron rest frame
 * \return cross section [mb]
 */
double deuteron_pion_inelastic(double pion_kinetic_energy);

/**
 * Parametrization of deuteron-nucleon inelastic cross section
 *
 * \param[in] N_kinetic_energy Nucleon kinetic energy [GeV]
 *            in the deuteron rest frame
 * \return cross section [mb]
 */
double deuteron_nucleon_inelastic(double N_kinetic_energy);

/**
 * Parametrization of deuteron-antinucleon inelastic cross section
 *
 * \param[in] aN_kinetic_energy [GeV] Anti-nucleon kinetic energy
 *             in the deuteron rest frame
 * \return cross section [mb]
 */
double deuteron_antinucleon_inelastic(double aN_kinetic_energy);

/**
 * K+ p total cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.8
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 *
 * \note \anchor KN_note In total parametrizations of KN processes,
 * if the interaction energy exceeds the bounds of the interpolation,
 * the last value available is returned, which is desired behavior.
 */
double kplusp_total(double mandelstam_s);

/**
 * K+ n total cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.8
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 *
 * \note See \ref KN_note "this note" about the return value.
 */
double kplusn_total(double mandelstam_s);

/**
 * K- n total cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.8
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 *
 * \note See \ref KN_note "this note" about return value.
 */
double kminusn_total(double mandelstam_s);

/**
 * K- p total cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.8
 *
 * \param[in] mandelstam_s the rest frame total energy squared [GeV^2]
 * \return the parametrized cross-section [mb]
 *
 * \note See \ref KN_note "this note" about return value.
 */
double kminusp_total(double mandelstam_s);

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

extern /*thread_local (see commit 897d0b8)*/ KaonNucleonRatios
    kaon_nucleon_ratios;

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
 * K- n <-> pi0 Sigma- cross section parametrization is
 * also handled with this.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double kminusn_piminussigma0(double sqrts);

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

/**
 * D⁰π⁺ elastic cross section (\iref{Abreu:2011ic}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D̄⁰π⁻ is also handled with
 * this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dzeropiplus_elastic(double sqrts);

/**
 * D⁰π⁺ -> D⁺π⁰ cross section (\iref{Abreu:2011ic}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D̄⁰π⁻ -> D⁻π⁰ is also handled
 * with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double Dzeropiplus_Dpluspizero(double sqrts);

/**
 * D⁰π⁻ elastic cross section (\iref{Abreu:2011ic}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D̄⁰π⁺ is also handled with
 * this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dzeropiminus_elastic(double sqrts);

/**
 * D⁰π⁰ elastic cross section (\iref{Abreu:2011ic}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D̄⁰π⁰ is also handled with
 * this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dzeropizero_elastic(double sqrts);

/**
 * D⁰π⁰ -> D⁺π⁻ cross section (\iref{Abreu:2011ic}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D̄⁰π⁰ -> D⁻π⁺ is also handled
 * with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double Dzeropizero_Dpluspiminus(double sqrts);

/**
 * D⁺π⁺ elastic cross section (\iref{Abreu:2011ic}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D⁻π⁻ is also handled with
 * this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dpluspiplus_elastic(double sqrts);

/**
 * D⁺π⁻ elastic cross section (\iref{Abreu:2011ic}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D⁻π⁺ is also handled with
 * this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dpluspiminus_elastic(double sqrts);

/**
 * D⁺π⁻ -> D⁰π⁰ cross section (\iref{Abreu:2011ic}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D⁻π⁺ -> D̄⁰π⁰ is also handled
 * with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double Dpluspiminus_Dzeropizero(double sqrts);

/**
 * D⁺π⁰ elastic cross section (\iref{Abreu:2011ic}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D⁻π⁰ is also handled with
 * this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dpluspizero_elastic(double sqrts);

/**
 * D⁺π⁰ -> D⁰π⁺ cross section (\iref{Abreu:2011ic}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D⁻π⁰ -> D̄⁰π⁻ is also handled
 * with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double Dpluspizero_Dzeropiplus(double sqrts);

/**
 * D⁺η elastic cross section (\iref{Tolos:2013kva}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D⁻η is also handled with this
 * function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dpluseta_elastic(double sqrts);

/**
 * D⁰η elastic cross section (\iref{Tolos:2013kva}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D̄⁰η is also handled with this
 * function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dzeroeta_elastic(double sqrts);

/**
 * D⁺K⁺ elastic cross section (\iref{Tolos:2013kva}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D⁻K⁻ is also handled with
 * this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> DplusKplus_elastic(double sqrts);

/**
 * D⁺K⁰ elastic cross section (\iref{Tolos:2013kva}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D⁻K̄⁰ is also handled with
 * this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> DplusKzero_elastic(double sqrts);

/**
 * D⁺K⁰ -> D⁰K⁺ cross section (\iref{Tolos:2013kva}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D⁻K̄⁰ -> D̄⁰K⁻ is also handled
 * with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double DplusKzero_DzeroKplus(double sqrts);

/**
 * D⁰K⁺ -> D⁺K⁰ cross section (\iref{Tolos:2013kva}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D̄⁰K⁻ -> D⁻K̄⁰ is also handled
 * with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double DzeroKplus_DplusKzero(double sqrts);

/**
 * D⁰K⁺ elastic cross section (\iref{Tolos:2013kva}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D̄⁰K⁻ is also handled with
 * this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> DzeroKplus_elastic(double sqrts);

/**
 * D⁰K⁰ elastic cross section (\iref{Tolos:2013kva}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D̄⁰K̄⁰ is also handled with
 * this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> DzeroKzero_elastic(double sqrts);

/**
 * D⁺K̄⁰ elastic cross section (\iref{Tolos:2013kva}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D⁻K⁰ is also handled with
 * this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> DplusKbarzero_elastic(double sqrts);

/**
 * D⁺K⁻ elastic cross section (\iref{Tolos:2013kva}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D⁻K⁺ is also handled with
 * this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> DplusKminus_elastic(double sqrts);

/**
 * D⁺K⁻ -> D⁰K̄⁰ cross section (\iref{Tolos:2013kva}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D⁻K⁺ -> D̄⁰K⁰ is also handled
 * with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double DplusKminus_DzeroKbarzero(double sqrts);

/**
 * D⁰K̄⁰ -> D⁺K⁻ cross section (\iref{Tolos:2013kva}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D̄⁰K⁰ -> D⁻K⁺ is also handled
 * with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double DzeroKbarzero_DplusKminus(double sqrts);

/**
 * D⁰K̄⁰ elastic cross section (\iref{Tolos:2013kva}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D̄⁰K⁰ is also handled with
 * this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> DzeroKbarzero_elastic(double sqrts);

/**
 * D⁰K⁻ elastic cross section (\iref{Tolos:2013kva}, data provided by Juan
 * Torres-Rincon). Charge conjugated cross section D̄⁰K⁺ is also handled with
 * this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> DzeroKminus_elastic(double sqrts);

/**
 * D*(2010)⁺π⁺ elastic cross section (closest reference \iref{Song:2015sfa},
 * data provided by Juan Torres-Rincon). Charge conjugated cross section
 * D*(2010)⁻π⁻ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dstarpluspiplus_elastic(double sqrts);

/**
 * D*(2010)⁺π⁻ elastic cross section (closest reference \iref{Song:2015sfa},
 * data provided by Juan Torres-Rincon). Charge conjugated cross section
 * D*(2010)⁻π⁺ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dstarpluspiminus_elastic(double sqrts);

/**
 * D*(2010)⁺π⁻ -> D*(2007)⁰π⁰ cross section (closest reference
 * \iref{Song:2015sfa}, data provided by Juan Torres-Rincon). Charge conjugated
 * cross section D*(2010)⁻π⁺ -> D̄*(2007)⁰π⁰ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double Dstarpluspiminus_Dstarzeropizero(double sqrts);

/**
 * D*(2010)⁺π⁰ elastic cross section (closest reference \iref{Song:2015sfa},
 * data provided by Juan Torres-Rincon). Charge conjugated cross section
 * D*(2010)⁻π⁰ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dstarpluspizero_elastic(double sqrts);

/**
 * D*(2010)⁺π⁰ -> D*(2007)⁰π⁺ cross section (closest reference
 * \iref{Song:2015sfa}, data provided by Juan Torres-Rincon). Charge conjugated
 * cross section D*(2010)⁻π⁰ -> D̄*(2007)⁰π⁻ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double Dstarpluspizero_Dstarzeropiplus(double sqrts);

/**
 * D*(2007)⁰π⁺ -> D*(2010)⁺π⁰ cross section (closest reference
 * \iref{Song:2015sfa}, data provided by Juan Torres-Rincon). Charge conjugated
 * cross section D̄*(2007)⁰π⁻ -> D*(2010)⁻π⁰ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double Dstarzeropiplus_Dstarpluspizero(double sqrts);

/**
 * D*(2007)⁰π⁺ elastic cross section (closest reference \iref{Song:2015sfa},
 * data provided by Juan Torres-Rincon). Charge conjugated cross section
 * D̄*(2007)⁰π⁻ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dstarzeropiplus_elastic(double sqrts);

/**
 * D*(2007)⁰π- elastic cross section (closest reference \iref{Song:2015sfa},
 * data provided by Juan Torres-Rincon). Charge conjugated cross section
 * D̄*(2007)⁰π⁺ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dstarzeropiminus_elastic(double sqrts);

/**
 * D*(2007)⁰π⁰ -> D*(2010)⁺π⁻ cross section (closest reference
 * \iref{Song:2015sfa}, data provided by Juan Torres-Rincon). Charge conjugated
 * cross section D̄*(2007)⁰π⁰ -> D*(2010)⁻π⁺ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double Dstarzeropizero_Dstarpluspiminus(double sqrts);

/**
 * D*(2007)⁰π⁰ elastic cross section (closest reference \iref{Song:2015sfa},
 * data provided by Juan Torres-Rincon). Charge conjugated cross section
 * D̄*(2007)⁰π⁰ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dstarzeropizero_elastic(double sqrts);

/**
 * D*(2010)⁺η elastic cross section (data provided by Juan Torres-Rincon).
 * Charge conjugated cross section D*(2010)⁻η is also handled with this
 * function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dstarpluseta_elastic(double sqrts);

/**
 * D*(2007)⁰η elastic cross section (data provided by Juan Torres-Rincon).
 * Charge conjugated cross section D̄*(2007)⁰η is also handled with this
 * function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dstarzeroeta_elastic(double sqrts);

/**
 * D*(2010)⁺K⁺ elastic cross section (closest reference \iref{Song:2015sfa},
 * data provided by Juan Torres-Rincon). Charge conjugated cross section
 * D*(2010)⁻K⁻ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> DstarplusKplus_elastic(double sqrts);

/**
 * D*(2010)⁺K⁰ elastic cross section (closest reference \iref{Song:2015sfa},
 * data provided by Juan Torres-Rincon). Charge conjugated cross section
 * D*(2010)⁻K̄⁰ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> DstarplusKzero_elastic(double sqrts);

/**
 * D*(2010)⁺K⁰ -> D*(2007)⁰K⁺ cross section (closest reference
 * \iref{Song:2015sfa}, data provided by Juan Torres-Rincon). Charge conjugated
 * cross section D*(2010)⁻K̄⁰ -> D̄*(2007)⁰K⁻ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double DstarplusKzero_DstarzeroKplus(double sqrts);

/**
 * D*(2007)⁰K⁺ -> D*(2010)⁺K⁰ cross section (closest reference
 * \iref{Song:2015sfa}, data provided by Juan Torres-Rincon). Charge conjugated
 * cross section D̄*(2007)⁰K⁻ -> D*(2010)⁻K̄⁰ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double DstarzeroKplus_DstarplusKzero(double sqrts);

/**
 * D*(2007)⁰K⁺ elastic cross section (closest reference \iref{Song:2015sfa},
 * data provided by Juan Torres-Rincon). Charge conjugated cross section
 * D̄*(2007)⁰K⁻ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> DstarzeroKplus_elastic(double sqrts);

/**
 * D*(2007)⁰K⁰ elastic cross section (closest reference \iref{Song:2015sfa},
 * data provided by Juan Torres-Rincon). Charge conjugated cross section
 * D̄*(2007)⁰K̄⁰ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> DstarzeroKzero_elastic(double sqrts);

/**
 * D*(2010)⁺K̄⁰ elastic cross section (closest reference \iref{Song:2015sfa},
 * data provided by Juan Torres-Rincon). Charge conjugated cross section
 * D*(2010)⁻K⁰ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> DstarplusKbarzero_elastic(double sqrts);

/**
 * D*(2010)⁺K⁻ elastic cross section (closest reference \iref{Song:2015sfa},
 * data provided by Juan Torres-Rincon). Charge conjugated cross section
 * D*(2010)⁻K⁺ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> DstarplusKminus_elastic(double sqrts);

/**
 * D*(2010)⁺K⁻ -> D*(2007)⁰K̄⁰ cross section (\iref{Tolos:2013kva}, data provided
 * by Juan Torres-Rincon). Charge conjugated cross section
 * D*(2010)⁻K⁺ -> D̄*(2007)⁰K⁰ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double DstarplusKminus_DstarzeroKbarzero(double sqrts);

/**
 * D*(2007)⁰K̄⁰ -> D*(2010)⁺K⁻ cross section (\iref{Tolos:2013kva}, data provided
 * by Juan Torres-Rincon). Charge conjugated cross section
 * D̄*(2007)⁰K⁰ -> D*(2010)⁻K⁺ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double DstarzeroKbarzero_DstarplusKminus(double sqrts);

/**
 * D*(2007)⁰K̄⁰ elastic cross section (closest reference \iref{Song:2015sfa},
 * data provided by Juan Torres-Rincon). Charge conjugated cross section
 * D̄*(2007)⁰K⁰ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> DstarzeroKbarzero_elastic(double sqrts);

/**
 * D*(2007)⁰K⁻ elastic cross section (closest reference \iref{Song:2015sfa},
 * data provided by Juan Torres-Rincon). Charge conjugated cross section
 * D̄*(2007)⁰K⁺ is also handled with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> DstarzeroKminus_elastic(double sqrts);

/**
 * D⁺n elastic cross section (\iref{Tolos:2013kva}), data provided by Juan
 * Torres-Rincon. Charge conjugated cross section D⁻n̄ is also handled with this
 * function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dplusn_elastic(double sqrts);

/**
 * D⁺n -> D⁰p cross section (\iref{Tolos:2013kva}), data provided by Juan
 * Torres-Rincon. Charge conjugated cross section D⁻n̄ -> D̄⁰p̄ is also handled
 * with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double Dplusn_Dzerop(double sqrts);

/**
 * D⁺p elastic cross section (\iref{Tolos:2013kva}), data provided by Juan
 * Torres-Rincon. Charge conjugated cross section D⁻p̄ is also handled with this
 * function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dplusp_elastic(double sqrts);

/**
 * D⁰n elastic cross section (\iref{Tolos:2013kva}), data provided by Juan
 * Torres-Rincon. Charge conjugated cross section D̄⁰n̄ is also handled with this
 * function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dzeron_elastic(double sqrts);

/**
 * D⁰p -> D⁺n cross section (\iref{Tolos:2013kva}), data provided by Juan
 * Torres-Rincon. Charge conjugated cross section D̄⁰p̄ -> D⁻n̄ is also handled
 * with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double Dzerop_Dplusn(double sqrts);

/**
 * D⁰p elastic cross section (\iref{Tolos:2013kva}), data provided by Juan
 * Torres-Rincon. Charge conjugated cross section D̄⁰p̄ is also handled with this
 * function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dzerop_elastic(double sqrts);

/**
 * D⁻n elastic cross section (\iref{Tolos:2013kva}), data provided by Juan
 * Torres-Rincon. Charge conjugated cross section D⁺n̄ is also handled with this
 * function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dminusn_elastic(double sqrts);

/**
 * D⁻p elastic cross section (\iref{Tolos:2013kva}), data provided by Juan
 * Torres-Rincon. Charge conjugated cross section D⁺p̄ is also handled with this
 * function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dminusp_elastic(double sqrts);

/**
 * D⁻p -> D̄⁰n cross section (\iref{Tolos:2013kva}), data provided by Juan
 * Torres-Rincon. Charge conjugated cross section D⁺p̄ -> D⁰n̄ is also handled
 * with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double Dminusp_Dbarzeron(double sqrts);

/**
 * D̄⁰n -> D⁻p cross section (\iref{Tolos:2013kva}), data provided by Juan
 * Torres-Rincon. Charge conjugated cross section D⁰n̄ -> D⁺p̄ is also handled
 * with this function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
double Dbarzeron_Dminusp(double sqrts);

/**
 * D̄⁰n elastic cross section (\iref{Tolos:2013kva}), data provided by Juan
 * Torres-Rincon. Charge conjugated cross section D⁰n̄ is also handled with this
 * function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dbarzeron_elastic(double sqrts);

/**
 * D̄⁰p elastic cross section (\iref{Tolos:2013kva}), data provided by Juan
 * Torres-Rincon. Charge conjugated cross section D⁰p̄ is also handled with this
 * function.
 *
 * \param[in] sqrts the rest frame total energy [GeV]
 * \return the parametrized cross-section [mb]
 */
std::optional<double> Dbarzerop_elastic(double sqrts);

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_PARAMETRIZATIONS_H_
