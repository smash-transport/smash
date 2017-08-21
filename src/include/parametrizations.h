/*
 *
 *    Copyright (c) 2013-2014
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

namespace Smash {

/** total hadronic cross sections at high energies parametrized in the 2016 PDG
 *  book(http://pdg.lbl.gov/2016/reviews/rpp2016-rev-cross-section-plots.pdf) */
double xs_high_energy(double mandelstam_s, bool is_opposite_charge, double ma,
                      double mb, double P, double R1, double R2);

/** pp total cross section at high energies */
double pp_high_energy(double mandelstam_s);

/** ppbar total cross section at high energies */
double ppbar_high_energy(double mandelstam_s);

/** np total cross section at high energies */
double np_high_energy(double mandelstam_s);

/** npbar total cross section at high energies */
double npbar_high_energy(double mandelstam_s);

/** pi+p total cross section at high energies */
double piplusp_high_energy(double mandelstam_s);

/** pi-p total cross section at high energies */
double piminusp_high_energy(double mandelstam_s);

/** pi+p elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
double piplusp_elastic(double mandelstam_s);

/** pi-p elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
double piminusp_elastic(double mandelstam_s);

/** pp elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
double pp_elastic(double mandelstam_s);

/** pp total cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
double pp_total(double mandelstam_s);

/** np elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
double np_elastic(double mandelstam_s);

/** np total cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
double np_total(double mandelstam_s);

/** ppbar elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
double ppbar_elastic(double mandelstam_s);

/** ppbar total cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
double ppbar_total(double mandelstam_s);

/** ppbar total cross section parametrization;
 * Used for reverse cross-section from detailed balance
 */
double ppbar_total(double mandelstam_s, double m_proj, double m_target);

/** K+ p elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
double kplusp_elastic(double mandelstam_s);

/** K+ n elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
double kplusn_elastic(double mandelstam_s);

/** K- p elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
double kminusp_elastic(double mandelstam_s);

/** K- n elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
double kminusn_elastic(double mandelstam_s);

/** K0 p elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
double k0p_elastic(double mandelstam_s);

/** K0 n elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
double k0n_elastic(double mandelstam_s);

/** Kbar0 p elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
double kbar0p_elastic(double mandelstam_s);

/** Kbar0 n elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
double kbar0n_elastic(double mandelstam_s);

/** K+ p inelastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
double kplusp_inelastic(double mandelstam_s);

/** K+ n inelastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
double kplusn_inelastic(double mandelstam_s);

/** Hash a pair of integers.
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

/** Isospin weights for inelastic K+ N channels.
 */
class KplusNRatios {
 private:
  mutable std::unordered_map<std::pair<uint64_t, uint64_t>, double, pair_hash>
      ratios_;

 public:
  /// Create an empty K+ N isospin ratio storage.
  KplusNRatios() : ratios_({}) {}

  /// Return the isospin ratio of the given K+ N reaction's cross section.
  ///
  /// On the first call all ratios are calculated.
  double get_ratio(const ParticleType& a, const ParticleType& b,
                   const ParticleType& c, const ParticleType& d) const;
};

extern /*thread_local (see #3075)*/ KplusNRatios kplusn_ratios;

/** K- p <-> Kbar0 n cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.9 */
double kminusp_kbar0n(double mandelstam_s);

/// K- p <-> pi- Sigma+ cross section parametrization
double kminusp_piminussigmaplus(double sqrts);

/// K- p <-> pi+ Sigma- cross section parametrization
double kminusp_piplussigmaminus(double sqrts);

/// K- p <-> pi0 Sigma0 cross section parametrization
double kminusp_pi0sigma0(double sqrts);

/// K- p <-> pi0 Lambda cross section parametrization
double kminusp_pi0lambda(double sqrts);

/// K- n <-> pi- Sigma0 cross section parametrization
double kminusn_piminussigma0(double sqrts);

/// K- n <-> pi0 Sigma- cross section parametrization
double kminusn_pi0sigmaminus(double sqrts);

/// K- n <-> pi- Lambda cross section parametrization
double kminusn_piminuslambda(double sqrts);

/// Lambda Lambda <-> Xi- p cross section parametrization
double lambdalambda_ximinusp(double sqrts_sqrts0, double p_N, double p_lambda);

/// Lambda Lambda <-> Xi0 n cross section parametrization
double lambdalambda_xi0n(double sqrts_sqrts0, double p_N, double p_lambda);

/// Lambda Sigma+ <-> Xi0 p cross section parametrization
double lambdasigmaplus_xi0p(double sqrts_sqrts0);

/// Lambda Sigma- <-> Xi- n cross section parametrization
double lambdasigmaminus_ximinusn(double sqrts_sqrts0);

/// Lambda Sigma0 <-> Xi- p cross section parametrization
double lambdasigma0_ximinusp(double sqrts_sqrts0);

/// Lambda Sigma0 <-> Xi0 n cross section parametrization
double lambdasigma0_xi0n(double sqrts_sqrts0);

/// Sigma0 Sigma0 <-> Xi- p cross section parametrization
double sigma0sigma0_ximinusp(double sqrts_sqrts0);

/// Sigma0 Sigma0 <-> Xi0 n cross section parametrization
double sigma0sigma0_xi0n(double sqrts_sqrts0);

/// Sigma+ Sigma- <-> Xi0 p cross section parametrization
double sigmaplussigmaminus_xi0p(double sqrts_sqrts0);

/// Sigma0 Sigma- <-> Xi- n cross section parametrization
double sigma0sigmaminus_ximinusn(double sqrts_sqrts0);

/// Sigma+ Sigma- <-> Xi- p cross section parametrization
double sigmaplussigmaminus_ximinusp(double sqrts_sqrts0);

/// Sigma+ Sigma- <-> Xi0 n cross section parametrization
double sigmaplussigmaminus_xi0n(double sqrts_sqrts0);

}  // namespace Smash

#endif  // SRC_INCLUDE_PARAMETRIZATIONS_H_
