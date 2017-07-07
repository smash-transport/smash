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
float xs_high_energy(double mandelstam_s, bool is_opposite_charge,
                     float ma, float mb, float P, float R1, float R2);

/** pp total cross section at high energies */
float pp_high_energy(double mandelstam_s);

/** ppbar total cross section at high energies */
float ppbar_high_energy(double mandelstam_s);

/** np total cross section at high energies */
float np_high_energy(double mandelstam_s);

/** npbar total cross section at high energies */
float npbar_high_energy(double mandelstam_s);

/** pi+p total cross section at high energies */
float piplusp_high_energy(double mandelstam_s);

/** pi-p total cross section at high energies */
float piminusp_high_energy(double mandelstam_s);

/** pi+p elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float piplusp_elastic(double mandelstam_s);

/** pi-p elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float piminusp_elastic(double mandelstam_s);

/** pp elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float pp_elastic(double mandelstam_s);

/** pp total cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float pp_total(double mandelstam_s);

/** np elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float np_elastic(double mandelstam_s);

/** np total cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float np_total(double mandelstam_s);


/** ppbar elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float ppbar_elastic(double mandelstam_s);

/** ppbar total cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float ppbar_total(double mandelstam_s);

/** ppbar total cross section parametrization;
 * Used for reverse cross-section from detailed balance
 */
float ppbar_total(double mandelstam_s, double m_proj, double m_target);

/** K+ p elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float kplusp_elastic(double mandelstam_s);

/** K+ n elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float kplusn_elastic(double mandelstam_s);

/** K- p elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float kminusp_elastic(double mandelstam_s);

/** K- n elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float kminusn_elastic(double mandelstam_s);

/** K0 p elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float k0p_elastic(double mandelstam_s);

/** K0 n elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float k0n_elastic(double mandelstam_s);

/** Kbar0 p elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float kbar0p_elastic(double mandelstam_s);

/** Kbar0 n elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float kbar0n_elastic(double mandelstam_s);

/** K+ p inelastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float kplusp_inelastic(double mandelstam_s);

/** K+ n inelastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float kplusn_inelastic(double mandelstam_s);

/** Hash a pair of integers.
 *
 * Note that symmetric pairs and permutations yield identical hashes with this
 * implementation.
 */
struct pair_hash {
    std::size_t operator () (const std::pair<uint64_t, uint64_t> &p) const {
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
  mutable std::unordered_map<std::pair<uint64_t, uint64_t>,
                             float, pair_hash> ratios_;
 public:
  /// Create an empty K+ N isospin ratio storage.
  KplusNRatios() : ratios_({}) {}

  /// Return the isospin ratio of the given K+ N reaction's cross section.
  ///
  /// On the first call all ratios are calculated.
  float get_ratio(const ParticleType& a, const ParticleType& b,
                  const ParticleType& c, const ParticleType& d) const;
};

extern thread_local KplusNRatios kplusn_ratios;

/** K- p <-> Kbar0 n cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.9 */
float kminusp_kbar0n(double mandelstam_s);

/// K- p <-> pi- Sigma+ cross section parametrization
float kminusp_piminussigmaplus(double sqrts);

/// K- p <-> pi+ Sigma- cross section parametrization
float kminusp_piplussigmaminus(double sqrts);

/// K- p <-> pi0 Sigma0 cross section parametrization
float kminusp_pi0sigma0(double sqrts);

/// K- p <-> pi0 Lambda cross section parametrization
float kminusp_pi0lambda(double sqrts);

/// K- n <-> pi- Sigma0 cross section parametrization
float kminusn_piminussigma0(double sqrts);

/// K- n <-> pi0 Sigma- cross section parametrization
float kminusn_pi0sigmaminus(double sqrts);

/// K- n <-> pi- Lambda cross section parametrization
float kminusn_piminuslambda(double sqrts);

/// Lambda Lambda <-> Xi- p cross section parametrization
float lambdalambda_ximinusp(double sqrts_sqrts0, double p_N, double p_lambda);

/// Lambda Lambda <-> Xi0 n cross section parametrization
float lambdalambda_xi0n(double sqrts_sqrts0, double p_N, double p_lambda);

/// Lambda Sigma+ <-> Xi0 p cross section parametrization
float lambdasigmaplus_xi0p(double sqrts_sqrts0);

/// Lambda Sigma- <-> Xi- n cross section parametrization
float lambdasigmaminus_ximinusn(double sqrts_sqrts0);

/// Lambda Sigma0 <-> Xi- p cross section parametrization
float lambdasigma0_ximinusp(double sqrts_sqrts0);

/// Lambda Sigma0 <-> Xi0 n cross section parametrization
float lambdasigma0_xi0n(double sqrts_sqrts0);

/// Sigma0 Sigma0 <-> Xi- p cross section parametrization
float sigma0sigma0_ximinusp(double sqrts_sqrts0);

/// Sigma0 Sigma0 <-> Xi0 n cross section parametrization
float sigma0sigma0_xi0n(double sqrts_sqrts0);

/// Sigma+ Sigma- <-> Xi0 p cross section parametrization
float sigmaplussigmaminus_xi0p(double sqrts_sqrts0);

/// Sigma0 Sigma- <-> Xi- n cross section parametrization
float sigma0sigmaminus_ximinusn(double sqrts_sqrts0);

/// Sigma+ Sigma- <-> Xi- p cross section parametrization
float sigmaplussigmaminus_ximinusp(double sqrts_sqrts0);

/// Sigma+ Sigma- <-> Xi0 n cross section parametrization
float sigmaplussigmaminus_xi0n(double sqrts_sqrts0);

}  // namespace Smash

#endif  // SRC_INCLUDE_PARAMETRIZATIONS_H_
