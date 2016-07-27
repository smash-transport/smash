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

namespace Smash {

/* pp elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float pp_elastic(double mandelstam_s);

/* pp total cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float pp_total(double mandelstam_s);

/* np elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float np_elastic(double mandelstam_s);

/* np total cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float np_total(double mandelstam_s);

/* ppbar elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float ppbar_elastic(double mandelstam_s);

/* ppbar total cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float ppbar_total(double mandelstam_s);

/* K+ p elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float kplusp_elastic(double mandelstam_s);

/* K+ n elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float kplusn_elastic(double mandelstam_s);

/* K- p elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float kminusp_elastic(double mandelstam_s);

/* K- n elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float kminusn_elastic(double mandelstam_s);

/* K0 p elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float k0p_elastic(double mandelstam_s);

/* K0 n elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float k0n_elastic(double mandelstam_s);

/* Kbar0 p elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float kbar0p_elastic(double mandelstam_s);

/* Kbar0 n elastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float kbar0n_elastic(double mandelstam_s);

/* K+ p inelastic cross section parametrization
 *
 * \fpPrecision Why \c double?
 */
float kplusp_inelastic(double mandelstam_s);

/// K- p <-> pi- Sigma+ cross section parametrization
float kminusp_piminussigmaplus(double sqrts);

/// K- p <-> pi+ Sigma- cross section parametrization
float kminusp_piplussigmaminus(double sqrts);

/// K- p <-> pi0 Sigma0 cross section parametrization
float kminusp_pi0sigma0(double sqrts);

/// K- p <-> pi0 Lambda cross section parametrization
float kminusp_pi0lambda(double sqrts);

/// Lambda Lambda <-> Xi- p cross section parametrization
float lambdalambda_ximinusp(double sqrts_sqrts0, double p_N, double p_lambda);

/// Lambda Lambda <-> Xi0 n cross section parametrization
float lambdalambda_xi0n(double sqrts_sqrts0, double p_N, double p_lambda);

/// Lambda Lambda <-> Xi N cross section parametrization
float lambdalambda_xiN(double sqrts_sqrts0, double p_N, double p_lambda);

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
