/*
 *    Copyright (c) 2015-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_TABULATION_H_
#define SRC_INCLUDE_TABULATION_H_

#include <functional>
#include <map>
#include <memory>
#include <vector>

#include "forwarddeclarations.h"
#include "integrate.h"
#include "kinematics.h"
#include "particletype.h"

namespace Smash {

enum class Extrapolation {
  Zero = 0,
  Const = 1,
  Linear = 2,
};

/**
 * A class for storing a one-dimensional lookup table of doubleing-point values.
 */
class Tabulation {
 public:
  /** Construct a new tabulation object.
   * \param x_min lower bound of tabulation domain
   * \param range range (x_max-x_min) of tabulation domain
   * \param num number of intervals (the number of tabulated points is actually
   * num+1)
   * \param f one-dimensional function f(x) which is supposed to be tabulated
   */
  Tabulation(double x_min, double range, int num,
             std::function<double(double)> f);
  /** Look up a value from the tabulation (without any interpolation, simply
   * using the closest tabulated value). If x is below the lower tabulation
   * bound we return 0, if it is above the upper bound we return the tabulated
   * value at the upper bound. */
  double get_value_step(double x) const;
  /** Look up a value from the tabulation using linear interpolation.
   * If x is above the upper bound, then by default we use linear extrapolation
   * of the two highest tabulated points. Optionally one can also extrapolate
   * with rightmost value or zero. Linear extrapolation is not an arbitrary
   * choice, in fact many functions tabulated in SMASH have a linear
   * asymptotics, e.g. rho(m) functions. */
  double get_value_linear(
      double x, Extrapolation extrapolation = Extrapolation::Linear) const;

 protected:
  // vector for storing tabulated values
  std::vector<double> values_;
  // lower bound and inverse step size 1/dx for tabulation
  const double x_min_, x_max_, inv_dx_;
};

/**
 * Spectral function integrand for GSL integration, with one resonance in the
 * final state (the second particle is stable).
 *
 * The integrand is \f$ A(m) p_{cm}^f \f$, where \f$ m \f$ is the
 * resonance mass, \f$ A(m) \f$ is the spectral function
 *  and \f$ p_{cm}^f \f$ is the center-of-mass momentum of the final state.
 *
 * \param[in] resonance_mass Actual mass of the resonance.
 * \param[in] sqrts Center-of-mass energy, i.e. sqrt of Mandelstam s.
 * \param[in] stable_mass mass of the stable particle in the final state
 * \param[in] type type of the resonance
 */
inline double spec_func_integrand_1res(double resonance_mass, double sqrts,
                                       double stable_mass,
                                       const ParticleType& type) {
  if (sqrts <= stable_mass + resonance_mass) {
    return 0.;
  }

  /* Integrand is the spectral function weighted by the CM momentum of the
   * final state. */
  return type.spectral_function(resonance_mass) *
         pCM(sqrts, stable_mass, resonance_mass);
}

/**
 * Spectral function integrand for GSL integration, with two resonances in the
 * final state.
 *
 * The integrand is \f$ A_1(m_1) A_2(m_2) p_{cm}^f \f$, where \f$ m_1 \f$ and
 * \f$ m_2 \f$ are the resonance masses, \f$ A_1 \f$ and \f$ A_2 \f$ are the
 * spectral functions and \f$ p_{cm}^f \f$ is the center-of-mass momentum of
 * the final state.
 *
 * \param[in] sqrts Center-of-mass energy, i.e. sqrt of Mandelstam s.
 * \param[in] res_mass_1 Actual mass of the first resonance.
 * \param[in] res_mass_2 Actual mass of the second resonance.
 * \param[in] t1 Type of the first resonance.
 * \param[in] t2 Type of the second resonance.
 */
inline double spec_func_integrand_2res(double sqrts, double res_mass_1,
                                       double res_mass_2,
                                       const ParticleType& t1,
                                       const ParticleType& t2) {
  if (sqrts <= res_mass_1 + res_mass_2) {
    return 0.;
  }

  /* Integrand is the product of the spectral function weighted by the
   * CM momentum of the final state. */
  return t1.spectral_function(res_mass_1) * t2.spectral_function(res_mass_2) *
         pCM(sqrts, res_mass_1, res_mass_2);
}

/**
 * Create a table for the spectral integral
 *  of a resonance and a stable particle.
 */
inline std::unique_ptr<Tabulation> spectral_integral_semistable(
    Integrator& integrate, const ParticleType& resonance,
    const ParticleType& stable, double range) {
  return make_unique<Tabulation>(resonance.min_mass_kinematic() + stable.mass(),
                                 range, 100, [&](double srts) {
                                   return integrate(
                                       resonance.min_mass_kinematic(),
                                       srts - stable.mass(), [&](double m) {
                                         return spec_func_integrand_1res(
                                             m, srts, stable.mass(), resonance);
                                       });
                                 });
}

/// Create a table for the spectral integral of two resonances.
inline std::unique_ptr<Tabulation> spectral_integral_unstable(
    Integrator2dCuhre& integrate2d, const ParticleType& res1,
    const ParticleType& res2, double range) {
  return make_unique<Tabulation>(
      res1.min_mass_kinematic() + res2.min_mass_kinematic(), range, 100,
      [&](double srts) {
        return integrate2d(
            res1.min_mass_kinematic(), srts - res2.min_mass_kinematic(),
            res2.min_mass_kinematic(), srts - res1.min_mass_kinematic(),
            [&](double m1, double m2) {
              return spec_func_integrand_2res(srts, m1, m2, res1, res2);
            });
      });
}

}  // namespace Smash

#endif  // SRC_INCLUDE_TABULATION_H_
