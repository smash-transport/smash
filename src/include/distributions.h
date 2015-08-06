/*
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_DISTRIBUTIONS_H_
#define SRC_INCLUDE_DISTRIBUTIONS_H_

namespace Smash {

/** Returns a Breit-Wigner distribution
 *
 * \param[in] mandelstam_s the Mandelstam s variable (available energy
 * squared - in GeV\f$^2\f$)
 * \param[in] resonance_mass resonance pole mass in GeV
 * \param[in] resonance_width resonance width in GeV
 *
 * \return \f$\frac{s \Gamma^2}{(s-m^2)^2 + s\Gamma^2}\f$
 *
 * \fpPrecision Why \c double?
 */
float breit_wigner(const double mandelstam_s, const float resonance_mass,
                    const float resonance_width);

/** Returns the Maxwell-Boltzmann distribution
 *
 * \todo rename this function to make clear what it is
 *
 * \param[in] energy \f$E\f$ (in GeV)
 * \param[in] momentum_sqr squared \f$p\f$ (in GeV\f$^2\f$)
 * \param[in] temperature \f$T\f$ (in GeV)
 *
 * \return \f$4\pi p^2 \exp{-\frac{E}{T}}\f$
 *
 * \fpPrecision Why \c double?
 */
double density_integrand(const double energy, const double momentum_sqr,
                         const double temperature);

/** samples a momentum from the Maxwell-Boltzmann distribution
 *
 * \todo rename this function to make clear what it is
 *
 * \see density_integrand
 * \param[in] temperature Temperature \f$T\f$
 * \param[in] mass Mass of the particle: \f$m = \sqrt{E^2 - p^2}\f$
 *
 * \return one possible momentum between mass and 50 GeV
 *
 * \fpPrecision Why \c double?
 */
double sample_momenta(const double temperature, const double mass);

double sample_momenta1(const double temperature, const double mass);



/**\param[in] lam +/-1 or 0 to determine the distribution type
 * lam=0,  for juttner distribution
 * lam=-1, for bose-einstein distribution
 * lam=1,  for fermi-dirac distribution
 * juttner_distribution_func is used in Adaptive rejection sampling to sample
 * momentum_radial
 * Notice: "const double &" is dangerous here because if we define:
 * float lam = 1.0;
 * and use juttner_distribution_func( p, m, T, muB, lam );
 * the lam will be not correct.
 */
double juttner_distribution_func(const double momentum_radial, \
        const double mass, const double temperature, const double  \
        baryon_chemical_potential, const double lam);



/** wookds-saxon distribution function  */
double woods_saxon_dist_func(const double r,  const double radius, \
        const double diffusion);


}  // namespace Smash

#endif  // SRC_INCLUDE_DISTRIBUTIONS_H_
