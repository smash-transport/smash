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
 */
float breit_wigner(const float mandelstam_s, const float resonance_mass,
                   const float resonance_width);

/**
 * Returns a Cauchy distribution (sometimes also called Lorentz or
 * non-relativistic Breit-Wigner distribution) with the given parameters.
 *
 * \param x Argument of the Cauchy function.
 * \param pole Pole parameter \f$ m_0 \f$ of the Cauchy function, i.e. location of the peak.
 * \param width Width parameter \f$ \Gamma \f$ of the Cauchy function, determining the sharpness of the peak.
 *
 * \return \f$ \frac{\Gamma}{\pi ((m-m_0)^2+\Gamma^2)}\f$
 */
float cauchy(float x, float pole, float width);

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

}  // namespace Smash

#endif  // SRC_INCLUDE_DISTRIBUTIONS_H_
