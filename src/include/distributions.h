/*
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_DISTRIBUTIONS_H_
#define SRC_INCLUDE_DISTRIBUTIONS_H_

namespace Smash {

/**
 * Returns a relativistic Breit-Wigner distribution. The normalization is such
 * that integrating over srts from 0 to inf yields one.
 *
 * \param[in] m Argument of the Breit-Wigner function (off-shell mass m in GeV)
 * \param[in] pole resonance pole mass \f$ m_0 \f$ in GeV
 * \param[in] width resonance width \f$ \Gamma \f$ in GeV
 *
 * \return \f$ \frac{2}{\pi} \frac{m^2\Gamma}{(m^2-m_0^2)^2 + m^2\Gamma^2} \f$
 */
float breit_wigner(float m, float pole, float width);

/**
 * Returns a non-relativistic Breit-Wigner distribution, which is essentially
 * a Cauchy distribution with half width.
 *
 * \param[in] m Argument of the Breit-Wigner function (off-shell mass m in GeV)
 * \param[in] pole resonance pole mass \f$ m_0 \f$ in GeV
 * \param[in] width resonance width \f$ \Gamma \f$ in GeV
 *
 * \return \f$ \frac{\Gamma/2}{\pi ((m-m_0)^2+\Gamma^2/4)}\f$
 */
float breit_wigner_nonrel(float m, float pole, float width);

/**
 * Returns a Cauchy distribution (sometimes also called Lorentz or
 * non-relativistic Breit-Wigner distribution) with the given parameters.
 * The normalization is such that integrating over x from -inf to inf yields one.
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
