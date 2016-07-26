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



/** samples a momentum from the Maxwell-Boltzmann (thermal) distribution
 * in a faster way, given by Pratt Scott
 *
 * \param[in] temperature Temperature \f$T\f$
 * \param[in] mass Mass of the particle: \f$m = \sqrt{E^2 - p^2}\f$
 *
 * \return one possible momentum
 *
 * \fpPrecision Why \c double?
 */
double sample_momenta_from_thermal(const double temperature, const double mass);



/** Relativistic Juttner distribution function
 * \param[in] momentum_radial \f$|\vec{p}|\f$ in units of [GeV]
 * \param[in] mass Mass of the particle: in units of [GeV]
 * \param[in] temperature Temperature of the system \f$T\f$ in units of [GeV]
 * \param[in] baryon_chemical_potential \f$n*\mu_{B}\f$ default=0
 * \param[in] lam +/-1 or 0 to determine the distribution type
 * lam=0,  for juttner distribution
 * lam=-1, for bose-einstein distribution
 * lam=1,  for fermi-dirac distribution
 * \return Unnormalized probability of relativistic juttner distribution
 */
double juttner_distribution_func(const double momentum_radial,
        const double mass, const double temperature, const double
        baryon_chemical_potential, const double lam);



/** 1D woods-saxon distribution function
 * \param[in] r Radial coordinates in the nucleus in units of [fm]
 * \param[in] radius Radius parameter of the nucleus in units of [fm]
 * \param[in] diffusion Diffusiveness parameter of the nucleus in units of [fm]
 * \return Unnormalized nucleon density in units [fm^{-3}]
 * in the nucleus at r  */
double woods_saxon_dist_func(const double r,  const double radius,
        const double diffusion);

double sample_noneq_photon_momenta (const double temperature);
}  // namespace Smash

#endif  // SRC_INCLUDE_DISTRIBUTIONS_H_
