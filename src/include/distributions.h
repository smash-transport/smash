/*
 *    Copyright (c) 2013-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_DISTRIBUTIONS_H_
#define SRC_INCLUDE_DISTRIBUTIONS_H_

namespace smash {

/**
 * Returns a relativistic Breit-Wigner distribution. The normalization is such
 * that the integration over \f$ \sqrt{s} \f$ from 0 to infinity yields one.
 *
 * \param[in] m Argument of the Breit-Wigner function (off-shell mass m in GeV)
 * \param[in] pole Resonance pole mass \f$ m_0 \f$ in GeV
 * \param[in] width Resonance width \f$ \Gamma \f$ in GeV
 *
 * \return \f$ \frac{2}{\pi} \frac{m^2\Gamma}{(m^2-m_0^2)^2 + m^2\Gamma^2} \f$
 */
double breit_wigner(double m, double pole, double width);

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
double breit_wigner_nonrel(double m, double pole, double width);

/**
 * Returns a Cauchy distribution (sometimes also called Lorentz or
 * non-relativistic Breit-Wigner distribution) with the given parameters.
 * The normalization is such that integrating over x from -infinity to
 * +infinity yields one.
 *
 * \param x Argument of the Cauchy function in GeV.
 * \param pole Pole parameter \f$ m_0 \f$ of the Cauchy function in GeV,
 *  i.e. location of the peak.
 * \param width Width parameter \f$ \Gamma \f$ of the Cauchy function in GeV,
 * determining the sharpness of the peak.
 *
 * \return \f$ \frac{\Gamma}{\pi ((x-m_0)^2+\Gamma^2)}\f$
 */
double cauchy(double x, double pole, double width);

/**
 * 1D Woods-Saxon distribution function
 * \todo(duplicate) Check nucleus and remove one of the Woods-Saxon
 * implementations, this one is actually only used in the test of the
 * adaptive rejection sampler
 * \param[in] r Radial coordinates in the nucleus in units of [fm]
 * \param[in] radius Radius parameter of the nucleus in units of [fm]
 * \param[in] diffusion Diffusiveness parameter of the nucleus in units of [fm]
 * \return Unnormalized nucleon density in units [fm^{-3}]
 * in the nucleus at r
 */
double woods_saxon_dist_func(const double r, const double radius,
                             const double diffusion);

/**
 * Returns the Maxwell-Boltzmann distribution
 *
 * \todo rename the following 4 functions to make clear what they are
 * and check if they are actually used
 *
 * \param[in] energy \f$E\f$ (in GeV)
 * \param[in] momentum_sqr squared \f$p\f$ (in GeV\f$^2\f$)
 * \param[in] temperature \f$T\f$ (in GeV)
 *
 * \return \f$4\pi p^2 \exp{-\frac{E}{T}}\f$
 */
double density_integrand(const double energy, const double momentum_sqr,
                         const double temperature);
/**
 * density_integrand_mass - off_equilibrium distribution for massive
 * particles
 *
 * \param[in] energy \f$E\f$ (in GeV)
 * \param[in] momentum_sqr squared \f$p\f$ (in GeV\f$^2\f$)
 * \param[in] temperature \f$T\f$ (in GeV)
 * \return \f[f=pe^{-\frac{\sqrt{m^2+p^2}}{T_0}}\f]
 */
double density_integrand_mass(const double energy, const double momentum_sqr,
                              const double temperature);
/**
 * density integrand - 1M_IC massless particles for expanding metric
 * initialization, see \iref{Bazow:2016oky}
 *
 * \param[in] energy \f$E\f$ (in GeV)
 * \param[in] momentum_sqr squared \f$p\f$ (in GeV\f$^2\f$)
 * \param[in] temperature \f$T\f$ (in GeV)
 * \return Value of function 1M_IC
 */
double density_integrand_1M_IC(const double energy, const double momentum_sqr,
                               const double temperature);

/**
 * density integrand - 2M_IC massless particles for expanding metric
 * initialization, see \iref{Bazow:2016oky}
 *
 * \param[in] energy \f$E\f$ (in GeV)
 * \param[in] momentum_sqr squared \f$p\f$ (in GeV\f$^2\f$)
 * \param[in] temperature \f$T\f$ (in GeV)
 * \return Value of function 2M_IC
 */
double density_integrand_2M_IC(const double energy, const double momentum_sqr,
                               const double temperature);

/**
 * Relativistic Juttner distribution function
 *
 * \param[in] momentum_radial \f$|\vec{p}|\f$ in units of [GeV]
 * \param[in] mass Mass of the particle: in units of [GeV]
 * \param[in] temperature Temperature of the system \f$T\f$ in units of [GeV]
 * \param[in] baryon_chemical_potential \f$n*\mu_{B}\f$ default = 0
 * \param[in] lam +/-1 or 0 to determine the distribution type
 * lam=0,  for Juttner distribution
 * lam=-1, for Bose-Einstein distribution
 * lam=1,  for Fermi-Dirac distribution
 * \return Unnormalized probability of relativistic Juttner distribution
 * \todo(unused) this is only used in the test of adaptive rejection sampler
 */
double juttner_distribution_func(const double momentum_radial,
                                 const double mass, const double temperature,
                                 const double baryon_chemical_potential,
                                 const double lam);
/**
 * Samples a momentum via rejection method from the non-equilibrium
 * distribution
 * \f[f=pe^{-\frac{\sqrt{m^2+p^2}}{T_0}}\f]
 *
 * \see density_integrand_mass
 * \param[in] temperature Temperature \f$T\f$ [GeV]
 * \param[in] mass Mass of the particle: \f$m = \sqrt{E^2 - p^2}\f$ [GeV]
 *
 * \return one possible momentum between mass and 50 GeV
 */
double sample_momenta_non_eq_mass(const double temperature, const double mass);

/**
 * Samples a momentum from the non-equilibrium distribution
 * 1M_IC from \iref{Bazow:2016oky}
 *
 * \see density_integrand_1M_IC
 * \param[in] temperature Temperature \f$T\f$ [GeV]
 * \param[in] mass Mass of the particle: \f$m = \sqrt{E^2 - p^2}\f$ [GeV]
 *
 * \return one possible momentum between mass and 50 GeV
 */
double sample_momenta_1M_IC(const double temperature, const double mass);

/**
 * Samples a momentum from the non-equilibrium distribution
 * 2M_IC from \iref{Bazow:2016oky}
 *
 * \see density_integrand_2M_IC
 * \param[in] temperature Temperature \f$T\f$ [GeV]
 * \param[in] mass Mass of the particle: \f$m = \sqrt{E^2 - p^2}\f$ [GeV]
 *
 * \return one possible momentum between mass and 50 GeV
 */
double sample_momenta_2M_IC(const double temperature, const double mass);

/**
 * Samples a momentum from the Maxwell-Boltzmann (thermal) distribution
 * in a faster way, given by Scott Pratt
 *
 * \param[in] temperature Temperature \f$T\f$ [GeV]
 * \param[in] mass Mass of the particle: \f$m = \sqrt{E^2 - p^2}\f$ [GeV]
 *
 * \return one possible momentum
 */
double sample_momenta_from_thermal(const double temperature, const double mass);

/**
 * Sample momenta according to the momentum distribution
 * in \iref{Bazow:2016oky}
 *
 * \param[in] temperature The temperature for the distribution [GeV]
 * \return Radial momentum
 */
double sample_momenta_IC_ES(const double temperature);
}  // namespace smash

#endif  // SRC_INCLUDE_DISTRIBUTIONS_H_
