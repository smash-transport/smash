/*
 *    Copyright (c) 2013-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_SMASH_DISTRIBUTIONS_H_
#define SRC_INCLUDE_SMASH_DISTRIBUTIONS_H_

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
 * Relativistic Juttner distribution function is just a convenience wrapper for
 * displaying Fermi, Bose, and Boltzmann distributions in one mathematical form.
 * \param[in] momentum_radial length of the momentum vector [GeV]
 * \param[in] mass (pole) mass of the particle species [GeV]
 * \param[in] temperature temperature of the system [GeV]
 * \param[in] effective_chemical_potential effective chemical potential of
 *            the system [GeV]
 * \param[in] statistics quantum statistics of the particles species
 *            (+1 for Fermi, -1 for Bose, 0 for Boltzmann)
 * \return the Juttner distribution function
 */
double juttner_distribution_func(double momentum_radial, double mass,
                                 double temperature,
                                 double effective_chemical_potential,
                                 double statistics);

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
 * in a faster way, given by Scott Pratt (see \iref{Pratt:2014vja})
 * APPENDIX: ALGORITHM FOR GENERATING PARTICLES
 * math trick: for \f$ x^{n-1}e^{-x} \f$ distribution, sample x by:
 * \f$ x = -ln(r_1 r_2 r_3 ... r_n) \f$
 * where \f$ r_i \f$ are uniform random numbers between [0,1)
 * for \f$ T/m > 0.6 \f$: \f$ p^2 e^{-E/T} = p^2 e^{-p/T} * e^{(p-E)/T} \f$,
 * where \f$ e^{(p-E)/T}\f$ is used as rejection weight.
 * Since \f$T/m > 0.6 \f$, \f$ e^{(p-E)/T}\f$ is close to 1.
 * for \f$ T/m < 0.6 \f$, there are many rejections
 * another manipulation is used:
 * \f$ p^2 e^{-E/T} dp = dE \frac{E}{p} p^2 e^{-E/T} \f$
 * \f$ = dK \frac{p}{E} (K+m)^2 e^{-K/T} e^{-m/T} \f$
 * \f$ = dK (K^2 + 2mK + m^2) e^{-K/T} \frac{p}{E}\f$
 *  where \f$ \frac{p}{E} \f$ is used as rejection weight.
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

#endif  // SRC_INCLUDE_SMASH_DISTRIBUTIONS_H_
