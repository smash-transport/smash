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


}  // namespace Smash

#endif  // SRC_INCLUDE_DISTRIBUTIONS_H_
