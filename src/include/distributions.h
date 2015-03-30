/*
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_DISTRIBUTIONS_H_
#define SRC_INCLUDE_DISTRIBUTIONS_H_

#include <gsl/gsl_sf_bessel.h>
#include <cmath>

#include "constants.h"

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
 * \param[in] momentum \f$p\f$ (in GeV)
 * \param[in] temperature \f$T\f$ (in GeV)
 *
 * \return \f$4\pi p^2 \exp{-\frac{E}{T}}\f$
 *
 * \fpPrecision Why \c double?
 */
double density_integrand(const double energy, const double momentum,
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

/** return number density from a Maxwell-Boltzmann distribution
 *
 * \todo rename this function to make clear what it is
 *
 * \see density_integrand
 *
 * \param[in] temperature Temperature
 * \param[in] mass Mass of the particle: \f$m = \sqrt{E^2 - p^2}\f$
 *
 * \return \f$\frac{1}{2\pi^2} m^2 T K_2\left(\frac{m}{T}\right)\f$
 *
 * where \f$K_2\f$ is the modified Bessel function of the second kind.
 *
 * \fpPrecision Why \c double?
 */
inline double number_density_maxwellboltzmann(double mass, double temperature) {
  return mass * mass * temperature * gsl_sf_bessel_Knu(2, mass / temperature)
    * 0.5 * M_1_PI * M_1_PI / hbarc / hbarc / hbarc;
}

}  // namespace Smash

#endif  // SRC_INCLUDE_DISTRIBUTIONS_H_
