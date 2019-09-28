/*
 *
 *    Copyright (c) 2013-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */


#ifndef SRC_INCLUDE_AGNIESZKA_SAMPLING_H_
#define SRC_INCLUDE_AGNIESZKA_SAMPLING_H_


// agnieszka included the following libraries:

#include <cmath>

// needed for solving for the distribution maximum
// and the chemical potential
#include <gsl/gsl_multiroots.h>
// needed to handle GSL vectors which are in declarations of functions
#include <gsl/gsl_vector.h>

namespace agnieszka {

  /**
   * ****************************************************************************
   *
   * This block is for:
   * Sampling radial momenta of given particle species from Bose, Boltzmann, or
   * Fermi distribution. The choice between the distributions is made based on 
   * the <statistics> variable:
   * -1 for Bose
   *  0 fot Boltzmann
   * +1 for Fermi
   * The information about degeneracy is not needed for sampling. The chemical
   * potential associated with the given particle species must be calculated 
   * before sampling.
   *
   * ****************************************************************************
   */

  /**
   * Sampler for radial momenta based on a given Juttner distribution function.
   * \param[in] mass (pole) mass of the particle species
   * \param[in] temperature temperature of the system
   * \param[in[ effective_chemical_potential effective chemical potential of 
   *            the system
   * \param[in] statistics quantum statistics of the particles species 
   *            (+1 for Fermi, -1 for Bose, 0 for Boltzmann)
   * \param[in] p_range the maximum value of the radial momentum sampled
   * \param[in] distribution_maximum the value of the distribution function at
   *            its maximum
   * return sampled momentum
   */
  double sample_momenta_from_Juttner
  (double mass,
   double temperature,
   double effective_chemical_potential,
   double statistics,
   double p_range,
   double distribution_maximum);


  
}  // namespace agnieszka

#endif  // SRC_INCLUDE_AGNIESZKA_SAMPLING_H_
