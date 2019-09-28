/*
 *
 *    Copyright (c) 2013-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */


#include <cmath>
#include <cstdlib>
#include <iostream>

#include <random>

// agnieszka included the following libraries:

#include "./include/smash/Agnieszka_distribution.h"
#include "./include/smash/Agnieszka_chemical_potential.h"

#include "./include/smash/Agnieszka_sampling.h"

// needed for solving for the distribution maximum
// and the chemical potential
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_integration.h>
// needed to handle GSL errors
//#include <gsl/gsl_vector.h>

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
   * This sampler is the simplest implementation of sampling based on sampling 
   * from a uniform distribution. [elaborate?]
   */
  double sample_momenta_from_Juttner
  (double mass,
   double temperature,
   double effective_chemical_potential, 
   double statistics,
   double p_range,
   double distribution_maximum)
  {
    std::random_device rd;
    std::mt19937 generate(rd());

    double sampled_momentum = 0.0;
    bool success = false;

    int iter = 0;

    while ( !success )
      {
	iter++;
	
        sampled_momentum = p_range *
	  std::generate_canonical<double, 1>(generate);
	double distribution_at_sampled_p =
	  sampled_momentum * sampled_momentum *
	  juttner_distribution_func(sampled_momentum,
				    mass,
				    temperature,
				    effective_chemical_potential,
				    statistics);
	double sampled_ratio = distribution_at_sampled_p/distribution_maximum;
	double accept_or_reject = std::generate_canonical<double, 1>(generate);


	if ( sampled_ratio > accept_or_reject )
	  {
	    success = true;
	  }
      }

    return sampled_momentum;
    
  }

  

}  // namespace agnieszka
