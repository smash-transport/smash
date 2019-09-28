/*
 *
 *    Copyright (c) 2013-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_AGNIESZKA_DISTRIBUTION_H_
#define SRC_INCLUDE_AGNIESZKA_DISTRIBUTION_H_


// agnieszka included the following libraries:

#include <cmath>

#include <gsl/gsl_multiroots.h>


namespace agnieszka {

  /* 
   * ****************************************************************************
   *
   * This block is for:
   * double juttner_distribution_func
   *
   * The name of this function is the same as in distributions.cc;
   * however, the distribution defined here is easier to handle numerically
   * than the one defined originally in SMASH. This is because 
   * exp(-x)/[ 1 + exp(-x) ] is numerically more stable than 1/[ exp(x) + 1 ] .
   *
   * ****************************************************************************
   */

  /**
   * Juttner distribution function is just a convenience wrapper for displaying
   * Fermi, Bose, and Boltzmann distributions in one mathematical form. 
   * \param[in] momentum_radial the length of the momentum vector
   * \param[in] mass (pole) mass of the particle species
   * \param[in] temperature temperature of the system
   * \param[in[ effective_chemical_potential effective chemical potential of 
   *            the system
   * \param[in] statistics quantum statistics of the particles species 
   *            (+1 for Fermi, -1 for Bose, 0 for Boltzmann)
   * \return the Juttner distribution function
   */
double juttner_distribution_func(double momentum_radial,
                                 double mass,
				 double temperature,
                                 double effective_chemical_potential,
                                 double statistics);

  /* 
   * ****************************************************************************
   *
   * This block is for:
   * Root equations and GSL procedure for finding the momentum for which the
   * maximum of a given Juttner distribution occurs. This is needed for a method 
   * of sampling the distribution function in which one samples uniformly below
   * the maximum of the distribution.
   *
   * ****************************************************************************
   */

  /**
   * Struct object that holds the parameters relevant to finding the momentum
   * for which the maximum of the distribution occurs.
   */
  struct ParametersForMaximumMomentumRootFinder
  {
    double mass;
    double temperature;
    double effective_chemical_potential;
    double statistics; 
  };

  
  /**
   * Root equation for finding the radial momentum at which the Juttner 
   * distribution function has its maximum.
   * \param[in] p radial momentum, i.e., the length of the momentum vector
   * \param[in] mass (pole) mass of the particle species
   * \param[in] temperature temperature of the system
   * \param[in[ effective_chemical_potential effective chemical potential of 
   *            the system 
   * \param[in] statistics quantum statistics of the particles species 
   *            (+1 for Fermi, -1 for Bose, 0 for Boltzmann)
   * \return the extremum equation for the maximum of the Juttner distribution
   */
  double p_max_root_equation
  (double p,
   double mass,
   double temperature,
   double effective_chemical_potential, 
   double statistics);

  /**
   * Root equation for finding the radial momentum at which the Juttner 
   * distribution function has its maximum, suited for the GSL root finding
   * procedure.
   * \param[in] roots_array an array holding the current best estimate of the 
   *            roots of the solved equation
   * \param[in] parameters refers to the parameters as provided in the GSL 
   *            root solving procedure
   * \param[in] function refers to the root equation(s) as provided in the GSL 
   *            root solving procedure (where it's called "function")
   * \return the root equation suited for GSL root finding procedure
   */
  int p_max_root_equation_for_GSL
  (const gsl_vector * roots_array,
   void * parameters,
   gsl_vector * function) ;

  /**
   * A GSL utility which allows for printing out the status of the solver
   * during the root finding procedure.
   * \param[in] iter variable keeping track of how many steps in the root
   *            solving procedure have been taken
   * \param[in] solver GSL solver object, which has acces to the current best
   *            estimate of the roots and the corresponding function values
   * \return message about the current state of the solver
   */
  void cout_state_p_max
  (unsigned int iter,
   gsl_multiroot_fsolver * solver);

  /**
   * A GSL root solver for finding the radial momentum value at which the
   * maximum of the given Juttner distribution function occurs. For the value of
   * the distribution at the maximu, one shoud use the function 
   * maximum_of_the_distribution(). 
   * \param[in] mass (pole) mass of the particle species
   * \param[in] temperature temperature of the system
   * \param[in[ effective_chemical_potential effective chemical potential of 
   *            the system
   * \param[in] statistics quantum statistics of the particles species 
   *            (+1 for Fermi, -1 for Bose, 0 for Boltzmann)
   * \param[in] p_max_initial_guess the initial guess for the value of the 
   *            solution
   * \param[in] solution_precision used for the precision with which the 
   *            solution is found
   * \param[out] p_max the solution (momentum for which the distribution takes 
   *             on the maximum value) stored in an array object
   */
  int find_maximum_of_the_distribution
  (double mass,
   double temperature,
   double effective_chemical_potential,
   double statistics,
   double p_max_initial_guess,
   double solution_precision,
   double* p_max);


  /**
   * A convenience wrapper for find_maximum_of_the_distribution(), yielding the
   * maximum value of the distribution corresponding to the momentum at maximum
   * found by find_maximum_of_the_distribution().
   * \param[in] mass (pole) mass of the particle species
   * \param[in] temperature temperature of the system
   * \param[in[ effective_chemical_potential effective chemical potential of 
   *            the system
   * \param[in] statistics quantum statistics of the particles species 
   *            (+1 for Fermi, -1 for Bose, 0 for Boltzmann)
   * \param[in] solution_precision used for the precision with which the 
   *            solution is found
   */
  double maximum_of_the_distribution
  (double mass,
   double temperature,
   double effective_chemical_potential,
   double statistics,
   double solution_precision);

}  // namespace agnieszka

#endif  // SRC_INCLUDE_AGNIESZKA_DISTRIBUTION_H_
