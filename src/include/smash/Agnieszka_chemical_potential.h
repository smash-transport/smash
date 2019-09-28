/*
 *
 *    Copyright (c) 2013-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_AGNIESZKA_CHEMICAL_POTENTIAL_H_
#define SRC_INCLUDE_AGNIESZKA_CHEMICAL_POTENTIAL_H_



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
   * Calculating the vector number density of one particle species from 
   * integrating the distribution function over the momentum space. The procedure
   * is needed for finding the effective chemical potential for that species. The
   * block includes the struct that holds the parameters, auxiliary functions, 
   * and the integration itself.
   *
   * ****************************************************************************
   */

  /**
   * Struct object that holds the parameters relevant to momentum integrals
   * where the integration is done for one particle species.
   */
  struct ParametersForThermalMomentumIntegralsOneSpecies
  {
    double degeneracy;
    double mass;
    double temperature;
    double effective_chemical_potential;
    double statistics;
  };

  /**
   * Vector number density integrand for obtaining density of one particle 
   * species; integrand written for integration on the (0,1) interval.
   * \param[in] x integration variable
   * \param[in] degeneracy degeneracy of the particle species
   * \param[in] mass (pole) mass of the particle species
   * \param[in] temperature temperature of the system
   * \param[in] effective_chemical_potential effective chemical potential of 
   *            the system
   * \param[in] statistics quantum statistics of the particles species 
   *            (+1 for Fermi, -1 for Bose, 0 for Boltzmann)
   * \return the complete integrand for the vector number density integration
   */
  double density_integrand_one_species_unit_range
  (double x,
   double degeneracy,
   double mass,
   double temperature,
   double effective_chemical_potential, 
   double statistics);


  /**
   * Vector number density integrand for obtaining density of one particle 
   * species as used in the GSL integration procedure (which requires that the 
   * integrand depends explicitly only on the <integration variable> and 
   * <parameters>).
   * \param[in] x integration variable
   * \param[in] parameters parameters provided in the struct/integration
   * \return the complete integrand for the vector number density integration 
   */
  double density_integrand_one_species_unit_range_for_GSL
  (double momentum,
   void * parameters);

   /**
   * Vector number density of one particle species, obtained through integrating
   * the distribution function over the momentum space.
   * \param[in] degeneracy degeneracy of the particle species
   * \param[in] mass (pole) mass of the particle species
   * \param[in] temperature temperature of the system
   * \param[in[ effective_chemical_potential effective chemical potential in 
   *            the system
   * \param[in] statistics quantum statistics of the particles species 
   *            (+1 for Fermi, -1 for Bose, 0 for Boltzmann)
   * \param[in] precision used for absolute error, relative error, and the lower
   *            limit of the integration interval
   * \return the vector number density for a given species of particles
   */
  double density_integration_one_species_unit_range
  (double degeneracy,
   double mass,
   double temperature,
   double effective_chemical_potential,
   double statistics,
   double precision);

  

  /** 
   * ****************************************************************************
   *
   * This block is for:
   * Struct, root equations, and procedure for finding the effective chemical 
   * potential for a given particle species. This chemical potential is NOT the 
   * equilibrium chemical potential of a system with multiple particle species.
   * If we were to calculate the equilibrium chemical potential in a system 
   * where there are multiple particle species, we would have to use the 
   * densities (baryonic as well as strangeness) of all particle species present
   * in the system. 
   * As it is, this is still a good estimate for sampling the particles. In
   * particular, the procedure is exact for a system composed of protons and 
   * neutrons only, as their mass is degenerate.
   * 
   * TODO: Give more of an explanation?
   * TODO2: Sometime in the future, attempt to implement calculating the real 
   * chemical potential (it would involve summing over all particle species 
   * present, which probably would mean passing the initial multiplicities etc.).
   *
   * ****************************************************************************
   */

  /**
   * Struct object that holds the parameters relevant to finding the effective
   * chemical potential of one particle species.
   */
  struct ParametersForEffectiveChemicalPotentialRootFinderOneSpecies
  {
    double degeneracy;
    double mass;
    double number_density;
    double temperature;
    double statistics;
    double precision;  
  };

  /**
   * Root equation for finding the value of the effective chemical potential
   * for one particle species.
   * \param[in] degeneracy degeneracy of the particle species
   * \param[in] mass (pole) mass of the particle species
   * \param[in] number_density number density of the particle species
   * \param[in] temperature temperature of the system
   * \param[in[ effective_chemical_potential effective chemical potential in 
   *            the system
   * \param[in] statistics quantum statistics of the particles species 
   *            (+1 for Fermi, -1 for Bose, 0 for Boltzmann)
   * \param[in] precision used for the precision of the number density integral
   * \return the root equation
   */
  double root_equation_for_effective_chemical_potential
  (double degeneracy,
   double mass,
   double number_density,
   double temperature,
   double effective_chemical_potential,
   double statistics,
   double precision);

  /**
   * Root equation for finding the value of the effective chemical potential
   * for one particle species, suited for the GSL root finding procedure. 
   * \param[in] roots_array an array holding the current best estimate of the 
   *            roots of the solved equation
   * \param[in] parameters refers to the parameters as provided in the GSL 
   *            root solving procedure
   * \param[in] function refers to the root equation(s) as provided in the GSL 
   *            root solving procedure (where it's called "function")
   * \return the root equation suited for GSL root finding procedure
   */
  int root_equation_effective_chemical_potential_for_GSL
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
  void cout_state_effective_chemical_potential
  (unsigned int iter,
   gsl_multiroot_fsolver * solver);
  
  /**
   * A GSL root solver for finding the effective chemical potential. In practice
   * one should use the convenience wrapper find_effective_chemical_potential(), 
   * which also performs sanity checks on the obtained solution.
   * \param[in] degeneracy degeneracy of the particle species
   * \param[in] mass (pole) mass of the particle species
   * \param[in] number_density number density of the particle species
   * \param[in] temperature temperature of the system
   * \param[in[ effective_chemical_potential effective chemical potential of 
   *            the system
   * \param[in] statistics quantum statistics of the particles species 
   *            (+1 for Fermi, -1 for Bose, 0 for Boltzmann)
   * \param[in] mu_initial_guess the initial guess for the value of the 
   *            effective chemical potential 
   * \param[in] solution_precision used for the precision with which the 
   *            solution is found; note that this precision also goes into the 
   *            precision of the integrals calculated in the process 
   * \param[out] effective_chemical_potential the solution stored in an array 
   *             object (we use an array for that in anticipation of generalizing 
   *             to multidimensional root finding, needed for example when scalar 
   *             interactions are present and effective mass has to be calculated 
   *             at the same time as the effective chemical potential)
   */
  int find_effective_chemical_potential
  (double degeneracy,
   double mass,
   double number_density,
   double temperature,
   double statistics,
   double mu_initial_guess,
   double solution_precision,
   double* effective_chemical_potential);


  /** 
   * ***************************************************************************
   *
   * This block is for:
   * Convenience wrapper for finding the effective chemical potential for a 
   * given particle species and performing sanity checks.
   *
   * ***************************************************************************
   */
  
  /**
   * The convenience wrapper performs the search for the effective chemical
   * potential, and then performs sanity checks against encountering the
   * Bose-Einstein condensate.
   * \param[in] degeneracy degeneracy of the particle species
   * \param[in] mass (pole) mass of the particle species
   * \param[in] number_density number density of the particle species
   * \param[in] temperature temperature of the system
   * \param[in[ effective_chemical_potential effective chemical potential of 
   *            the system
   * \param[in] statistics quantum statistics of the particles species 
   *            (+1 for Fermi, -1 for Bose, 0 for Boltzmann)
   * \param[in] mu_initial_guess the initial guess for the value of the 
   *            effective chemical potential 
   * \param[in] solution_precision used for the precision with which the 
   *            solution is found; note that this precision also goes into the 
   *            precision of the integrals calculated in the process 
   * \param[out] effective_chemical_potential the solution stored in an array 
   *             object (we use an array for that in anticipation of generalizing 
   *             to multidimensional root finding, needed for example when scalar 
   *             interactions are present and effective mass has to be calculated 
   *             at the same time as the effective chemical potential)
   */
  double effective_chemical_potential
  (double degeneracy,
   double mass,
   double number_density,
   double temperature,
   double statistics,
   double solution_precision);
  
  
}  // namespace agnieszka

#endif  // SRC_INCLUDE_AGNIESZKA_CHEMICAL_POTENTIAL_H_
