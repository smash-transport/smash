/*
 *
 *    Copyright (c) 2020-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_QUANTUMSAMPLING_H_
#define SRC_INCLUDE_SMASH_QUANTUMSAMPLING_H_

#include <map>

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include "smash/pdgcode.h"
#include "smash/random.h"



namespace smash {

/*
 * This block is for:
 * Root equations and GSL procedure for finding the momentum for which the
 * maximum of a given Juttner distribution occurs. This is needed for a method
 * of sampling the distribution function in which one samples uniformly below
 * the maximum of the distribution.
 */
  
/**
 * This class:
 * - Tabulates chemical potential of quantum distribution given density (WRONG)
 * - Tabulates maxima of a Juttner distribution for these chemical potentials
 * - Samples Juttner distribution. This is the main intent of this class,
 *   while previous points are auxiliary calculations for it.
 */

class QuantumSampling {
  public:
    QuantumSampling(const std::map<PdgCode, int> &initial_multiplicities,
		    double volume, double temperature);

    /**
     * Struct object that holds the parameters relevant to finding the momentum
     * for which the maximum of the distribution occurs.
     */
    struct ParametersForMaximumMomentumRootFinder {
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
    static double p_max_root_equation(double p, double mass, double temperature,
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
    static int p_max_root_equation_for_GSL(const gsl_vector* roots_array,
					   void* parameters,
					   gsl_vector* function);

    /**
     * A GSL utility which allows for printing out the status of the solver
     * during the root finding procedure.
     * \param[in] iter variable keeping track of how many steps in the root
     *            solving procedure have been taken
     * \param[in] solver GSL solver object, which has acces to the current best
     *            estimate of the roots and the corresponding function values
     * \return message about the current state of the solver
     */
    static void print_state_p_max(unsigned int iter,
				 gsl_multiroot_fsolver* solver);

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
    static int find_maximum_of_the_distribution
    (double mass, double temperature, double effective_chemical_potential,
     double statistics, double p_max_initial_guess, double solution_precision,
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
    static double maximum_of_the_distribution(double mass, double temperature,
				       double effective_chemical_potential,
				       double statistics,
				       double solution_precision);



    /*
     * Sampling radial momenta of given particle species from Boltzmann, Bose, or
     * Fermi distribution. This sampler uses the simplest rejection sampling.
     */
    double sample(const PdgCode pdg);
  private:
    /// Tabulated effective chemical potentials for every particle species
    std::map<PdgCode, double> effective_chemical_potentials_;
    /// Tabulated distribution function maxima for every particle species
    std::map<PdgCode, double> distribution_function_maximums_;
    /// Volume [fm^3] in which particles sre sampled
    const double volume_;
    /// Temperature [GeV]
    const double temperature_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_QUANTUMSAMPLING_H_
