/*
 *
 *    Copyright (c) 2013-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/chemical_potential.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_multiroots.h>
#include "smash/constants.h"
#include "smash/distributions.h"
#include "smash/integrate.h"

namespace smash {

/*
 * This block is for:
 * Calculating the vector number density of one particle species from
 * integrating the distribution function over the momentum space. The procedure
 * is needed for finding the effective chemical potential for that species. The
 * block includes the auxiliary functions and the integration itself.
 */

double density_integration_one_species_unit_range(
    double degeneracy, double mass, double temperature,
    double effective_chemical_potential, double statistics, double precision) {
  Integrator integrator = Integrator(1E6);
  // Integration precision should be better than solving precision
  const double integration_precision = 0.1 * precision;
  integrator.set_precision(integration_precision, integration_precision);
  const auto integral = integrator(0.0, 1.0, [&](double x)
    {
      return ((1.0 - x) * (1.0 - x) / (x * x * x * x)) *
	juttner_distribution_func((1.0 - x) / x, mass, temperature,
				  effective_chemical_potential, statistics);
    });
  const double integrated_number_density =
    ( degeneracy/(2.0 * M_PI * M_PI) ) * integral;

  return integrated_number_density;
}

/*
 * This block is for:
 * Root equations, and procedure for finding the effective chemical potential
 * for a given particle species. This chemical potential is NOT the equilibrium
 * chemical potential of a system with multiple particle species. If we were to
 * calculate the equilibrium chemical potential in a system  where there are
 * multiple particle species, we would have to use the densities (baryonic as
 * well as strangeness) of all particle species present in the system.
 * As it is, this is still a good estimate for sampling the particles. In
 * particular, the procedure is exact for a system composed of protons and
 * neutrons only, as their mass is degenerate.
 *
 * TODO: Give more of an explanation?
 * TODO2: Sometime in the future, attempt to implement calculating the real
 * chemical potential (it would involve summing over all particle species
 * present, which probably would mean passing the initial multiplicities etc.).
 */

double root_equation_effective_chemical_potential(
    double degeneracy, double mass, double number_density, double temperature,
    double effective_chemical_potential, double statistics, double precision) {
  
  const double integrated_number_density =
    density_integration_one_species_unit_range(
      degeneracy, mass, temperature, effective_chemical_potential, statistics,
      precision);

  return number_density - integrated_number_density;
}

int root_equation_effective_chemical_potential_for_GSL(
    const gsl_vector *roots_array, void *parameters, gsl_vector *function) {
  struct ParametersForEffectiveChemicalPotentialRootFinderOneSpecies *par =
      static_cast<
          ParametersForEffectiveChemicalPotentialRootFinderOneSpecies *>(
          parameters);

  const double degeneracy = (par->degeneracy);
  const double mass = (par->mass);
  const double number_density = (par->number_density);
  const double temperature = (par->temperature);
  const double statistics = (par->statistics);
  const double precision = (par->precision);

  const double effective_chemical_potential = gsl_vector_get(roots_array, 0);

  gsl_vector_set(function, 0,
                 root_equation_effective_chemical_potential(
                     degeneracy, mass, number_density, temperature,
                     effective_chemical_potential, statistics, precision));

  return GSL_SUCCESS;
}

void cout_state_effective_chemical_potential(unsigned int iter,
                                             gsl_multiroot_fsolver *solver) {
  printf(
      "\n***\nfind_effective_chemical_potential(): iter = %3u \t"
      "x = % .3f \t"
      "f(x) = % .3e \n",
      iter, gsl_vector_get(solver->x, 0), gsl_vector_get(solver->f, 0));
}

int find_effective_chemical_potential(double degeneracy, double mass,
                                      double number_density, double temperature,
                                      double statistics,
                                      double mu_initial_guess,
                                      double solution_precision,
                                      double *effective_chemical_potential) {
  const gsl_multiroot_fsolver_type *Solver_name;
  gsl_multiroot_fsolver *Root_finder;

  int status;
  size_t iter = 0;

  const size_t problem_dimension = 1;

  struct ParametersForEffectiveChemicalPotentialRootFinderOneSpecies
      parameters = {degeneracy,  mass,       number_density,
                    temperature, statistics, solution_precision};

  gsl_multiroot_function EffectiveChemicalPotential = {
      &root_equation_effective_chemical_potential_for_GSL, problem_dimension,
      &parameters};

  double roots_array_initial[problem_dimension] = {mu_initial_guess};

  gsl_vector *roots_array = gsl_vector_alloc(problem_dimension);
  gsl_vector_set(roots_array, 0, roots_array_initial[0]);

  Solver_name = gsl_multiroot_fsolver_hybrids;
  Root_finder = gsl_multiroot_fsolver_alloc(Solver_name, problem_dimension);
  gsl_multiroot_fsolver_set(Root_finder, &EffectiveChemicalPotential,
                            roots_array);

  // cout_state_effective_chemical_potential (iter, Root_finder);

  do {
    iter++;

    status = gsl_multiroot_fsolver_iterate(Root_finder);

    // cout_state_effective_chemical_potential (iter, Root_finder);

    /*
     * Check if the solver is stuck
     */
    if (status) {
      const double conversion_factor =
          smash::hbarc * smash::hbarc * smash::hbarc;
      std::cout << "\n\nThe GSL solver"
                << "\nfind_effective_chemical_potential"
                << "\nis stuck!"
                << "\n\nInput parameters:"
                << "\n            degeneracy = " << degeneracy
                << "\n            mass [GeV] = " << mass
                << "\nnumber_density [fm^-3] = "
                << number_density / (conversion_factor)
                << "\n     temperature [GeV] = " << temperature
                << "\n            statistics = " << statistics
                << "\n    solution_precision = " << solution_precision << "\n\n"
                << std::endl;
      std::cout << "The program cannot continue without calculating the "
                << "chemical potential. \nTry adjusting the initial "
                << "guess or the solution precision.\n\n\n"
                << std::endl;
      throw std::runtime_error(
          "Chemical potential calculation returned no result.\n\n");
      continue;
    }

    status = gsl_multiroot_test_residual(Root_finder->f, solution_precision);

    if (status == GSL_SUCCESS) {
      effective_chemical_potential[0] = gsl_vector_get(Root_finder->x, 0);
      /*
       * Make sure that the calculated chemical potential reproduces
       * the input number density within accepted error.
       */
      double nDensityCheck = density_integration_one_species_unit_range(
          degeneracy, mass, temperature, effective_chemical_potential[0],
          statistics, solution_precision);

      if (abs((number_density - nDensityCheck)) > 1e-4) {
        const double conversion_factor =
            smash::hbarc * smash::hbarc * smash::hbarc;
        std::cout << "\n\nThe calculated chemical potential does not "
                  << "reproduce the input number density!"
                  << "\nnumber_density [GeV^3] = " << number_density
                  << "\n nDensityCheck [GeV^3] = " << nDensityCheck
                  << "\nnumber_density [fm^-3] = "
                  << number_density / (conversion_factor)
                  << "\n nDensityCheck [fm^-3] = "
                  << nDensityCheck / (conversion_factor) << "\n\n"
                  << std::endl;
        std::cout << "Input parameters:"
                  << "\n            degeneracy = " << degeneracy
                  << "\n            mass [GeV] = " << mass
                  << "\nnumber_density [fm^-3] = "
                  << number_density / (conversion_factor)
                  << "\n     temperature [GeV] = " << temperature
                  << "\n            statistics = " << statistics
                  << "\n    solution_precision = " << solution_precision
                  << "\n\n"
                  << std::endl;
        std::cout << "The program cannot continue without a correct "
                  << "calculation of the chemical potential. \nTry "
                  << "adjusting the initial guess or the solution "
                  << "comparison precision."
                  << "\n(in chemical_potential.cc)\n\n\n"
                  << std::endl;
        throw std::runtime_error(
            "Chemical potential calculation does not reproduce the "
            "input number density.\n\n");
      }
    }

  } while (status == GSL_CONTINUE && iter < 100000);

  gsl_multiroot_fsolver_free(Root_finder);
  gsl_vector_free(roots_array);

  return 0;
}

/*
 * This block is for:
 * A convenience wrapper for finding the effective chemical potential for a
 * given particle species and performing sanity checks.
 */

double effective_chemical_potential(double degeneracy, double mass,
                                    double number_density, double temperature,
                                    double statistics,
                                    double solution_precision) {
  /*
   * The initial guess for the GSL chemical potential finding
   * procedure; the guess has to be different for Bose and Fermi
   * particles.
   */
  double initial_guess = 0.0;
  if (statistics == -1.0) {
    /*
     * Such initial guess allows one to sample really close to the Bose
     * condensate region; at the same time, it works fine in other regions as
     * well.
     */
    initial_guess = 0.999999 * mass;
  } else {
    initial_guess = 1.05 * mass;
  }

  /*
   * Array for the GSL chemical potential finder.
   * Use of array is justified here as we're planning on adding scalar
   * interactions in the (distant) future, which will lead to multidimensional
   * root finding; then the result is best passed as an array.
   */
  double chemical_potential[1];
  chemical_potential[0] = 0.0;

  /*
   * Call the GSL effective chemical potential finding procedure, which yields
   * the value of the effective chemical potential.
   */
  find_effective_chemical_potential(degeneracy, mass, number_density,
                                    temperature, statistics, initial_guess,
                                    solution_precision, chemical_potential);

  /*
   * Sanity checks are performed for the obtained value of the chemical
   * potential .
   */

  /*
   * Check if the obtained value of the effective chemical potential indicates
   * that the number density is such that for bosons one has encountered the
   * Bose-Einstein condensate. We need mu < mass for the Bose distribution NOT
   * to be in the Bose-Einstein condensate region.
   */
  if ((statistics == -1) && (chemical_potential[0] > mass)) {
    const double conversion_factor = smash::hbarc * smash::hbarc * smash::hbarc;
    std::cout << "\n\nThe calculated chemical potential indicates the "
              << "input number density is such that Bose-Einstein "
              << "condensate is encountered." << std::endl;
    std::cout << "\n\nInput parameters:"
              << "\n            degeneracy = " << degeneracy
              << "\n            mass [GeV] = " << mass
              << "\nnumber_density [fm^-3] = "
              << number_density / (conversion_factor)
              << "\n     temperature [GeV] = " << temperature
              << "\n            statistics = " << statistics
              << "\n    solution_precision = " << solution_precision << "\n\n"
              << std::endl;
    std::cout << "The program cannot proceed with a pathological "
              << "distribution function leading to the Bose-Einstein "
              << "condensate. \nTry increasing the temperature or "
              << "decreasing the number density.\n\n\n"
              << std::endl;
    throw std::runtime_error(
        "Chemical potential calculation indicates that the Bose-Einstein "
        "condensate is produced.");
  }

  /*
   * Check if the distribution is pathological, i.e., has negative values
   * for low momenta (Bose statistics, in the Bose-Einstein condensate regime,
   *  may produce pathological distributions for which f(p) < 0 ).
   */
  /// value of small_p in GeV; value is chosen to be small;
  /// possibly an arbitrary choice
  const double small_p = 0.005;
  const double distribution_at_small_p = juttner_distribution_func(
      small_p, mass, temperature, chemical_potential[0], statistics);
  if (distribution_at_small_p < 0.0) {
    const double conversion_factor = smash::hbarc * smash::hbarc * smash::hbarc;
    std::cout << "\n\nNegative values of the distribution function at "
              << "small momenta indicate that Bose-Einstein "
              << "condensate is encountered." << std::endl;
    std::cout << "\n\nInput parameters:"
              << "\n            degeneracy = " << degeneracy
              << "\n            mass [GeV] = " << mass
              << "\nnumber_density [fm^-3] = "
              << number_density / (conversion_factor)
              << "\n     temperature [GeV] = " << temperature
              << "\n            statistics = " << statistics
              << "\n    solution_precision = " << solution_precision
              << "\n\nf(p=" << small_p << ") = " << distribution_at_small_p
              << "\n\n"
              << std::endl;
    std::cout << "The program cannot proceed with a pathological "
              << "distribution function leading to the Bose-Einstein "
              << "condensate. \nTry increasing the temperature or "
              << "decreasing the number density.\n\n\n"
              << std::endl;
    throw std::runtime_error(
        "Negative values of the distribution function at low momenta "
        "indicate that the Bose-Einstein condensate is produced.");
  }

  double effective_chemical_potential = chemical_potential[0];
  return effective_chemical_potential;
}

}  // namespace smash
