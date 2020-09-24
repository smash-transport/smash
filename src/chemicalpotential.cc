/*
 *
 *    Copyright (c) 2020-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/chemicalpotential.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_multiroots.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include "smash/constants.h"
#include "smash/distributions.h"
#include "smash/integrate.h"
#include "smash/logging.h"

namespace smash {

/**
 * Vector number density of one particle species.
 */
double ChemicalPotentialSolver::density_one_species(
    double degeneracy, double mass, double temperature,
    double effective_chemical_potential, double statistics,
    Integrator *integrator) {
  const auto integral = (*integrator)(0.0, 1.0, [&](double x) {
    return ((1.0 - x) * (1.0 - x) / (x * x * x * x)) *
           juttner_distribution_func((1.0 - x) / x, mass, temperature,
                                     effective_chemical_potential, statistics);
  });
  const double integrated_number_density =
      (degeneracy / (2.0 * M_PI * M_PI)) * integral;

  return integrated_number_density;
}

/*
 * Root equations and procedure for finding the effective chemical potential
 * for a given particle species.
 */

double ChemicalPotentialSolver::root_equation_effective_chemical_potential(
    double degeneracy, double mass, double number_density, double temperature,
    double effective_chemical_potential, double statistics,
    Integrator *integrator) {
  const double integrated_number_density =
      density_one_species(degeneracy, mass, temperature,
                          effective_chemical_potential, statistics, integrator);

  return number_density - integrated_number_density;
}

int ChemicalPotentialSolver::root_equation_effective_chemical_potential_for_GSL(
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
  Integrator *integrator = (par->integrator);

  const double effective_chemical_potential = gsl_vector_get(roots_array, 0);

  gsl_vector_set(function, 0,
                 root_equation_effective_chemical_potential(
                     degeneracy, mass, number_density, temperature,
                     effective_chemical_potential, statistics, integrator));

  return GSL_SUCCESS;
}

void ChemicalPotentialSolver::print_state_effective_chemical_potential(
    unsigned int iter, gsl_multiroot_fsolver *solver) {
  printf(
      "\n***\nfind_effective_chemical_potential(): iter = %3u \t"
      "x = % .3f \t"
      "f(x) = % .3e \n",
      iter, gsl_vector_get(solver->x, 0), gsl_vector_get(solver->f, 0));
}

int ChemicalPotentialSolver::find_effective_chemical_potential(
    double degeneracy, double mass, double number_density, double temperature,
    double statistics, double mu_initial_guess, double solution_precision,
    Integrator *integrator, double *effective_chemical_potential) {
  const gsl_multiroot_fsolver_type *Solver_name;
  gsl_multiroot_fsolver *Root_finder;

  int status;
  size_t iter = 0;

  const size_t problem_dimension = 1;

  struct ParametersForEffectiveChemicalPotentialRootFinderOneSpecies
      parameters = {degeneracy,  mass,       number_density,
                    temperature, statistics, integrator};

  gsl_multiroot_function EffectiveChemicalPotential = {
      &(ChemicalPotentialSolver::
            root_equation_effective_chemical_potential_for_GSL),
      problem_dimension, &parameters};

  double roots_array_initial[problem_dimension] = {mu_initial_guess};

  gsl_vector *roots_array = gsl_vector_alloc(problem_dimension);
  gsl_vector_set(roots_array, 0, roots_array_initial[0]);

  Solver_name = gsl_multiroot_fsolver_hybrids;
  Root_finder = gsl_multiroot_fsolver_alloc(Solver_name, problem_dimension);
  gsl_multiroot_fsolver_set(Root_finder, &EffectiveChemicalPotential,
                            roots_array);

  // print_state_effective_chemical_potential (iter, Root_finder);

  do {
    iter++;

    status = gsl_multiroot_fsolver_iterate(Root_finder);

    // print_state_effective_chemical_potential (iter, Root_finder);

    /*
     * Check if the solver is stuck
     */
    if (status) {
      const double conversion_factor =
          smash::hbarc * smash::hbarc * smash::hbarc;
      logg[LogArea::Main::id].warn(
          "\n\nThe GSL solver\nfind_effective_chemical_potential\nis stuck!"
          "\n\nInput parameters:"
          "\n            degeneracy = ",
          degeneracy, "\n            mass [GeV] = ", mass,
          "\nnumber_density [fm^-3] = ", number_density / (conversion_factor),
          "\n     temperature [GeV] = ", temperature,
          "\n            statistics = ", statistics,
          "\n    solution_precision = ", solution_precision,
          "\n\n"
          "Initialization cannot continue without calculating the "
          "chemical potential. \nTry adjusting the initial "
          "guess or the solution precision.\n\n\n");
      throw std::runtime_error(
          "Chemical potential calculation returned no result.\n\n");
      continue;
    }

    status = gsl_multiroot_test_residual(Root_finder->f, solution_precision);

    if (status == GSL_SUCCESS) {
      effective_chemical_potential[0] = gsl_vector_get(Root_finder->x, 0);
    }
  } while (status == GSL_CONTINUE && iter < 100000);

  gsl_multiroot_fsolver_free(Root_finder);
  gsl_vector_free(roots_array);

  return 0;
}

/*
 * A convenience wrapper for finding the effective chemical potential for a
 * given particle species and performing sanity checks.
 */

double ChemicalPotentialSolver::effective_chemical_potential(
    double degeneracy, double mass, double number_density, double temperature,
    double statistics, double solution_precision) {
  assert(temperature > 0);  // zero temperature case is not considered here
  assert(degeneracy > 0);
  assert(mass > 0);
  /*
   * The initial guess for the GSL chemical potential finding
   * procedure; the guess has to be different for Bose and Fermi
   * particles.
   */
  double initial_guess = 0.0;
  if (statistics == -1.0) {
    /*
     * Such initial guess allows one to sample really close to the Bose
     * condensate region; at the same time, it works fine in other regions.
     */
    initial_guess = 0.999999 * mass;
  } else {
    initial_guess = 1.05 * mass;
  }

  /// Integration precision should be better than solving precision
  const double integration_precision = 0.1 * solution_precision;
  integrator_.set_precision(integration_precision, integration_precision);

  /*
   * Array for the GSL chemical potential finder.
   * Use of array is motivated by anticipating a generalization or overloading
   * of this solver to include effective mass calculation, which will lead to
   * multidimensional root finding, in which case the result is best passed as
   * an array.
   */
  double chemical_potential[1];
  chemical_potential[0] = 0.0;

  /*
   * Call the GSL effective chemical potential finding procedure.
   */
  find_effective_chemical_potential(
      degeneracy, mass, number_density, temperature, statistics, initial_guess,
      solution_precision, &integrator_, chemical_potential);

  /*
   * Sanity checks are performed for the obtained value of the chemical
   * potential .
   */

  /*
   * Check that the calculated chemical potential reproduces the input number
   * density within acceptable error.
   */
  // precision for number density in units of fm^-3
  const double precision_fm = 1e-3;
  const double conversion_factor = smash::hbarc * smash::hbarc * smash::hbarc;
  double nDensityCheck =
      density_one_species(degeneracy, mass, temperature, chemical_potential[0],
                          statistics, &integrator_);

  if (abs((number_density - nDensityCheck) / conversion_factor) >
      precision_fm) {
    logg[LogArea::Main::id].warn(
        "\n\nThe calculated chemical potential = ", chemical_potential[0],
        " [GeV] does not reproduce the input number density to within ",
        precision_fm,
        "!"
        "\nnumber_density [GeV^3] = ",
        number_density, "\n nDensityCheck [GeV^3] = ", nDensityCheck,
        "\nnumber_density [fm^-3] = ", number_density / (conversion_factor),
        "\n nDensityCheck [fm^-3] = ", nDensityCheck / (conversion_factor),
        "\n\n"
        "Input parameters:"
        "\n            degeneracy = ",
        degeneracy, "\n            mass [GeV] = ", mass,
        "\nnumber_density [fm^-3] = ", number_density / (conversion_factor),
        "\n     temperature [GeV] = ", temperature,
        "\n            statistics = ", statistics,
        "\n    solution_precision = ", solution_precision,
        "\n\n"
        "Initialization cannot continue without a correct "
        "calculation of the chemical potential."
        "\nTry adjusting the initial guess or the solution precision.\n\n\n");
    throw std::runtime_error(
        "Chemical potential calculation does not reproduce the "
        "input number density.\n\n");
  }

  /*
   * Check if the obtained value of the effective chemical potential indicates
   * that the number density is such that for bosons one has encountered the
   * Bose-Einstein condensate. We need mu < mass in order for the system NOT to
   * be in the Bose-Einstein condensate region.
   */
  if ((statistics == -1) && (chemical_potential[0] > mass)) {
    logg[LogArea::Main::id].warn(
        "\n\nThe calculated chemical potential indicates the "
        "input number density is such that Bose-Einstein "
        "condensate is encountered."
        "Input parameters:"
        "\n            degeneracy = ",
        degeneracy, "\n            mass [GeV] = ", mass,
        "\nnumber_density [fm^-3] = ", number_density / (conversion_factor),
        "\n     temperature [GeV] = ", temperature,
        "\n            statistics = ", statistics,
        "\n    solution_precision = ", solution_precision,
        "\n\n"
        "Initialization cannot proceed with a pathological "
        "distribution function leading to the Bose-Einstein "
        "condensate. \nTry increasing the temperature or "
        "decreasing the number density.\n\n\n");
    throw std::runtime_error(
        "Chemical potential calculation indicates that the Bose-Einstein "
        "condensate is produced.");
  }

  /*
   * Check if the distribution is pathological, i.e., has negative values
   * for low momenta (Bose statistics, in the Bose-Einstein condensate regime,
   *  may produce pathological distributions for which f(p) < 0 ).
   */
  /// value of small_p in GeV; value is chosen to be small (arbitrary choice)
  const double small_p = 0.005;
  const double distribution_at_small_p = juttner_distribution_func(
      small_p, mass, temperature, chemical_potential[0], statistics);
  if (distribution_at_small_p < 0.0) {
    logg[LogArea::Main::id].warn(
        "\n\nNegative values of the distribution function at "
        "small momenta indicate that Bose-Einstein "
        "condensate is encountered."
        "Input parameters:"
        "\n            degeneracy = ",
        degeneracy, "\n            mass [GeV] = ", mass,
        "\nnumber_density [fm^-3] = ", number_density / (conversion_factor),
        "\n     temperature [GeV] = ", temperature,
        "\n            statistics = ", statistics,
        "\n    solution_precision = ", solution_precision, "\n\nf(p=", small_p,
        ") = ", distribution_at_small_p,
        "\n\n"
        "Initialization cannot proceed with a pathological "
        "distribution function leading to the Bose-Einstein "
        "condensate. \nTry increasing the temperature or "
        "decreasing the number density.\n\n\n");
    throw std::runtime_error(
        "Negative values of the distribution function at low momenta "
        "indicate that the Bose-Einstein condensate is produced.");
  }

  return chemical_potential[0];
}

}  // namespace smash
