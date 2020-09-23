/*
 *
 *    Copyright (c) 2020-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/quantumsampling.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_multiroots.h>
#include <cstdlib>
#include <iostream>
#include <random>

#include "smash/chemicalpotential.h"
#include "smash/constants.h"
#include "smash/distributions.h"
#include "smash/logging.h"
#include "smash/particletype.h"

namespace smash {

/*
 * Root equations and GSL procedure for finding the momentum for which the
 * maximum of a given Juttner distribution occurs. This is needed for a method
 * of sampling the distribution function in which one samples uniformly below
 * the maximum of the distribution.
 */

double QuantumSampling::p_max_root_equation(double p, double mass,
                                            double temperature,
                                            double effective_chemical_potential,
                                            double statistics) {
  const double Ekin = std::sqrt(p * p + mass * mass);
  const double term1 =
      2 * (1 + statistics * std::exp(-(Ekin - effective_chemical_potential) /
                                     temperature));
  const double term2 = (p * p) / (temperature * Ekin);

  return term1 - term2;
}

int QuantumSampling::p_max_root_equation_for_GSL(const gsl_vector *roots_array,
                                                 void *parameters,
                                                 gsl_vector *function) {
  struct ParametersForMaximumMomentumRootFinder *par =
      static_cast<struct ParametersForMaximumMomentumRootFinder *>(parameters);

  const double mass = (par->mass);
  const double temperature = (par->temperature);
  const double effective_chemical_potential =
      (par->effective_chemical_potential);
  const double statistics = (par->statistics);

  const double p_radial = gsl_vector_get(roots_array, 0);

  gsl_vector_set(function, 0,
                 p_max_root_equation(p_radial, mass, temperature,
                                     effective_chemical_potential, statistics));

  return GSL_SUCCESS;
}

void QuantumSampling::print_state_p_max(unsigned int iter,
                                        gsl_multiroot_fsolver *solver) {
  printf(
      "\n***\nfind_maximum_of_the_distribution(): iter = %3u \t"
      "x = % .3f \t"
      "f(x) = % .3e \n",
      iter, gsl_vector_get(solver->x, 0), gsl_vector_get(solver->f, 0));
}

int QuantumSampling::find_maximum_of_the_distribution(
    double mass, double temperature, double effective_chemical_potential,
    double statistics, double p_max_initial_guess, double solution_precision,
    double *p_max) {
  const gsl_multiroot_fsolver_type *Solver_name;
  gsl_multiroot_fsolver *Root_finder;

  int status;
  size_t iter = 0;
  size_t initial_guess_update = 0;

  const size_t problem_dimension = 1;

  struct ParametersForMaximumMomentumRootFinder parameters = {
      mass, temperature, effective_chemical_potential, statistics};

  gsl_multiroot_function MaximumOfDistribution = {
      &p_max_root_equation_for_GSL, problem_dimension, &parameters};

  double roots_array_initial[1] = {p_max_initial_guess};

  gsl_vector *roots_array = gsl_vector_alloc(problem_dimension);
  gsl_vector_set(roots_array, 0, roots_array_initial[0]);

  Solver_name = gsl_multiroot_fsolver_hybrids;
  Root_finder = gsl_multiroot_fsolver_alloc(Solver_name, problem_dimension);
  gsl_multiroot_fsolver_set(Root_finder, &MaximumOfDistribution, roots_array);

  // print_state_p_max (iter, Root_finder);

  do {
    iter++;

    status = gsl_multiroot_fsolver_iterate(Root_finder);

    // print_state_p_max (iter, Root_finder);

    /*
     * Check if the solver is stuck
     */
    if (status) {
      if (initial_guess_update < 100) {
        // we're updating the starting point
        initial_guess_update++;
        p_max_initial_guess = initial_guess_update * 0.05;
        roots_array_initial[0] = p_max_initial_guess;
        gsl_vector_set(roots_array, 0, roots_array_initial[0]);
        gsl_multiroot_fsolver_set(Root_finder, &MaximumOfDistribution,
                                  roots_array);
        iter = 0;
      } else {
        logg[LogArea::Distributions::id].warn(
            "\n\nThe GSL solver\nfind_maximum_of_the_distribution\nis stuck!"
            "\n\nInput parameters:"
            "\n                  mass [GeV] = ",
            mass, "\n           temperature [GeV] = ", temperature,
            "\neffective_chemical_potential = ", effective_chemical_potential,
            "\n                  statistics = ", statistics,
            "\n          solution_precision = ", solution_precision,
            "\n\n"
            "Initialization cannot sample the momenta without "
            "calculating the distribution maximum."
            "\nTry adjusting the initial guess (which is "
            "looped over in the GSL procedure) or the "
            "solution precision."
            "\nUncomment print_state_p_max to check solver progress.\n\n\n");
        throw std::runtime_error(
            "QuantumSampling::find_maximum_of_the_distribution returned "
            "no result.\n\n");
        continue;
      }
    }

    status = gsl_multiroot_test_residual(Root_finder->f, solution_precision);

    if (status == GSL_SUCCESS) {
      p_max[0] = gsl_vector_get(Root_finder->x, 0);
    }

  } while (status == GSL_CONTINUE && iter < 100000);

  gsl_multiroot_fsolver_free(Root_finder);
  gsl_vector_free(roots_array);

  return 0;
}

double QuantumSampling::maximum_of_the_distribution(
    double mass, double temperature, double effective_chemical_potential,
    double statistics, double solution_precision) {
  /*
   * Momentum at which the distribution function has its maximum.
   */
  double p_max[1];
  p_max[0] = 0.0;

  /*
   * Initial guess for the value of p_max, in GeV. This value is
   * looped over within the GSL solver, so that many initial guesses
   * are probed if the solution is not found.
   */
  double initial_guess_p_max = 0.050;

  /*
   * Calling the GSL distribution maximum finder
   */
  find_maximum_of_the_distribution(
      mass, temperature, effective_chemical_potential, statistics,
      initial_guess_p_max, solution_precision, p_max);

  double distribution_function_maximum =
      p_max[0] * p_max[0] *
      juttner_distribution_func(p_max[0], mass, temperature,
                                effective_chemical_potential, statistics);

  return distribution_function_maximum;
}

/*
 * Initializing the QuantumSampling object triggers calculation of the
 * chemical potential and distribution maximum for all species present.
 */
QuantumSampling::QuantumSampling(
    const std::map<PdgCode, int> &initial_multiplicities, double volume,
    double temperature)
    : volume_(volume), temperature_(temperature) {
  /*
   * This is the precision which we expect from the solution; note that
   * solution precision also goes into the precision of calculating the
   * integrals involved etc. Recommended precision is at least 1e-7.
   */
  constexpr double solution_precision = 1e-8;

  for (const auto &pdg_and_mult : initial_multiplicities) {
    const PdgCode pdg = pdg_and_mult.first;
    const int number_of_particles = pdg_and_mult.second;
    const double V_in_GeV = volume_ / (hbarc * hbarc * hbarc);
    const double number_density = number_of_particles / V_in_GeV;
    const double spin_degeneracy = pdg.spin_degeneracy();
    // '+' for fermions, '-' for bosons
    const double quantum_statistics = (pdg.spin() % 2 == 0) ? -1.0 : 1.0;
    const ParticleType &ptype = ParticleType::find(pdg);
    const double particle_mass = ptype.mass();
    double chemical_potential = 0.0;
    ChemicalPotentialSolver mu_solver;
    // Calling the wrapper for the GSL chemical potential finder
    chemical_potential = mu_solver.effective_chemical_potential(
        spin_degeneracy, particle_mass, number_density, temperature_,
        quantum_statistics, solution_precision);
    effective_chemical_potentials_[pdg] = chemical_potential;
    const double distribution_function_maximum = maximum_of_the_distribution(
        particle_mass, temperature_, chemical_potential, quantum_statistics,
        solution_precision);
    distribution_function_maximums_[pdg] = distribution_function_maximum;
  }
}

/*
 * Sampling radial momenta of given particle species from Bose, Boltzmann, or
 * Fermi distribution. The choice between the distributions is made based on
 * the <statistics> variable:
 * -1 for Bose
 *  0 fot Boltzmann
 * +1 for Fermi
 *
 * This sampler is the simplest implementation of sampling based on sampling
 * from a uniform distribution.
 */
double QuantumSampling::sample(const PdgCode pdg) {
  const ParticleType &ptype = ParticleType::find(pdg);
  const double mass = ptype.mass();
  const double mu = effective_chemical_potentials_.find(pdg)->second;
  const double distr_max = distribution_function_maximums_.find(pdg)->second;
  /*
   * The variable maximum_momentum denotes the "far right" boundary of the
   * sampled region; we assume that no particle has momentum larger than 10 GeV
   */
  constexpr double maximum_momentum = 10.0;  // in [GeV]
  const double statistics = (pdg.spin() % 2 == 0) ? -1.0 : 1.0;
  double sampled_momentum = 0.0, sampled_ratio = 0.0;

  do {
    sampled_momentum = random::uniform(0.0, maximum_momentum);
    double distribution_at_sampled_p =
        sampled_momentum * sampled_momentum *
        juttner_distribution_func(sampled_momentum, mass, temperature_, mu,
                                  statistics);
    sampled_ratio = distribution_at_sampled_p / distr_max;
  } while (random::canonical() > sampled_ratio);

  return sampled_momentum;
}

}  // namespace smash
