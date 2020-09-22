/*
 *
 *    Copyright (c) 2020-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_CHEMICALPOTENTIAL_H_
#define SRC_INCLUDE_SMASH_CHEMICALPOTENTIAL_H_


// needed for solving for the distribution maximum
// and the chemical potential
#include <gsl/gsl_multiroots.h>
// needed to handle GSL vectors which are in declarations of functions
#include <gsl/gsl_vector.h>

namespace smash {

/**
 * Vector number density of one particle species, obtained through integrating
 * the Juttner distribution function over the momentum space.
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
double density_integration_one_species_unit_range(
    double degeneracy, double mass, double temperature,
    double effective_chemical_potential, double statistics, double precision);

/**
 * This block is for:
 * Struct, root equations, and procedure for finding the effective chemical
 * potential for a given particle species. This chemical potential is NOT the
 * equilibrium chemical potential of a system with multiple particle species.
 * Rather, it is the effective chemical potential of a particle species 
 * calculated as if only that species is present. (In particular, this applies
 * also to particles and antiparticles of the same species.)
 * Note that the procedure is exact for a system composed of protons and
 * neutrons only, as their mass is degenerate.
 */

/**
 * Struct object that holds the parameters relevant to finding the effective
 * chemical potential of one particle species.
 */
struct ParametersForEffectiveChemicalPotentialRootFinderOneSpecies {
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
double root_equation_effective_chemical_potential(
    double degeneracy, double mass, double number_density, double temperature,
    double effective_chemical_potential, double statistics, double precision);

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
int root_equation_effective_chemical_potential_for_GSL(
    const gsl_vector* roots_array, void* parameters, gsl_vector* function);

/**
 * A GSL utility which allows for printing out the status of the solver
 * during the root finding procedure.
 * \param[in] iter variable keeping track of how many steps in the root
 *            solving procedure have been taken
 * \param[in] solver GSL solver object, which has acces to the current best
 *            estimate of the roots and the corresponding function values
 * \return message about the current state of the solver
 */
void print_state_effective_chemical_potential(unsigned int iter,
                                             gsl_multiroot_fsolver* solver);

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
int find_effective_chemical_potential(double degeneracy, double mass,
                                      double number_density, double temperature,
                                      double statistics,
                                      double mu_initial_guess,
                                      double solution_precision,
                                      double* effective_chemical_potential);

/**
 * Convenience wrapper for finding the effective chemical potential for a
 * given particle species and performing sanity checks on the result.
 * 
 * The convenience wrapper performs the search for the effective chemical
 * potential, and then performs sanity checks, such as against encountering the
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
double effective_chemical_potential(double degeneracy, double mass,
                                    double number_density, double temperature,
                                    double statistics,
                                    double solution_precision);

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_CHEMICALPOTENTIAL_H_
