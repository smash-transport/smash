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

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

#include <smash/integrate.h>

namespace smash {

/**
 * A class which encapsulates a GSL algorithm for finding the effective
 * chemical potential and supporting functions.
 */
class ChemicalPotentialSolver {
 public:
  /**
   * Vector number density of one particle species, obtained through integrating
   * the Juttner distribution function over the momentum space.
   * \param[in] degeneracy degeneracy g of the particle species
   * \param[in] mass (pole) mass m of the particle species [GeV]
   * \param[in] temperature temperature T of the system [GeV]
   * \param[in] effective_chemical_potential effective chemical potential of the
   *            system [GeV]
   * \param[in] statistics quantum statistics of the particles species
   *            (+1 for Fermi, -1 for Bose, 0 for Boltzmann)
   * \param[in] integrator wrapper for gsl numerical integration
   * \return the vector number density for a given species of particles [GeV^3]
   */
  static double density_one_species(double degeneracy, double mass,
                                    double temperature,
                                    double effective_chemical_potential,
                                    double statistics, Integrator* integrator);

  /**
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
  struct ParametersForChemPotSolver {
    /// degeneracy g of the particle species
    double degeneracy;

    /// (pole) mass m of the particle species
    double mass;

    /// number density n of the particle species [GeV^3]
    double number_density;

    /// temperature T of the system [GeV]
    double temperature;

    /**
     * statistics quantum statistics of the particles species
     * (+1 for Fermi, -1 for Bose, 0 for Boltzmann)
     */
    double statistics;

    /// wrapper for gsl numerical integration
    Integrator* integrator;
  };

  /**
   * Root equation for finding the value of the effective chemical potential
   * for one particle species.
   * \param[in] degeneracy degeneracy g of the particle species
   * \param[in] mass (pole) mass m of the particle species
   * \param[in] number_density number density n of the particle species [GeV^3]
   * \param[in] temperature temperature T of the system [GeV]
   * \param[in] effective_chemical_potential effective chemical potential of the
   *            system [GeV]
   * \param[in] statistics quantum statistics of the particles species
   *            (+1 for Fermi, -1 for Bose, 0 for Boltzmann)
   * \param[in] integrator wrapper for gsl numerical integration
   * \return the root equation
   */
  static double root_equation_effective_chemical_potential(
      double degeneracy, double mass, double number_density, double temperature,
      double effective_chemical_potential, double statistics,
      Integrator* integrator);

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
  static int root_equation_effective_chemical_potential_for_GSL(
      const gsl_vector* roots_array, void* parameters, gsl_vector* function);

  /**
   * A GSL utility which allows for printing out the status of the solver
   * during the root finding procedure.
   * \param[in] iter variable keeping track of how many steps in the root
   *            solving procedure have been taken
   * \param[in] solver GSL solver object, which has access to the current best
   *            estimate of the roots and the corresponding function values
   * \return message about the current state of the solver
   */
  static void print_state_effective_chemical_potential(
      unsigned int iter, gsl_multiroot_fsolver* solver);

  /**
   * A GSL root solver for finding the effective chemical potential. In practice
   * one should use the convenience wrapper effective_chemical_potential(),
   * which also performs sanity checks on the obtained solution.
   *
   * Solver of the equation for chemical potential \f$ \mu \f$
   * \f[ n = \frac{g}{2 \pi^2} \int
   *          \frac{p^2 dp}{ e^{(\sqrt{p^2 + m^2} - \mu)/T} \pm 1}
   * \f] given the density \f$ n \f$, and temperature \f$ T \f$.
   * \param[in] degeneracy degeneracy g of the particle species
   * \param[in] mass m (pole) mass m of the particle species [GeV]
   * \param[in] number_density number density n of the particle species n
   *            [GeV^3]
   * \param[in] temperature temperature T of the system in GeV
   * \param[in] statistics quantum statistics of the particles species
   *            (+1 for Fermi, -1 for Bose, 0 for Boltzmann)
   * \param[in] mu_initial_guess the initial guess for the value of the
   *            effective chemical potential [GeV]
   * \param[in] solution_precision precision with which the solution is found;
   *            note that this precision also goes into the precision of
   *            integrals calculated in the process
   * \param[in] integrator wrapper for gsl numerical integration
   * \param[out] effective_chemical_potential the solution stored in an array
   *             object (we use an array in anticipation of generalizing to
   *             multidimensional root finding, needed for example when scalar
   *             interactions are present and effective mass is calculated at
   *             the same time as the effective chemical potential)
   */
  static int find_effective_chemical_potential(
      double degeneracy, double mass, double number_density, double temperature,
      double statistics, double mu_initial_guess, double solution_precision,
      Integrator* integrator, double* effective_chemical_potential);

  /**
   * Convenience wrapper for finding the effective chemical potential for a
   * given particle species and performing sanity checks on the result.
   *
   * \param[in] degeneracy degeneracy g of the particle species
   * \param[in] mass (pole) mass m of the particle species [GeV]
   * \param[in] number_density number density n of the particle species [GeV^3]
   * \param[in] temperature temperature T of the system [GeV]
   * \param[in] statistics quantum statistics of the particles species
   *            (+1 for Fermi, -1 for Bose, 0 for Boltzmann)
   * \param[in] solution_precision precision with which the solution is found;
   *            note that this precision also goes into the precision of
   *            integrals calculated in the process
   * \return effective chemical potential mu, the solution [GeV]
   */
  double effective_chemical_potential(double degeneracy, double mass,
                                      double number_density, double temperature,
                                      double statistics,
                                      double solution_precision);

 private:
  /**
   * A wrapper for gsl numerical integration
   */
  Integrator integrator_ = Integrator(1E6);
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_CHEMICALPOTENTIAL_H_
