/*
 *
 *    Copyright (c) 2016-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#ifndef SRC_INCLUDE_HADGAS_EOS_H_
#define SRC_INCLUDE_HADGAS_EOS_H_

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>

#include <array>
#include <string>
#include <vector>

#include "constants.h"
#include "particletype.h"

namespace smash {

// Forward declaration of HadronGasEos - it is used in EosTable
class HadronGasEos;

/**
 * A class to hold, compute and access tabulated EoS.
 */

class EosTable {
 public:
  /**
   * Sets up a table p/T/muB/mus versus (e, nb), where e is energy density,
   * nb is net baryon density, p - pressure, T - temperature, muB -
   * net baryon chemical potential, muS - net strangeness potential.
   * Net strangeness density and isospin projection density are assumed to be 0
   * (Note that corresponding chemical potentials are still non-zero,
   * because muB != 0).
   *
   * After calling this constructor the table is allocated, but it is
   * still empty. To compute values call compile_table.
   *
   * \param[in] de step in energy density [GeV/fm^4]
   * \param[in] dnb step in net baryon density [GeV/fm^3]
   * \param[in] n_e number of steps in energy density
   * \param[in] n_b number of steps in net baryon density
   *
   * Entry at (ie, inb) corresponds to energy density and net baryon density
   * (e, nb) = (ie*de, inb*dnb) [GeV/fm^4, GeV/fm^3].
   */
  EosTable(double de, double dnb, size_t n_e, size_t n_b);
  struct table_element {
    double p;
    double T;
    double mub;
    double mus;
  };
  /**
   * Computes the actual content of the table (for EosTable description see
   * documentation of the constructor).
   *
   * \param[in] eos equation of state
   * \param[in] eos_savefile_name name of the file to save tabulated equation
   *            of state
   */
  void compile_table(HadronGasEos& eos,
                     const std::string& eos_savefile_name = "hadgas_eos.dat");
  /**
   * Obtain interpolated p/T/muB/muS from the tabulated equation of state
   * given energy density and net baryon density.
   *
   * \param[in] e energy density
   * \param[in] nb net baryon density
   * \param[out] res structure, that contains p/T/muB/muS
   */
  void get(table_element& res, double e, double nb) const;

 private:
  /// proper index in a 1d vector, where the 2d table is stored
  size_t index(size_t ie, size_t inb) const { return ie * n_nb_ + inb; }
  /// Storage for the tabulated equation of state
  std::vector<table_element> table_;
  double de_;
  double dnb_;
  size_t n_e_;
  size_t n_nb_;
};

/**
 * Class to handle the equation of state (EoS) of the hadron gas, consisting
 * of all hadrons included into SMASH. This implementation deals with ideal
 * Boltzmann gas and allows to compute:
 *   - energy density \f$\epsilon\f$, pressure \f$p\f$, density \f$n\f$,
 *     net baryon density \f$n_B\f$ and net strangeness \f$n_S\f$ as a
 *     function of temperature \f$T\f$, baryon chemical potential \f$\mu_B\f$
 *     and strange chemical potential \f$\mu_S\f$.
 *   - Temperature and chemical potentials given energy-, net baryon- and
 *     net strangeness density. This requires solving a system of
 *     nonlinear equations.
 */
class HadronGasEos {
 public:
  explicit HadronGasEos(const bool tabulate = false);
  ~HadronGasEos();

  /**
   * \brief Compute energy density.
   *
   * Grand-canonical Boltzmann ideal gas, consisting of all hadrons in SMASH:
   * \f[ \epsilon = \sum \frac{g_i m_i^2 T^2}{2\pi^2(\hbar c)^3}
   *               exp \left(\frac{\mu_B B_i + \mu_S S_i}{T} \right) \times
   *               \left[ 3 K_2\left( \frac{m_i}{T}\right) +
   *               \frac{m_i}{T} K_1\left( \frac{m_i}{T}\right)\right]
   * \f]
   *
   * \param[in] T temperature [GeV]
   * \param[in] mub baryon chemical potential [GeV]
   * \param[in] mus strangeness chemical potential [GeV]
   * \return energy density e [GeV/fm\f$^3\f$]
   */
  static double energy_density(double T, double mub, double mus);

  /**
   * \brief Compute particle number density.
   *
   * Grand-canonical Boltzmann ideal gas, consisting of all hadrons in SMASH:
   * \f[ n = \sum \frac{g_i m_i^2 T}{2\pi^2(\hbar c)^3}
   *               exp \left(\frac{\mu_B B_i + \mu_S S_i}{T} \right)
   *               K_2\left( \frac{m_i}{T}\right)
   * \f]
   *
   * \param[in] T temperature [GeV]
   * \param[in] mub baryon chemical potential [GeV]
   * \param[in] mus strangeness chemical potential [GeV]
   * \return particle number density n [fm\f$^{-3}\f$]
   */
  static double density(double T, double mub, double mus);

  /**
   * Compute pressure \f$ p = n T \f$.
   *
   * \param[in] T temperature [GeV]
   * \param[in] mub baryon chemical potential [GeV]
   * \param[in] mus strangeness chemical potential [GeV]
   * \return pressure p [GeV/fm\f$^{-3}\f$]
   */
  static double pressure(double T, double mub, double mus) {
    return T * density(T, mub, mus);
  }

  /**
   * \brief Compute net baryon density.
   *
   * Grand-canonical Boltzmann ideal gas, consisting of all hadrons in SMASH:
   * \f[ n_B = \sum B_i \frac{g_i m_i^2 T}{2\pi^2(\hbar c)^3}
   *               exp \left(\frac{\mu_B B_i + \mu_S S_i}{T} \right)
   *               K_2\left( \frac{m_i}{T}\right)
   * \f]
   *
   * \param[in] T temperature [GeV]
   * \param[in] mub baryon chemical potential [GeV]
   * \param[in] mus strangeness chemical potential [GeV]
   * \return net baryon density density \f$n_B\f$ [fm\f$^{-3}\f$]
   */
  static double net_baryon_density(double T, double mub, double mus);

  /**
   * \brief Compute net strangeness density.
   *
   * Grand-canonical Boltzmann ideal gas, consisting of all hadrons in SMASH:
   * \f[ n_S = \sum S_i \frac{g_i m_i^2 T}{2\pi^2(\hbar c)^3}
   *               exp \left(\frac{\mu_B B_i + \mu_S S_i}{T} \right)
   *               K_2\left( \frac{m_i}{T}\right)
   * \f]
   *
   * \param[in] T temperature [GeV]
   * \param[in] mub baryon chemical potential [GeV]
   * \param[in] mus strangeness chemical potential [GeV]
   * \return net strangeness density density \f$n_S\f$ [fm\f$^{-3}\f$]
   */
  static double net_strange_density(double T, double mub, double mus);

  /**
   * \brief Compute partial density of one hadron sort.
   *
   * Grand-canonical Boltzmann ideal gas:
   * \f[ n =  \frac{g m^2 T}{2\pi^2(\hbar c)^3}
   *               exp \left(\frac{\mu_B B + \mu_S S}{T} \right)
   *               K_2\left( \frac{m}{T}\right)
   * \f]
   *
   * \param[in] ptype the hadron sort, for which partial density is computed
   * \param[in] T temperature [GeV]
   * \param[in] mub baryon chemical potential [GeV]
   * \param[in] mus strangeness chemical potential [GeV]
   * \return partial density of the given hadron sort \f$n\f$ [fm\f$^{-3}\f$]
   */
  static double partial_density(const ParticleType& ptype, double T, double mub,
                                double mus);
  /**
   * \brief Sample resonance mass in a thermal medium
   *
   * Samples mass from the distribution
   * \f[ dN/dm \sim A(m) m^2 K_2\left( \frac{m}{T}\right) \f]
   * For stable particles always returns pole mass.
   * \param[in] ptype the hadron sort, for which mass is sampled
   * \param[in] beta inverse temperature 1/T [1/GeV]
   * \return sampled mass
   */
  static double sample_mass_thermal(const ParticleType &ptype, double beta);
  /**
   * Compute temperature and chemical potentials given energy-,
   * net baryon-, net strangeness density and an inital approximation.
   *
   * \param[in] e energy density [GeV/fm\f$^3\f$]
   * \param[in] nb net baryon density [fm\f$^{-3}\f$]
   * \param[in] ns net strangeness density [fm\f$^{-3}\f$]
   * \param[in] initial_approximation (T [GeV], mub [GeV], mus [GeV])
   *        to use as starting point
   * \return array of 3 values: temperature, baryon chemical potential
   *         and strange chemical potential
   */
  std::array<double, 3> solve_eos(double e, double nb, double ns,
                                  std::array<double, 3> initial_approximation);

  /**
   * Compute temperature and chemical potentials given energy-,
   * net baryon- and net strangeness density without an inital approximation.
   *
   * \param[in] e energy density [GeV/fm\f$^3\f$]
   * \param[in] nb net baryon density [fm\f$^{-3}\f$]
   * \param[in] ns net strangeness density [fm\f$^{-3}\f$]
   * \return array of 3 values: temperature, baryon chemical potential
   *         and strange chemical potential
   */
  std::array<double, 3> solve_eos(double e, double nb, double ns) {
    return solve_eos(e, nb, ns, solve_eos_initial_approximation(e, nb));
  }

  /**
   * Compute a reasonable initial approximation for solve_eos.
   *
   * \param[in] e energy density [GeV/fm\f$^3\f$]
   * \param[in] nb net baryon density [fm\f$^{-3}\f$]
   * \return array of 3 values: temperature, baryon chemical potential
   *         and strange chemical potential
   */
  std::array<double, 3> solve_eos_initial_approximation(double e, double nb);

  /**
   * Compute strangeness chemical potential, requiring that net strangeness = 0
   *
   * \param[in] mub baryon chemical potential [GeV]
   * \param[in] T temperature [GeV]
   * \return strangeness chemical potential [GeV]
   */
  static double mus_net_strangeness0(double T, double mub);

  /// Get the element of eos table
  void from_table(EosTable::table_element& res, double e, double nb) const {
    eos_table_.get(res, e, nb);
  }

  /// Check if a particle belongs to the EoS
  static bool is_eos_particle(const ParticleType& ptype) {
    return ptype.is_hadron() && ptype.pdgcode().charmness() == 0;
  }

  /// Create an EoS table or not?
  bool is_tabulated() const { return tabulate_; }

 private:
  /// A structure for passing equation parameters to the gnu library
  struct rparams {
    /// energy density
    double e;
    /// net baryon density
    double nb;
    /// net strange density
    double ns;
  };

  /// Another structure for passing energy density to the gnu library
  struct eparams {
    /// energy density
    double edens;
  };

  /**
   * Function used to avoid duplications in density calculations.
   * \param[in] m_over_T mass to temperature ratio \f$ m/T \f$
   * \param[in] mu_over_T chemical potential to temperature ratio \f$ \mu/T \f$
   * \return calculated \f$ (m/T)^2 exp(\mu/T) K_2(m/T) \f$
   */
  static double scaled_partial_density_auxiliary(double m_over_T,
                                                 double mu_over_T);
  /**
   * Compute (unnormalized) density of one hadron sort - helper functions
   * used to reduce code duplication.
   *
   * \param[in] ptype the hadron sort, for which partial density is computed
   * \param[in] beta inverse temperature [1/GeV]
   * \param[in] mub baryon chemical potential [GeV]
   * \param[in] mus strangeness chemical potential [GeV]
   * \return partial (unnormalized) density of the given hadron sort
   *         \f$n\f$ [fm\f$^{-3}\f$]
   */
  static double scaled_partial_density(const ParticleType& ptype, double beta,
                                       double mub, double mus);

  /// Interface EoS equations to be solved to gnu library
  static int set_eos_solver_equations(const gsl_vector* x, void* params,
                                      gsl_vector* f);

  /// \see set_eos_solver_equations()
  static double e_equation(double T, void* params);

  /**
   * Helpful printout, useful for debugging if gnu equation solving goes crazy
   *
   * \param[in] iter current value of iterator
   * \return debug output string with iter, x and f(x) from solver
   */
  std::string print_solver_state(size_t iter) const;

  /// Constant factor, that appears in front of many thermodyn. expressions
  static constexpr double prefactor_ =
      0.5 * M_1_PI * M_1_PI / (hbarc * hbarc * hbarc);

  /// Precision of equation solving
  static constexpr double tolerance_ = 1.e-8;

  /// Number of equations in the system of equations to be solved
  static constexpr size_t n_equations_ = 3;

  /// EOS Table to be used
  EosTable eos_table_ = EosTable(1.e-2, 1.e-2, 900, 900);

  /**
   * Variables used by gnu equation solver. They are stored here to allocate
   * and deallocate memory for them only once. It is expected that this class
   * will be used for solving the EoS many times, so multiple allocations and
   * frees are unwanted.
   */
  gsl_vector* x_;

  /// \see x_
  gsl_multiroot_fsolver* solver_;

  /// Create an EoS table or not?
  const bool tabulate_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_HADGAS_EOS_H_
