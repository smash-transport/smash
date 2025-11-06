/*
 *
 *    Copyright (c) 2016-2020,2022,2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#ifndef SRC_INCLUDE_SMASH_HADGAS_EOS_H_
#define SRC_INCLUDE_SMASH_HADGAS_EOS_H_

#include <array>
#include <string>
#include <vector>

#include "gsl/gsl_multiroots.h"
#include "gsl/gsl_roots.h"
#include "gsl/gsl_vector.h"

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
   * Sets up a table p/T/muB/mus/muQ versus (e, nb, nq), where e - energy
   * density, nb - net baryon density, nq - net charge density, p - pressure, T
   * - temperature, muB - net baryon chemical potential, muS - net strangeness
   * potential, muQ - net charge chemical potential. Net strangeness density and
   * isospin projection density are assumed to be 0 (Note that the corresponding
   * chemical potential is still non-zero, because muB != 0).
   *
   * After calling this constructor the table is allocated, but it is
   * still empty. To compute values call compile_table.
   *
   * \param[in] de step in energy density [GeV/fm^4]
   * \param[in] dnb step in net baryon density [GeV/fm^3]
   * \param[in] dq step in net charge density [GeV/Gev^3]
   * \param[in] n_e number of steps in energy density
   * \param[in] n_b number of steps in net baryon density
   * \param[in] n_q number of steps in net charge density
   *
   * Entry at (ie, inb, inq) corresponds to energy density and net baryon
   * density (e, nb, nq) = (ie*de, inb*dnb, inq*dnq) [GeV/fm^4, GeV/fm^3].
   */
  EosTable(double de, double dnb, double dq, size_t n_e, size_t n_b,
           size_t n_q);
  /// Define the data structure for one element of the table.
  struct table_element {
    /// Pressure
    double p;
    /// Temperature
    double T;
    /// Net baryochemical potential
    double mub;
    /// Net strangeness potential
    double mus;
    /// Net charge chemical potential
    double muq;
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
   * Obtain interpolated p/T/muB/muS/muQ from the tabulated equation of state
   * given energy density, net baryon density and net charge density
   *
   * \param[in] e energy density
   * \param[in] nb net baryon density
   * \param[in] nq net charge density
   * \param[out] res structure, that contains p/T/muB/muS/muQ
   */
  void get(table_element& res, double e, double nb, double nq) const;

 private:
  /// proper index in a 1d vector, where the 3d table is stored
  size_t index(size_t ie, size_t inb, size_t inq) const {
    return n_q_ * (ie * n_nb_ + inb) + inq;
  }
  /// Storage for the tabulated equation of state
  std::vector<table_element> table_;
  /// Step in energy density
  double de_;
  /// Step in net-baryon density
  double dnb_;
  /// Step in net-charge density
  double dq_;
  /// Number of steps in energy density
  size_t n_e_;
  /// Number of steps in net-baryon density
  size_t n_nb_;
  /// Number of steps in net-charge density
  size_t n_q_;
};

/**
 * Class to handle the equation of state (EoS) of the hadron gas, consisting
 * of all hadrons included in SMASH. This implementation deals with an ideal
 * Boltzmann gas and allows to compute:
 *   - energy density \f$\epsilon\f$, pressure \f$p\f$, density \f$n\f$,
 *     net baryon density \f$n_B\f$, net strangeness \f$n_S\f$ and net charge
 *     density \f$n_Q\f$ as a function of temperature \f$T\f$, baryon chemical
 *     potential \f$\mu_B\f$, strange chemical potential \f$\mu_S\f$ and charge
 *     chemical potential \f$\mu_Q\f$.
 *   - Temperature and chemical potentials given energy-, net baryon- and
 *     net strangeness density. This requires solving a system of
 *     nonlinear equations.
 */
class HadronGasEos {
 public:
  /**
   *  Constructor of HadronGasEos
   *  \param[in] tabulate Whether the equation of state should be tabulated
   *             Tabulation takes time once (typically around 5 minutes), but
   *             makes the further usage of the class much faster. Tabulated
   *             values are saved in a file and loaded at the next run.
   *  \param[in] account_for_widths Whether equation of state should account
   *             for resonance spectral functions. Normally one wants to do it,
   *             if HadronGasEos is used for density calculations,
   *             for example in the box initialization. However, it is not
   *             recommended to account for spectral functions if EoS is
   *             tabulated (tabulate = true), because this makes tabulation
   *             incredibly slow. Therefore, for HadronGasEos to be used
   *             in thermalizer, this option has to be false. Also note that
   *             presently width account is not implemented for energy density
   *             calculation.
   */
  HadronGasEos(bool tabulate, bool account_for_widths);
  ~HadronGasEos();

  /**
   * \brief Compute energy density.
   *
   * Grand-canonical Boltzmann ideal gas, consisting of all hadrons in SMASH:
   * \f[ \epsilon = \sum \frac{g_i m_i^2 T^2}{2\pi^2(\hbar c)^3}
   *               exp \left(\frac{\mu_B B_i + \mu_S S_i + \mu_Q Q_i}{T} \right)
   * \times \left[ 3 K_2\left( \frac{m_i}{T}\right) + \frac{m_i}{T} K_1\left(
   * \frac{m_i}{T}\right)\right] \f]
   *
   * \param[in] T temperature [GeV]
   * \param[in] mub baryon chemical potential [GeV]
   * \param[in] mus strangeness chemical potential [GeV]
   * \param[in] muq charge chemical potential [GeV]
   * \return energy density e [GeV/fm\f$^3\f$]
   */
  static double energy_density(double T, double mub, double mus, double muq);

  /**
   * \brief Compute particle number density.
   *
   * Grand-canonical Boltzmann ideal gas, consisting of all hadrons in SMASH:
   * \f[ n = \sum \frac{g_i m_i^2 T}{2\pi^2(\hbar c)^3}
   *               exp \left(\frac{\mu_B B_i + \mu_S S_i + \mu_Q Q_i}{T} \right)
   *               K_2\left( \frac{m_i}{T}\right)
   * \f]
   *
   * \param[in] T temperature [GeV]
   * \param[in] mub baryon chemical potential [GeV]
   * \param[in] mus strangeness chemical potential [GeV]
   * \param[in] muq charge chemical potential [GeV]
   * \param[in] account_for_resonance_widths if false, pole masses are used;
   *            if true, then integration over spectral function is included
   * \return particle number density n [fm\f$^{-3}\f$]
   */
  static double density(double T, double mub, double mus, double muq,
                        bool account_for_resonance_widths = false);

  /**
   * Compute pressure \f$ p = n T \f$.
   *
   * \param[in] T temperature [GeV]
   * \param[in] mub baryon chemical potential [GeV]
   * \param[in] mus strangeness chemical potential [GeV]
   * \param[in] muq charge chemical potential [GeV]
   * \param[in] account_for_resonance_widths if false, pole masses are used;
   *            if true, then integration over spectral function is included
   * \return pressure p [GeV/fm\f$^{-3}\f$]
   */
  static double pressure(double T, double mub, double mus, double muq,
                         bool account_for_resonance_widths = false) {
    return T * density(T, mub, mus, muq, account_for_resonance_widths);
  }

  /**
   * \brief Compute net baryon density.
   *
   * Grand-canonical Boltzmann ideal gas, consisting of all hadrons in SMASH:
   * \f[ n_B = \sum B_i \frac{g_i m_i^2 T}{2\pi^2(\hbar c)^3}
   *               exp \left(\frac{\mu_B B_i + \mu_S S_i + \mu_Q Q_i}{T} \right)
   *               K_2\left( \frac{m_i}{T}\right)
   * \f]
   *
   * \param[in] T temperature [GeV]
   * \param[in] mub baryon chemical potential [GeV]
   * \param[in] mus strangeness chemical potential [GeV]
   * \param[in] muq charge chemical potential [GeV]
   * \param[in] account_for_resonance_widths if false, pole masses are used;
   *            if true, then integration over spectral function is included
   * \return net baryon density \f$n_B\f$ [fm\f$^{-3}\f$]
   */
  static double net_baryon_density(double T, double mub, double mus, double muq,
                                   bool account_for_resonance_widths = false);

  /**
   * \brief Compute net strangeness density.
   *
   * Grand-canonical Boltzmann ideal gas, consisting of all hadrons in SMASH:
   * \f[ n_S = \sum S_i \frac{g_i m_i^2 T}{2\pi^2(\hbar c)^3}
   *               exp \left(\frac{\mu_B B_i + \mu_S S_i + \mu_Q Q_i}{T} \right)
   *               K_2\left( \frac{m_i}{T}\right)
   * \f]
   *
   * \param[in] T temperature [GeV]
   * \param[in] mub baryon chemical potential [GeV]
   * \param[in] mus strangeness chemical potential [GeV]
   * \param[in] muq charge chemical potential [GeV]
   * \param[in] account_for_resonance_widths if false, pole masses are used;
   *            if true, then integration over spectral function is included
   * \return net strangeness density density \f$n_S\f$ [fm\f$^{-3}\f$]
   */
  static double net_strange_density(double T, double mub, double mus,
                                    double muq,
                                    bool account_for_resonance_widths = false);

  /**
   * \brief Compute net charge density.
   *
   * Grand-canonical Boltzmann ideal gas, consisting of all hadrons in SMASH:
   * \f[ n_Q = \sum Q_i \frac{g_i m_i^2 T}{2\pi^2(\hbar c)^3}
   *               exp \left(\frac{\mu_B B_i + \mu_S S_i + \mu_Q Q_i}{T} \right)
   *               K_2\left( \frac{m_i}{T}\right)
   * \f]
   *
   * \param[in] T temperature [GeV]
   * \param[in] mub baryon chemical potential [GeV]
   * \param[in] mus strangeness chemical potential [GeV]
   * \param[in] muq charge chemical potential [GeV]
   * \param[in] account_for_resonance_widths if false, pole masses are used;
   *            if true, then integration over spectral function is included
   * \return net charge density density \f$n_S\f$ [fm\f$^{-3}\f$]
   */
  static double net_charge_density(double T, double mub, double mus, double muq,
                                   bool account_for_resonance_widths = false);

  /**
   * \brief Compute partial density of one hadron sort.
   *
   * Grand-canonical Boltzmann ideal gas:
   * \f[ n =  \frac{g m^2 T}{2\pi^2(\hbar c)^3}
   *               exp \left(\frac{\mu_B B + \mu_S S + \mu_Q Q_i}{T} \right)
   *               K_2\left( \frac{m}{T}\right)
   * \f]
   *
   * \param[in] ptype the hadron sort, for which partial density is computed
   * \param[in] T temperature [GeV]
   * \param[in] mub baryon chemical potential [GeV]
   * \param[in] mus strangeness chemical potential [GeV]
   * \param[in] muq charge chemical potential [GeV]
   * \param[in] account_for_resonance_widths if false, pole masses are used;
   *            if true, then integration over spectral function is included
   * \return partial density of the given hadron sort \f$n\f$ [fm\f$^{-3}\f$]
   */
  static double partial_density(const ParticleType& ptype, double T, double mub,
                                double mus, double muq,
                                bool account_for_resonance_widths = false);
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
  static double sample_mass_thermal(const ParticleType& ptype, double beta);
  /**
   * Compute temperature and chemical potentials given energy-,
   * net baryon-, net strangeness- and net charge density and an
   * inital approximation.
   *
   * \param[in] e energy density [GeV/fm\f$^3\f$]
   * \param[in] nb net baryon density [fm\f$^{-3}\f$]
   * \param[in] ns net strangeness density [fm\f$^{-3}\f$]
   * \param[in] nq net charge density [fm\f$^{-3}\f$]
   * \param[in] initial_approximation (T [GeV], mub [GeV], mus [GeV])
   *        to use as starting point
   * \return array of 4 values: temperature, baryon chemical potential,
   *          strange chemical potential and charge chemical potential
   */
  std::array<double, 4> solve_eos(double e, double nb, double ns, double nq,
                                  std::array<double, 4> initial_approximation);

  /**
   * Compute temperature and chemical potentials given energy-,
   * net baryon-, net strangeness- and net charge density without an
   * inital approximation.
   *
   * \param[in] e energy density [GeV/fm\f$^3\f$]
   * \param[in] nb net baryon density [fm\f$^{-3}\f$]
   * \param[in] ns net strangeness density [fm\f$^{-3}\f$]
   * \param[in] nq net charge density [fm\f$^{-3}\f$]
   * \return array of 4 values: temperature, baryon chemical potential
   *         and strange chemical potential and charge
   */
  std::array<double, 4> solve_eos(double e, double nb, double ns, double nq) {
    return solve_eos(e, nb, ns, nq, solve_eos_initial_approximation(e, nb, nq));
  }

  /**
   * Compute a reasonable initial approximation for solve_eos.
   *
   * \param[in] e energy density [GeV/fm\f$^3\f$]
   * \param[in] nb net baryon density [fm\f$^{-3}\f$]
   * \param[in] nq net charge density [fm\f$^{-3}\f$]
   * \return array of 3 values: temperature, baryon chemical potential
   *         and strange chemical potential
   */
  std::array<double, 4> solve_eos_initial_approximation(double e, double nb,
                                                        double nq);

  /**
   * Compute strangeness chemical potential, requiring that net strangeness = 0
   *
   * \param[in] T temperature [GeV]
   * \param[in] mub baryon chemical potential [GeV]
   * \param[in] muq charge chemical potential [GeV]
   * \return strangeness chemical potential [GeV]
   */
  static double mus_net_strangeness0(double T, double mub, double muq);

  /// Get the element of eos table
  void from_table(EosTable::table_element& res, double e, double nb,
                  double nq) const {
    eos_table_.get(res, e, nb, nq);
  }

  /// Check if a particle belongs to the EoS
  static bool is_eos_particle(const ParticleType& ptype) {
    return ptype.is_hadron() && !ptype.pdgcode().is_heavy_flavor();
  }

  /// Create an EoS table or not?
  bool is_tabulated() const { return tabulate_; }

  /// If resonance spectral functions are taken into account
  bool account_for_resonance_widths() const {
    return account_for_resonance_widths_;
  }

 private:
  /// A structure for passing equation parameters to the gnu library
  struct rparams {
    /// energy density
    double e;
    /// net baryon density
    double nb;
    /// net strange density
    double ns;
    /// net charge density
    double nq;
    /// use pole masses of resonances, or integrate over spectral functions
    bool account_for_width;
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
   * \param[in] muq charge chemical potential [GeV]
   * \param[in] account_for_width Take hadron spectral functions into account
   *            or not. When taken into account, they result in a considerable
   *            slow down.
   * \return partial (unnormalized) density of the given hadron sort
   *         \f$n\f$ [fm\f$^{-3}\f$]
   */
  static double scaled_partial_density(const ParticleType& ptype, double beta,
                                       double mub, double mus, double muq,
                                       bool account_for_width = false);

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
  static constexpr size_t n_equations_ = 4;

  /// EOS Table to be used
  EosTable eos_table_ = EosTable(1.e-1, 1.e-1, 1.e-1, 90, 90, 90);

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

  /// Use pole masses of resonances or integrate over spectral functions
  const bool account_for_resonance_widths_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_HADGAS_EOS_H_
