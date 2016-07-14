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

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <array>
#include <vector>

#include "constants.h"
#include "particletype.h"


namespace Smash {

// Forward declaration of HadronGasEos - it is used in EosTable
class HadronGasEos;

/**
 * A class to hold, compute and access tabulated EoS.
 */
class EosTable {
 public:
  EosTable(double de, double dnb, int n_e, int n_b);
  struct table_element {
    double p;
    double T;
    double mub;
    double mus;
  };
  void compile_table(HadronGasEos &eos,
    const std::string& eos_savefile_name = "hadgas_eos.dat");
  void get(table_element& res, double e, double nb) const;

 private:
  int index(int ie, int inb) const { return ie*n_nb_ + inb; }
  std::vector<table_element> table_;
  double de_;
  double dnb_;
  int n_e_;
  int n_nb_;
};

/**
 * \class HadronGasEos
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
  /// List particle species included in the eos
  static const ParticleTypePtrList list_eos_particles();
  /**
   * \brief Compute energy density
   * Grand-canonical Boltzmann ideal gas, consisting of all hadrons in SMASH:
   * \f[ \epsilon = \sum \frac{g_i m_i^2 T^2}{2\pi^2(\hbar c)^3}
   *               exp \left(\frac{\mu_B B_i + \mu_S S_i}{T} \right) \times
   *               \left[ 3 K_2\left( \frac{m_i}{T}\right) +
   *               \frac{m_i}{T} K_1\left( \frac{m_i}{T}\right)\right]
   * \f]
   *
   * \param T temperature [GeV]
   * \param mub baryon chemical potential [GeV]
   * \param mus strangeness chemical potential [GeV]
   *
   * \return energy density e [GeV/fm\f$^3\f$]
   */
  static double energy_density(double T, double mub, double mus);
  /**
   * \brief Compute particle number density
   * Grand-canonical Boltzmann ideal gas, consisting of all hadrons in SMASH:
   * \f[ n = \sum \frac{g_i m_i^2 T}{2\pi^2(\hbar c)^3}
   *               exp \left(\frac{\mu_B B_i + \mu_S S_i}{T} \right)
   *               K_2\left( \frac{m_i}{T}\right)
   * \f]
   *
   * \param T temperature [GeV]
   * \param mub baryon chemical potential [GeV]
   * \param mus strangeness chemical potential [GeV]
   *
   * \return particle number density n [fm\f$^{-3}\f$]
   */
  static double density(double T, double mub, double mus);
  /// Compute pressure \f$ p = n T \f$, same arguments as above
  static double pressure(double T, double mub, double mus) {
    return T * density(T, mub, mus);
  }
  /**
   * \brief Compute net baryon density
   * Grand-canonical Boltzmann ideal gas, consisting of all hadrons in SMASH:
   * \f[ n_B = \sum B_i \frac{g_i m_i^2 T}{2\pi^2(\hbar c)^3}
   *               exp \left(\frac{\mu_B B_i + \mu_S S_i}{T} \right)
   *               K_2\left( \frac{m_i}{T}\right)
   * \f]
   *
   * \param T temperature [GeV]
   * \param mub baryon chemical potential [GeV]
   * \param mus strangeness chemical potential [GeV]
   *
   * \return net baryon density density \f$n_B\f$ [fm\f$^{-3}\f$]
   */
  static double net_baryon_density(double T, double mub, double mus);
  /**
   * \brief Compute net strangeness density
   * Grand-canonical Boltzmann ideal gas, consisting of all hadrons in SMASH:
   * \f[ n_S = \sum S_i \frac{g_i m_i^2 T}{2\pi^2(\hbar c)^3}
   *               exp \left(\frac{\mu_B B_i + \mu_S S_i}{T} \right)
   *               K_2\left( \frac{m_i}{T}\right)
   * \f]
   *
   * \param T temperature [GeV]
   * \param mub baryon chemical potential [GeV]
   * \param mus strangeness chemical potential [GeV]
   *
   * \return net strangeness density density \f$n_S\f$ [fm\f$^{-3}\f$]
   */
  static double net_strange_density(double T, double mub, double mus);
  /**
   * \brief Compute partial density of one hadron sort
   * Grand-canonical Boltzmann ideal gas:
   * \f[ n =  \frac{g m^2 T}{2\pi^2(\hbar c)^3}
   *               exp \left(\frac{\mu_B B + \mu_S S}{T} \right)
   *               K_2\left( \frac{m}{T}\right)
   * \f]
   *
   * \param T temperature [GeV]
   * \param mub baryon chemical potential [GeV]
   * \param mus strangeness chemical potential [GeV]
   *
   * \return partial density of the given hadron sort \f$n\f$ [fm\f$^{-3}\f$]
   */
  static double partial_density(const ParticleType& ptype,
                                double T, double mub, double mus);
  /**
   * Compute temperature and chemical potentials given energy-,
   * net baryon- and net strangeness density.
   *
   * \param e energy density [GeV/fm\f$^3\f$]
   * \param nb net baryon density [fm\f$^{-3}\f$]
   * \param ns net strangeness density [fm\f$^{-3}\f$]
   * \param initial_approximation (T [GeV], mub [GeV], mus [GeV])
   *        to use as starting point, optional, default = (0.15, 0.5, 0.05)
   *
   * \return array of 3 values: temperature, baryon chemical potential
   *         and strange chemical potential
   */
  std::array<double, 3> solve_eos(double e, double nb, double ns,
             std::array<double, 3> initial_approximation = {0.15, 0.0, 0.0});
  /// Compute strange chemical potential, requiring that net strangeness = 0
  static double mus_net_strangeness0(double T, double mub);
  /// Get the element of eos table
  void from_table(EosTable::table_element& res, double e, double nb) const {
    eos_table_.get(res, e, nb);
  }
  static bool is_eos_particle(const ParticleType& ptype) {
    return ptype.is_hadron() && ptype.pdgcode().charmness() == 0;
  }
  bool is_tabulated() const { return tabulate_; }

 private:
  /// A structure for passing equation parameters to the gnu library
  struct rparams {
    double e;
    double nb;
    double ns;
  };
  /**
   * Compute (unnormalized) density of one hadron sort - helper function
   * used to reduce code duplication.
   */
  static double scaled_partial_density(const ParticleType& ptype,
                                       double beta, double mub, double mus);
  /// Interfaces EoS equations to be solved to gnu library
  static int set_eos_solver_equations(const gsl_vector* x,
                                      void* params, gsl_vector* f);
  /// Helpful printout, useful for debugging if gnu equation solving goes crazy
  void print_solver_state(size_t iter) const;
  /// Constant factor, that appears in front of many thermodyn. expressions
  static constexpr double prefactor_ = 0.5*M_1_PI*M_1_PI/(hbarc*hbarc*hbarc);
  /// Precision of equation solving
  static constexpr double tolerance_ = 1.e-5;
  /// Number of equations in the system of equations to be solved
  static constexpr size_t n_equations_ = 3;
  EosTable eos_table_ = EosTable(1.e-3, 1.e-3, 900, 900);
  /**
   * Variables used by gnu equation solver. They are stored here to allocate
   * and deallocate memory for them only once. It is expected that this class
   * will be used for solving the EoS many times, so multiple allocations and
   * frees are unwanted.
   */
  gsl_vector *x_;
  gsl_multiroot_fsolver *solver_;
  /// Create an EoS table or not?
  const bool tabulate_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_HADGAS_EOS_H_
