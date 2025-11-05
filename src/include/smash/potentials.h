/*
 *
 *    Copyright (c) 2014-2015,2017-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_POTENTIALS_H_
#define SRC_INCLUDE_SMASH_POTENTIALS_H_

#include <tuple>
#include <utility>
#include <vector>

#include "gsl/gsl_multiroots.h"
#include "gsl/gsl_vector.h"

#include "configuration.h"
#include "density.h"
#include "forwarddeclarations.h"
#include "particledata.h"
#include "rootsolver.h"
#include "threevector.h"

namespace smash {
static constexpr int LPotentials = LogArea::Potentials::id;
/**
 * A class that stores parameters of potentials, calculates
 * potentials and their gradients. Potentials are responsible
 * for long-range interactions and stand in the left part of
 * Boltzmann equation. Short-range interactions are taken into
 * account in the right part of it - in the collision term.
 */
class Potentials {
 public:
  /**
   * Potentials constructor.
   *
   * \param[in] conf Configuration which contains the switches
   *            determining whether to turn on the Skyrme or the
   *            symmetry potentials, and the coefficients controlling
   *            how strong the potentials are.
   * \param[in] parameters Struct that contains the gaussian smearing factor
   *            \f$\sigma\f$, the distance cutoff \f$r_{\rm cut}\f$ and
   *            the testparticle number needed for the density calculation.
   */
  Potentials(Configuration conf, const DensityParameters &parameters);
  /// Standard destructor
  virtual ~Potentials();

  /**
   * Calculates the gradient of the single-particle energy (including
   * potentials) in the calculation frame in MeV/fm
   *
   * \param jB_lattice Pointer to the baryon density lattice
   * \param position Position of the particle of interest in fm
   * \param momentum Momentum of the particle of interest in GeV
   * \param mass Mass of the particle of interest in GeV
   * \param plist List of all particles
   * \return ThreeVector gradient of the single particle energy in the
   * calculation frame in MeV/fm
   */
  ThreeVector single_particle_energy_gradient(DensityLattice *jB_lattice,
                                              const ThreeVector &position,
                                              const ThreeVector &momentum,
                                              double mass,
                                              ParticleList &plist) const {
    const std::array<double, 3> dr = (jB_lattice)
                                         ? jB_lattice->cell_sizes()
                                         : std::array<double, 3>{0.1, 0.1, 0.1};
    ThreeVector result, position_left, position_right;
    DensityOnLattice jmu_left, jmu_right;
    FourVector net_4current_left, net_4current_right;
    for (int i = 0; i < 3; i++) {
      position_left = position;
      position_left[i] -= dr[i];
      position_right = position;
      position_right[i] += dr[i];

      if (jB_lattice && jB_lattice->value_at(position_left, jmu_left)) {
        net_4current_left = jmu_left.jmu_net();
      } else if (use_potentials_outside_lattice_) {
        auto current = current_eckart(position_left, plist, param_,
                                      DensityType::Baryon, false, true);
        net_4current_left = std::get<1>(current);
      } else {
        return {0., 0., 0.};
      }

      if (jB_lattice && jB_lattice->value_at(position_right, jmu_right)) {
        net_4current_right = jmu_right.jmu_net();
      } else if (use_potentials_outside_lattice_) {
        auto current = current_eckart(position_right, plist, param_,
                                      DensityType::Baryon, false, true);
        net_4current_right = std::get<1>(current);
      } else {
        return {0., 0., 0.};
      }
      result[i] =
          (calculation_frame_energy(momentum, net_4current_right, mass) -
           calculation_frame_energy(momentum, net_4current_left, mass)) /
          (2 * dr[i]);
    }
    return result;
  }

  /**
   * Evaluates the single-particle energy (including the potential) of a
   * particle at a given position and momentum in the calculation frame
   *
   * \param[in] momentum Momentum of interest in GeV
   * \param[in] jmu_B Baryon current density at pos
   * \param[in] mass mass of the particle of interest
   * \return the energy of a particle in the calculation frame
   **/
  double calculation_frame_energy(const ThreeVector &momentum,
                                  const FourVector &jmu_B, double mass) const {
    std::function<double(double)> root_equation = [momentum, jmu_B, mass,
                                                   this](double energy) {
      return root_eq_potentials(energy, momentum, jmu_B, mass, skyrme_a_,
                                skyrme_b_, skyrme_tau_, mom_dependence_C_,
                                mom_dependence_Lambda_);
    };
    RootSolver1D root_solver{root_equation};
    constexpr std::size_t max_number_root_solver_iterations = 100000;
    const double initial_guess = std::sqrt(mass * mass + momentum * momentum);
    const std::array<double, 4> interval_half_widths = {0.05, 0.5, 5.0, 50.0};
    for (const double half_width : interval_half_widths) {
      const std::pair<double, double> x_range = {initial_guess - half_width,
                                                 initial_guess + half_width};
      const auto calc_frame_energy = root_solver.try_find_root(
          x_range.first, x_range.second, max_number_root_solver_iterations);
      if (calc_frame_energy) {
        return *calc_frame_energy;
      } else {
        logg[LPotentials].debug()
            << "Did not find a root for potentials in the interval ["
            << x_range.first << " GeV ," << x_range.second << " GeV].";
      }
    }
    /* If the above attempts failed, try a scan in a larger interval with a fine
     * resolution. This is a compromise between trying to be sure to find a root
     * and not introduce a too large overhead. As plan B it should be fine. */
    constexpr double half_interval_width = 100;  // GeV
    constexpr double scanning_resolution = 0.1;  // GeV
    const double upper_bound = initial_guess + half_interval_width;
    const double lower_bound = initial_guess - half_interval_width;
    std::pair<double, double> x_range = {initial_guess - scanning_resolution,
                                         initial_guess + scanning_resolution};
    logg[LPotentials].debug()
        << "Trying to find a root for potentials around " << initial_guess
        << " GeV in the range [" << lower_bound << ", " << upper_bound
        << "] GeV\nExploring interval with a resolution of "
        << scanning_resolution << " GeV, starting with x_range: ["
        << x_range.first << ", " << x_range.second << "] GeV";
    while (true) {
      const auto calc_frame_energy = root_solver.try_find_root(
          x_range.first, x_range.second, max_number_root_solver_iterations);
      if (calc_frame_energy) {
        return *calc_frame_energy;
      } else {
        const bool is_possible_to_go_right = (x_range.second < upper_bound);
        const bool is_possible_to_go_left = (x_range.first > lower_bound);
        // The following logic might be compacted, but it would be much harder
        // to read and follow, hence accept here a trivial duplication of code
        if (!is_possible_to_go_left && !is_possible_to_go_right) {
          break;
        } else if (!is_possible_to_go_left) {
          x_range.second += scanning_resolution;
        } else if (!is_possible_to_go_right) {
          x_range.first -= scanning_resolution;
        } else {
          const bool go_right = random::uniform_int(0, 1);
          if (go_right) {
            x_range.second += scanning_resolution;
          } else {
            x_range.first -= scanning_resolution;
          }
        }
        logg[LPotentials].trace() << "x_range: [" << x_range.first << ", "
                                  << x_range.second << "] GeV";
      }
    }
    logg[LPotentials].debug()
        << "Did not find any sub-range with a root scanning the interval ["
        << x_range.first << ", " << x_range.second << "] GeV";
    logg[LPotentials].error(
        "Failed to find root for momentum-dependent potentials.");
    throw std::runtime_error("Unable to continue simulation.");
  }

  /**
   * Evaluates Skyrme potential given a baryon density.
   *
   * \param[in] baryon_density Baryon density \f$\rho\f$ evaluated in the
   *            local rest frame in fm\f$^{-3}\f$.
   * \return Skyrme potential \f[U_B=10^{-3}\times\frac{\rho}{|\rho|}
   *         (A\frac{\rho}{\rho_0}+B(\frac{\rho}{\rho_0})^\tau)\f] in GeV
   */
  double skyrme_pot(const double baryon_density) const {
    return skyrme_pot(baryon_density, skyrme_a_, skyrme_b_, skyrme_tau_);
  }

  /**
   * Evaluates symmetry potential given baryon isospin density.
   *
   * \note The second term is neglected if \f$\gamma\f$ is not specified in the
   * config. \param[in] baryon_isospin_density The difference between the proton
   * and the neutron density in the local rest frame in fm\f$^{-3}\f$.
   * \param[in] baryon_density
   * \return Symmetry potential \f[U_I=2\times 10^{-3}S_{\rm sym}
   *         \frac{\rho_{I_3}}{\rho_0}
   *         + \left[12.3\left(\frac{\rho_B}{\rho_0}\right)^{2/3}
   *         + 20\left(\frac{\rho_B}{\rho_0}\right)^\gamma\right]
   *         \left(\frac{\rho_{I_3}}{\rho_B}\right)^2\f] in GeV
   */
  double symmetry_pot(const double baryon_isospin_density,
                      const double baryon_density) const;

  /**
   * Calculate the factor \f$S(\rho)\f$ in the symmetry potential.
   *
   * \param[in] baryon_density baryon density
   * \return factor S in symmetry potenial
   */
  double symmetry_S(const double baryon_density) const;

  /**
   * Evaluates the FourVector potential in the VDF model given the rest frame
   * density and the computational frame baryon current.
   *
   * \param[in] rhoB rest frame baryon density, in fm\f$^{-3}\f$
   * \param[in] jmuB_net net baryon current in the computational frame, in
   *            fm\f$^{-3}\f$
   * \return VDF potential \f[A^{\mu} = 10^{-3}\times
   *         \sum_i C_i \left(\frac{\rho}{\rho_0}\right)^{b_i - 2}
   * \frac{j^{\mu}}{\rho_0}\f] in GeV
   */
  FourVector vdf_pot(double rhoB, const FourVector jmuB_net) const;

  /**
   * Evaluates potential (Skyrme with optional Symmetry or VDF) at point r.
   * For Skyrme and Symmetry options, potential is always taken in the local
   * Eckart rest frame, but point r is in the computational frame.
   *
   * \param[in] r Arbitrary space point where potential is calculated
   * \param[in] plist List of all particles to be used in \f$j^{\mu}\f$
   *            calculation. If the distance between particle and calculation
   *            point r, \f$ |r-r_i| > r_{cut} \f$ then particle input
   *            to density will be ignored.
   * \param[in] acts_on Type of particle on which potential is going to act.
   *            It gives the charges (or more precisely, the scaling factors)
   *		of the particle moving in the potential field.
   * \return Total potential energy acting on the particle: for Skyrme and
   *         Symmetry potentials, \f[U_{\rm tot} =Q_BU_B+2I_3U_I\f] in GeV,
   *         while for the VDF potential \f[U_{\rm tot} =Q_B A^0\f] in GeV,
   *         where \f$Q_B\f$ is the baryon charge scaled by the ratio of the
   *         light (u, d) quark to the total quark number and \f$I_3\f$ is the
   *         third compnent of the isospin.
   */
  double potential(const ThreeVector &r, const ParticleList &plist,
                   const ParticleType &acts_on) const;

  /**
   * Evaluates the scaling factor of the forces acting on the particles.
   *
   * The forces are equal to the product of the scaling factor and the gradient
   * of the potential. We need these scaling factors to describe the motions of
   * the hyperons as well as the anti-particles in the potentials. For Lambda
   * and Sigma, since they carry 2 light (u or d) quarks, they are affected by
   * 2/3 of the Skyrme force. Xi carries 1 light quark, it is affected by 1/3 of
   * the Skyrme force. Omega carries no light quark, so it's not affected by the
   * Skyrme force. Anti-baryons are affected by the force as large as the force
   * acting on baryons but with an opposite direction.
   *
   * \param[in] data Type of particle on which potential is going to act.
   * \return (\f$Q_B(1-\frac{|Q_S|}{3}), Q_B\f$) where \f$Q_B\f$ is the baryon
   *         charge and \f$Q_S\f$ is the strangeness.
   */
  static std::pair<double, int> force_scale(const ParticleType &data);

  /**
   * Evaluates the electric and magnetic components of the skyrme force.
   *
   * \param[in] rhoB Eckart baryon density [fm\f$^{-3}\f$].
   * \param[in] grad_j0B Gradient of baryon density [fm\f$^{-4}\f$]. This
   *            density is evaluated in the computational frame.
   * \param[in] dvecjB_dt Time derivative of the vector baryon current density
   *            [fm\f$^{-4}\f$
   * \param[in] curl_vecjB Curl of the baryon vector current
   *            density [fm\f$^{-4}\f$
   * \return (\f$E_B, B_B\f$), where
   *         \f[
   *         E_B = - V_B^\prime(\rho^\ast)(\boldsymbol{\nabla}\rho_B
   *               + \partial_t\,\mathbf{j}_B)
   *         \f]
   *         is the electro component of Skyrme force and
   *         \f[
   *         B_B = V_B^\prime(\rho^\ast) \boldsymbol{\nabla}\times\mathbf{j}_B
   *         \f]
   *         is the magnetic component of the Skyrme force
   *         with \f$\rho^\ast\f$ being the Eckart baryon density.
   */
  std::pair<ThreeVector, ThreeVector> skyrme_force(
      const double rhoB, const ThreeVector grad_j0B,
      const ThreeVector dvecjB_dt, const ThreeVector curl_vecjB) const;

  /**
   * Evaluates the electric and magnetic components of the symmetry force.
   *
   * \param[in] rhoI3 Relative isospin 3 density.
   * \param[in] grad_j0I3 Gradient of I3/I density [fm\f$^{-4}\f$]. This
   *            density is evaluated in the computational frame.
   * \param[in] dvecjI3_dt Time derivative of the I3/I vector current density
   *            [fm\f$^{-4}\f$]
   * \param[in] curl_vecjI3 Curl of the I3/I vector current density
   *            [fm\f$^{-4}\f$]
   * \param[in] rhoB Net-baryon density in the rest frame
   * \param[in] grad_j0B  Gradient of the net-baryon density in the
   *            computational frame
   * \param[in] dvecjB_dt Time derivative of the net-baryon vector current
   *            density
   * \param[in] curl_vecjB Curl of the net-baryon vector current density
   * \return (\f$E_{I_3}, B_{I_3}\f$) [GeV/fm], where
   *         \f[
   *         \mathbf{E}
   *         = - \frac{\partial V^\ast}{\partial\rho_{I_3}^\ast}
   *             (\boldsymbol{\nabla}\rho_{I_3} + \partial_t \mathbf{j}_{I_3})
   *           - \frac{\partial V^\ast}{\partial\rho_B^\ast}
   *             (\boldsymbol{\nabla}\rho_B + \partial_t\,\mathbf{j}_B)
   *         \f]
   *         is the electrical component of symmetry force and
   *         \f[
   *         \mathbf{B}
   *         = \frac{\partial V^\ast}{\rho_{I_3}^\ast}
   *           \boldsymbol{\nabla}\times\mathbf{j}_{I_3}
   *           + \frac{\partial V^\ast}{\rho_B^\ast}
   *           \boldsymbol{\nabla}\times\mathbf{j}_B
   *         \f]
   *         is the magnetic component of the symmetry force
   *         with \f$\rho^\ast\f$ being the respective Eckart density.
   */
  std::pair<ThreeVector, ThreeVector> symmetry_force(
      const double rhoI3, const ThreeVector grad_j0I3,
      const ThreeVector dvecjI3_dt, const ThreeVector curl_vecjI3,
      const double rhoB, const ThreeVector grad_j0B,
      const ThreeVector dvecjB_dt, const ThreeVector curl_vecjB) const;

  /**
   * Integrand for calculating the electric field.
   *
   * The field is calculated via
   * \f[
   * \mathbf{E}(\mathbf{r})
   * = \int\frac{(\mathbf{r}-\mathbf{r}^\prime)\rho(\mathbf{r}^\prime)}
   *            {|\mathbf{r}-\mathbf{r}^\prime|^3}d^3r^\prime
   * \f]
   *
   * \param[in] pos position vector to be integrated over
   * \param[in] charge_density electric charge density at position pos
   * \param[in] point position where to calculate the field
   */
  static ThreeVector E_field_integrand(ThreeVector pos,
                                       DensityOnLattice &charge_density,
                                       ThreeVector point) {
    ThreeVector dr = point - pos;
    if (dr.abs() < really_small) {
      return {0., 0., 0.};
    }
    return elementary_charge * charge_density.rho() * dr /
           std::pow(dr.abs(), 3);
  }

  /**
   * Integrand for calculating the magnetic field using the Biot-Savart formula.
   *
   * \param[in] pos position vector to be integrated over
   * \param[in] charge_density electric charge density and current
   * \param[in] point position where the magnetic field will be calculated
   */
  static ThreeVector B_field_integrand(ThreeVector pos,
                                       DensityOnLattice &charge_density,
                                       ThreeVector point) {
    ThreeVector dr = point - pos;
    if (dr.abs() < really_small) {
      return {0., 0., 0.};
    }
    return elementary_charge *
           charge_density.jmu_net().threevec().cross_product(dr) /
           std::pow(dr.abs(), 3);
  }
  /**
   * Evaluates the electric and magnetic components of force in the VDF model
   * given the derivatives of the baryon current \f$j^{\mu}\f$.
   *
   * \param[in] rhoB rest frame baryon density in fm\f$^{-3}\f$
   * \param[in] drhoB_dt time derivative of the rest frame density
   * \param[in] grad_rhoB gradient of the rest frame density
   * \param[in] gradrhoB_cross_vecjB cross product of the gradient of the rest
   *            frame density and the 3-vector baryon current density
   * \param[in] j0B computational frame baryon density in fm\f$^{-3}\f$
   * \param[in] grad_j0B gradient of the computational frame baryon density
   * \param[in] vecjB 3-vector baryon current
   * \param[in] dvecjB_dt time derivative of the computational frame 3-vector
   *            baryon current
   * \param[in] curl_vecjB curl of the 3-vector baryon current
   * \return (\f$E_{VDF},
   *         B_{VDF}\f$) [GeV/fm],
   *         where
   *         \f[
   *         \mathbf{E}_{VDF}
   *         = - F_1  \big[(\boldsymbol{\nabla} \rho)\,j^0 +
   *                       (\partial_t \rho)\,\mathbf{j}\big]
   *           - F_2  (\boldsymbol{\nabla} j^0 + \partial_t\,\mathbf{j})
   *         \f]
   *         is the electrical component of VDF force and
   *         \f[
   *         \mathbf{B}_{VDF}
   *         =   F_1 (\boldsymbol{\nabla} \rho) \times \mathbf{j}
   *           + F_2  \boldsymbol{\nabla} \times \mathbf{j}
   *         \f]
   *         is the magnetic component of the VDF force, with
   *         \f{aligned}
   *         F_1 &= \sum_i C_i (b_i - 2)
   *                       \frac{\rho^{b_i - 3}}{\rho_0^{b_i - 1}}\\
   *         F_2 &= \sum_i C_i \frac{\rho^{b_i - 2}}{\rho_0^{b_i - 1}} \;,
   *         \f}
   *         where \f$\rho_0\f$ is the saturation density.
   */
  std::pair<ThreeVector, ThreeVector> vdf_force(
      double rhoB, const double drhoB_dt, const ThreeVector grad_rhoB,
      const ThreeVector gradrhoB_cross_vecjB, const double j0B,
      const ThreeVector grad_j0B, const ThreeVector vecjB,
      const ThreeVector dvecjB_dt, const ThreeVector curl_vecjB) const;

  /**
   * Evaluates the electric and magnetic components of force in the VDF force
   * given the derivatives of the VDF mean-field \f$A^\mu\f$.
   *
   * \param[in] grad_A_0 gradient of the zeroth component of the field A^mu
   * \param[in] dA_dt time derivative of the field A^mu
   * \param[in] curl_vecA curl of the vector component of the field A^mu
   * \return (\f$E_{VDF}, B_{VDF}\f$) [GeV/fm], where
   *         \f[
   *         \mathbf{E}_{VDF} = - \boldsymbol{\nabla} A^0 - \partial_t\mathbf{A}
   *         \f]
   *         is the electrical component of VDF force and
   *         \f[
   *         \mathbf{B}_{VDF} = \boldsymbol{\nabla} \times \mathbf{A}
   *         \f]
   *         is the magnetic component of the VDF force.
   */
  std::pair<ThreeVector, ThreeVector> vdf_force(
      const ThreeVector grad_A_0, const ThreeVector dA_dt,
      const ThreeVector curl_vecA) const;

  /**
   * Evaluates the electric and magnetic components of the forces at point r.
   * Point r is in the computational frame.
   *
   * \param[in] r Arbitrary space point where potential gradient is calculated
   * \param[in] plist List of all particles to be used in \f$j^{\mu}\f$
   *            calculation. If the distance between particle and calculation
   *            point r, \f$ |r-r_i| > r_{cut} \f$ then particle input
   *            to density will be ignored.
   * \return (\f$E_B, B_B, E_{I_3}, B_{I_3}\f$) [GeV/fm], where
   *          \f$E_B\f$: the electric component of the Skyrme or VDF force,
   *          \f$B_B\f$: the magnetic component of the Skyrme or VDF force,
   *          \f$E_{I_3}\f$: the electric component of the symmetry force,
   *          \f$B_{I_3}\f$: the magnetic component of the symmetry force
   */
  virtual std::tuple<ThreeVector, ThreeVector, ThreeVector, ThreeVector>
  all_forces(const ThreeVector &r, const ParticleList &plist) const;

  /// \return Is Skyrme potential on?
  virtual bool use_skyrme() const { return use_skyrme_; }
  /// \return Is symmetry potential on?
  virtual bool use_symmetry() const { return use_symmetry_; }
  /// \return Is Coulomb potential on?
  virtual bool use_coulomb() const { return use_coulomb_; }
  /// \return Use momentum-dependent part of the potential?
  virtual bool use_momentum_dependence() const {
    return use_momentum_dependence_;
  }

  /// \return Skyrme parameter skyrme_a, in MeV
  double skyrme_a() const { return skyrme_a_; }
  /// \return Skyrme parameter skyrme_b, in MeV
  double skyrme_b() const { return skyrme_b_; }
  /// \return Skyrme parameter skyrme_tau
  double skyrme_tau() const { return skyrme_tau_; }
  /// \return Skyrme parameter S_pot, in MeV
  double symmetry_S_pot() const { return symmetry_S_Pot_; }

  /// \return Is VDF potential on?
  virtual bool use_vdf() const { return use_vdf_; }
  /// \return Value of the saturation density used in the VDF potential
  double saturation_density() const { return saturation_density_; }
  /// \return Vector of the VDF coefficients \f$C_i\f$, coefficients_
  const std::vector<double> &coeffs() const { return coeffs_; }
  /// \return Vector of the VDF exponents \f$b_i\f$, powers_
  const std::vector<double> &powers() const { return powers_; }
  /// \return Number of terms in the VDF potential
  int number_of_terms() const { return powers_.size(); }

  /// \return cutoff radius in ntegration for coulomb potential in fm
  double coulomb_r_cut() const { return coulomb_r_cut_; }

  /**
   * \return Wether to take potentials into account for particles outside
   * of the lattice
   */
  bool use_potentials_outside_lattice() const {
    return use_potentials_outside_lattice_;
  }

 private:
  /**
   * Struct that contains the gaussian smearing width \f$\sigma\f$,
   * the distance cutoff \f$r_{\rm cut}\f$ and the testparticle number
   * needed for the density calculation.
   */
  const DensityParameters param_;

  /// Skyrme potential on/off
  bool use_skyrme_;

  /// Symmetry potential on/off
  bool use_symmetry_;

  /// Coulomb potential on/Off
  bool use_coulomb_;

  /// VDF potential on/off
  bool use_vdf_;

  /// Momentum-dependent part on/off
  bool use_momentum_dependence_;

  /**
   * Parameter of skyrme potentials:
   * the coefficient in front of \f$\frac{\rho}{\rho_0}\f$ in GeV
   */
  double skyrme_a_;

  /**
   * Parameters of skyrme potentials:
   * the coefficient in front of \f$(\frac{\rho}{\rho_0})^\tau\f$ in GeV
   */
  double skyrme_b_;

  /**
   * Parameters of skyrme potentials:
   * the power index.
   */
  double skyrme_tau_;

  /**
   * Parameter Lambda of the momentum-dependent part of the potentials
   * given in 1/fm
   */
  double mom_dependence_Lambda_;

  /**
   * Parameter C of the momentum-dependent part of the potentials
   * given in MeV
   */
  double mom_dependence_C_;

  /// Parameter S_Pot in the symmetry potential in MeV
  double symmetry_S_Pot_;

  /**
   * Whether the baryon density dependence of the symmetry potential is
   * included
   */
  bool symmetry_is_rhoB_dependent_ = false;
  /**
   * Power \f$ \gamma \f$ in formula for \f$ S(\rho) \f$:
   * \f[ S(\rho)=12.3\,\mathrm{MeV}\times
   * \left(\frac{\rho}{\rho_0}\right)^{2/3}+20\,\mathrm{MeV}\times
   * \left(\frac{\rho}{\rho_0}\right)^\gamma \f]
   */
  double symmetry_gamma_;

  /// Cutoff in integration for coulomb potential
  double coulomb_r_cut_;

  /// Wether potentials should be included outside of the lattice
  bool use_potentials_outside_lattice_;

  /**
   * Saturation density of nuclear matter used in the VDF potential; it may
   * vary between different parameterizations.
   */
  double saturation_density_;
  /// Parameters of the VDF potential: coefficients \f$C_i\f$, in GeV
  std::vector<double> coeffs_;
  /// Parameters of the VDF potential: exponents \f$b_i\f$
  std::vector<double> powers_;

  /**
   * Calculate the derivative of the symmetry potential with respect to
   * the isospin density in GeV * fm^3
   * \f[ \frac{\partial V_\mathrm{sym}}{\partial \rho_{I_3}}
   * = 2\frac{S_\mathrm{Pot}}{\rho_0}
   * + \frac{2\rho_{I_3}\left[12.3\left(\frac{\rho_B}{\rho_0}\right)^{2/3}
   * + 20 \left(\frac{\rho_B}{\rho_0}\right)^\gamma\right]}{\rho_B^2} \f]
   *
   * \note The isospin 3 density here is actually the density of I3 / I.
   *
   * \param[in] rhoB net baryon density
   * \param[in] rhoI3 isospin density
   * \return partial derivative of the symmetry potenital with respect to the
   * isospin density.
   */
  double dVsym_drhoI3(const double rhoB, const double rhoI3) const;

  /**
   * Calculate the derivative of the symmetry potential with respect to the
   * net baryon density in GeV * fm^3
   * \f[ \frac{\partial V_\mathrm{sym}}{\partial \rho_B} =
   * \left(\frac{\rho_{I_3}}{\rho_B}\right)^2
   * \left[\frac{8.2}{\rho_0}\left(\frac{\rho_B}{\rho_0}\right)^{-1/3}
   * + \frac{20\gamma}{\rho_B}\left(\frac{\rho_B}{\rho_0}\right)^\gamma\right]
   * -2\frac{\rho_{I_3}^2}{\rho_B^3}
   * \left[12.3\left(\frac{\rho_B}{\rho_0}\right)^{2/3}
   * + 20\left(\frac{\rho_B}{\rho_0}\right)^\gamma\right]\f]
   *
   * \note The isospin 3 density here is actually the density of I3 / I
   *
   * \param[in] rhoB net baryon density
   * \param[in] rhoI3 isospin density
   * \return partial derivative of the symmetry potenital with respect to the
   *         net baryon density.
   */
  double dVsym_drhoB(const double rhoB, const double rhoI3) const;

  /**
   * Single particle Skyrme potential in MeV
   *
   * \param baryon_density net baryon density in the local rest-frame in 1/fm^3
   * \param A Skyrme parameter A in MeV
   * \param B Skyrme parameter B in MeV
   * \param tau Skyrme parameter tau
   * \return Single particle Skyrme potential in MeV
   */
  static double skyrme_pot(const double baryon_density, const double A,
                           const double B, const double tau);

  /**
   * Root equation used to determine the energy in the calculation frame
   *
   * It is the difference between energy (including potential) squared minus
   * momentum squared in the calculation frame and in the local rest-frame. The
   * equation should be zero due to Lorentz invariance but a root finder is
   * required to determine the calculation frame energy such that this is indeed
   * the case.
   *
   * \param[in] energy_calc Energy in the calculation frame at which the
   *                        equation should be evaluated in GeV
   * \param[in] momentum_calc Momentum of
   *                          the particle in calculation frame of interest
   *                          in GeV
   * \param[in] jmu Baryon current fourvector at the position of the particle
   * \param[in] m Mass of the particle of interest in GeV
   * \param[in] A Skyrme parameter A in MeV
   * \param[in] B Skyrme parameter B in MeV
   * \param[in] tau Skyrme parameter tau
   * \param[in] C Parameter C of the momentum dependent part of the potential
   *              in MeV
   * \param[in] Lambda Parameter Lambda of the momentum-dependent term of
   *                   the potential in 1/fm
   * \return effective mass squared in calculation frame
   *         minus effective mass in rest_frame in GeV^2
   */
  static double root_eq_potentials(double energy_calc,
                                   const ThreeVector &momentum_calc,
                                   const FourVector &jmu, double m, double A,
                                   double B, double tau, double C,
                                   double Lambda) {
    // get velocity for boost to the local rest frame
    double rho_LRF = jmu.abs();
    ThreeVector beta_LRF = jmu.x0() > really_small ? jmu.threevec() / jmu.x0()
                                                   : ThreeVector(0, 0, 0);
    // get momentum in the local rest frame
    FourVector pmu_calc = FourVector(energy_calc, momentum_calc);
    FourVector pmu_LRF = beta_LRF.abs() > really_small
                             ? pmu_calc.lorentz_boost(beta_LRF)
                             : pmu_calc;
    double p_LRF = pmu_LRF.threevec().abs();
    double energy_LRF = std::sqrt(m * m + p_LRF * p_LRF) +
                        skyrme_pot(rho_LRF, A, B, tau) +
                        momentum_dependent_part(p_LRF, rho_LRF, C, Lambda);
    const double result = energy_calc * energy_calc - momentum_calc.sqr() -
                          (energy_LRF * energy_LRF - p_LRF * p_LRF);
    logg[LPotentials].debug()
        << "root equation for potentials called with E_calc=" << energy_calc
        << " p_calc=" << momentum_calc << " jmu=" << jmu << " m=" << m
        << " tau=" << tau << " A=" << A << " B=" << B << " C=" << C
        << " Lambda=" << Lambda << " and the root equation is " << result;
    return result;
  }

  /**
   * Momentum dependent term of the potential
   *
   * To be added to the momentum independent part.
   *
   * \param[in] momentum Absolute momentum of the particle of interest in the
   *                     local rest-frame in GeV
   * \param[in] rho Baryon density in the Eckart frame in 1/fm^3
   * \param[in] C Parameter C of the momentum dependent part of the
   *              potential in MeV
   * \param[in] Lambda Parameter Lambda of the momentum-dependent part of the
   *                   potential in 1/fm
   * \return momentum dependent part of the potential in GeV
   */
  static double momentum_dependent_part(double momentum, double rho, double C,
                                        double Lambda) {
    /* We assume here that the distribution function is the one of cold
     * nuclear matter, which consists only of protons and neutrons.
     * That is, the degeneracy factor simply equals nucleon spin degeneracy
     * times nucleon isospin degeneracy, even though all baryons contribute
     * to the density and the potential is applied to all baryons. */
    int g = 4;
    const double fermi_momentum =
        std::cbrt(6. * M_PI * M_PI * rho / g);  // in 1/fm
    momentum = momentum / hbarc;                // convert to 1/fm
    if (unlikely(momentum < really_small)) {
      return mev_to_gev * g * C / (M_PI * M_PI * nuclear_density) *
             (Lambda * Lambda * fermi_momentum -
              std::pow(Lambda, 3) * std::atan(fermi_momentum / Lambda));
    }
    const std::array<double, 7> temp = {
        2 * g * C * M_PI * std::pow(Lambda, 3) /
            (std::pow(2 * M_PI, 3) * nuclear_density),
        (fermi_momentum * fermi_momentum + Lambda * Lambda -
         momentum * momentum) /
            (2 * momentum * Lambda),
        std::pow(momentum + fermi_momentum, 2) + Lambda * Lambda,
        std::pow(momentum - fermi_momentum, 2) + Lambda * Lambda,
        2 * fermi_momentum / Lambda,
        (momentum + fermi_momentum) / Lambda,
        (momentum - fermi_momentum) / Lambda};
    const double result =
        temp[0] * (temp[1] * std::log(temp[2] / temp[3]) + temp[4] -
                   2 * (std::atan(temp[5]) - std::atan(temp[6])));
    return mev_to_gev * result;
  }
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_POTENTIALS_H_
