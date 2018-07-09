/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_POTENTIALS_H_
#define SRC_INCLUDE_POTENTIALS_H_

#include <utility>
#include <vector>

#include "configuration.h"
#include "density.h"
#include "forwarddeclarations.h"
#include "particledata.h"
#include "threevector.h"

#ifdef BUILD_TESTS
#define VIRTUAL_FOR_TESTS virtual
#else
#define VIRTUAL_FOR_TESTS
#endif

namespace smash {

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
   *            symmetry potentials, and the coefficents controlling
   *            how strong the potentials are.
   * \param[in] parameters Struct that contains the gaussian smearing factor
   *            \f$\sigma\f$, the distance cutoff \f$r_{\rm cut}\f$ and
   *            the testparticle number needed for the density calculation.
   */
  Potentials(Configuration conf, const DensityParameters &parameters);
  VIRTUAL_FOR_TESTS
  /// Standard destructor
  ~Potentials();

  /**
   * Evaluates skyrme potential given a baryon density.
   *
   * \param[in] baryon_density Baryon density \f$\rho\f$ evaluated in the
   *            local rest frame in fm\f$^{-3}\f$.
   * \return Skyrme potential \f[U_B=10^{-3}\times\frac{\rho}{|\rho|}
   *         (A\frac{\rho}{\rho_0}+B(\frac{\rho}{\rho_0})^\tau)\f] in GeV
   */
  double skyrme_pot(const double baryon_density) const;

  /**
   * Evaluates symmetry potential given baryon isospin density.
   *
   * \param[in] baryon_isospin_density The difference between the proton and
   *            the neutron density in the local rest frame in fm\f$^{-3}\f$.
   * \return Symmetry potential \f[U_I=2\times 10^{-3}S_{\rm sym}
   *         \frac{\rho_n-\rho_p}{\rho_0}\f] in GeV
   */
  double symmetry_pot(const double baryon_isospin_density) const;

  /**
   * Evaluates potential at point r. Potential is always taken in the local
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
   * \return Total potential energy acting on the particle: \f[U_{\rm tot}
   *         =Q_BU_B+2I_3U_I\f] in GeV, where \f$Q_B\f$ is the baryon charge
   *	     scaled by the ratio of the light (u, d) quark to the total quark
   *         number and \f$I_3\f$ is the third compnent of the isospin.
   */
  VIRTUAL_FOR_TESTS
  double potential(const ThreeVector &r, const ParticleList &plist,
                   const ParticleType &acts_on) const;

  /**
   * Evaluates the scaling factor of the forces acting on the particles. The
   * forces are equal to the product of the scaling factor and the gradient of
   * the potential. We need these scaling factors to describe the motions of
   * the hyperons as well as the anti-particles in the potentials. For Lambda
   * and Sigma, since they carry 2 light (u or d) quarks, they are affected by
   * 2/3 of the Skyrme force. Xi carries 1 light quark, it is affected by 1/3
   * of the Skyrme force. Omega carries no light quark, so it's not affected by
   * the Skyrme force. Anti-baryons are affected by the force as large as the
   * force acting on baryons but with an opposite direction.
   *
   * \param[in] data Type of particle on which potential is going to act.
   * \return (\f$Q_B(1-\frac{Q_S}{3}), 2I_3\f$) where \f$Q_B\f$ is the baryon
   *         charge, \f$Q_S\f$ is the strangeness, and \f$I_3\f$ is the third
   *         component of the isospin.
   */
  std::pair<double, int> force_scale(const ParticleType &data) const;

  /**
   * Evaluates potential gradient at point r. Potential is always taken in
   * the local Eckart rest frame, but point r is in the computational frame.
   *
   * \param[in] r Arbitrary space point where potential gradient is calculated
   * \param[in] plist List of all particles to be used in \f$j^{\mu}\f$
   *            calculation. If the distance between particle and calculation
   *            point r, \f$ |r-r_i| > r_{cut} \f$ then particle input
   *            to density will be ignored.
   * \return (\f$\nabla U_B, \nabla U_I\f$), where \f$U_B\f$ and \f$U_I\f$ are
   *         the Skyrme and symmetry potentials correspondingly.
   */
  VIRTUAL_FOR_TESTS
  std::pair<ThreeVector, ThreeVector> potential_gradient(
      const ThreeVector &r, const ParticleList &plist) const;

  /// \return Is Skyrme potential on?
  VIRTUAL_FOR_TESTS
  bool use_skyrme() const { return use_skyrme_; }
  /// \return Is symmetry potential on?
  VIRTUAL_FOR_TESTS
  bool use_symmetry() const { return use_symmetry_; }

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

  /// coefficent in front of the symmetry term.
  double symmetry_s_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_POTENTIALS_H_
