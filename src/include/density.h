/*
 *
 *    Copyright (c) 2013-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#ifndef SRC_INCLUDE_DENSITY_H_
#define SRC_INCLUDE_DENSITY_H_

#include <utility>
#include <vector>
#include <iostream>

#include "forwarddeclarations.h"
#include "fourvector.h"
#include "particledata.h"
#include "pdgcode.h"
#include "threevector.h"

namespace Smash {

  /** Allows to choose which kind of density to calculate.
   *  The baryon density is necessary for the Skyrme potential.
   *  For the symmetry potential one needs to know the isospin density.
   */
  enum class DensityType {
    particle = 0,
    baryon = 1,
    baryonic_isospin = 2,
    pion = 3,
    none = 4,
  };

  std::ostream& operator<<(std::ostream& os, DensityType dt);

  /** Get the factor that determines how much a particle contributes to the
   *  density type that is computed.
   *
   *  \param pdg PDG code of particle to be tested
   *  \param dens_type The density type
   *
   *  \return The corresponding factor (0 if the particle doesn't
   *          contribute at all).
   */
  float density_factor(const PdgCode pdg, DensityType dens_type);

  /**
   * Implements gaussian smearing for any quantity.
   * Computes smearing factor taking Lorentz contraction into account.
   * Integral of unnormalized smearing factor over space should be
   *  \f[ (2 \pi \sigma^2)^{3/2}. \f]. Division over norm is splitted
   *  for efficiency: it is not nice to recalculate the same constant
   *  norm at every call.
   *
   *  Returns smearing factor itself and optionally also its gradient.
   *
   * \param[in] r vector from the particle to the point of interest
   * \param[in] p particle momentum to account for Lorentz contraction
   * \param[in] two_sigma_sqr 2*sigma^2, sigma - width of gaussian smearing
   * \param[in] r_cut_sqr radius, where gaussian is cut, squared
   * \param[in] compute_gradient option, true - compute gradient, false - no.
   */
  std::pair<double, ThreeVector> unnormalized_smearing_factor(
                       const ThreeVector &r, const FourVector &p,
                       const double two_sigma_sqr, const double r_cut_sqr,
                       const bool compute_gradient = false);
  /**
   * Norm of the smearing function, \f[ (2 \pi \sigma^2)^{3/2}. \f].
   *
   * \param[in] two_sigma_sqr 2*sigma^2, sigma - width of gaussian smearing
   */
  inline double smearing_factor_norm(const double two_sigma_sqr) {
    const double tmp = two_sigma_sqr * M_PI;
    return tmp * std::sqrt(tmp);
  }

  /**
   * Norm of the smearing factor gradient, \f[ (2 \pi \sigma^2)^{3/2} *
   * (2 \sigma^2). \f].
   *
   * \param[in] two_sigma_sqr 2*sigma^2, sigma - width of gaussian smearing
   */
  inline double smearing_factor_grad_norm(const double two_sigma_sqr) {
    const double tmp = two_sigma_sqr * M_PI;
    return tmp * std::sqrt(tmp) * 0.5 * two_sigma_sqr;
  }

  /** Calculates Eckart rest frame density and optionally its gradient.
   *  \f[j^{\mu} = (\sqrt{2\pi} \sigma )^{-3} \sum_{i=1}^N C_i u^{\mu}_i
   *  exp \left(- \frac{(\vec r -\vec r_i + \frac{\gamma_i^2}{1 + \gamma_i}
   *  \vec \beta_i (\vec \beta_i, \vec r - \vec r_i))^2}{2\sigma^2} \right)\f]
   *  \f[ \rho^{Eckart} = \sqrt{j^{\mu} j_{\mu}} \f]
   *  Here \f$ C_i \f$ is a corresponding value of "charge". If baryon
   *  current option is selected then \f$ C_i \f$ is 1 for baryons,
   *  -1 for antibaryons and 0 otherwize. For proton/neutron
   *  current \f$ C_i = 1\f$ for proton/neutron and 0 otherwize.
   *
   *  For gradient:
   *  \f[ \frac{d\rho_{Eck}}{d \vec r} = \frac{\frac{dj^{\mu}}{d \vec r}
   *  j_{\mu}}{\sqrt{j^{\mu}j_{\mu}}} \f]
   *
   * \param[in] r Arbitrary space point where 4-current is calculated
   * \param[in] plist List of all particles to be used in \f$j^{\mu}\f$
   *            calculation. If the distance between particle and calculation
   *            point r, \f$ |r-r_i| > 6 \sigma \f$ then particle input
   *            to density will be ignored.
   * \param[in] gs_sigma Width of the gaussian (\f$ \sigma \f$),
   *  which represents particle Wigner density.
   * \param[in] dens_type type of four-currect to be calculated:
   *            baryon, proton or neutron options are currently available
   * \param[in] ntest Number of test-particles
   * \param[in] compute_gradient true - compute gradient, false - no
   * \fpPrecision Density gradient is returned as double, because it is
   *   ThreeVector and ThreeVector currently comes only as double.
   *   Density itself is double for uniformity: if gradient is double,
   *   density should also be.
   */
  std::pair<double, ThreeVector> rho_eckart(const ThreeVector &r,
                    const ParticleList &plist,
                    double gs_sigma, DensityType dens_type, int ntest,
                    bool compute_gradient);
  /// convenience overload of the above
  std::pair<double, ThreeVector> rho_eckart(const ThreeVector &r,
                    const Particles &plist,
                    double gs_sigma, DensityType dens_type, int ntest,
                    bool compute_gradient);

/** A class for time-efficient (time-memory trade-off) calculation of density
 *  on the lattice. It holds two FourVectors - positive and negative
 *  summands of 4-current, and the density itself. Four-currents are
 *  additive by particles, density is not. That is why such structure is
 *  used. It is efficient to calculate additive jmu's in one loop over
 *  particles and then calculate density from jmu's. Splitting into
 *  positive and negative parts of jmu is necessary to avoid
 *  problems with the definition of Eckart rest frame.
 */
class DensityOnLattice {
 public:
  /// Default constructor
  DensityOnLattice() : jmu_pos_(FourVector()),
                       jmu_neg_(FourVector()),
                       density_(0.0) {}
  /// Adds particle to 4-current: \f$j^{\mu} += p^{\mu}/p^0 * factor \f$
  void add_particle(const ParticleData &part, double factor) {
    if (factor > 0.0) {
      jmu_pos_ += FourVector(factor, part.velocity() * factor);
    } else {
      jmu_neg_ += FourVector(factor, part.velocity() * factor);
    }
  }
  /// Computes density from jmu
  void compute_density(const double norm) {
    density_ = (jmu_pos_.abs() - jmu_neg_.abs()) / norm;
  }
  /// Getter for density
  double density() const { return density_; }
  /// Getter for jmu_pos
  FourVector jmu_pos() const { return jmu_pos_; }
  /// Getter for jmu_neg
  FourVector jmu_neg() const { return jmu_neg_; }
 private:
  /// Positive and negative parts of four-current
  FourVector jmu_pos_, jmu_neg_;
  /// The density itself
  double density_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_DENSITY_H_
