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

#include "experimentparameters.h"
#include "forwarddeclarations.h"
#include "fourvector.h"
#include "lattice.h"
#include "particledata.h"
#include "pdgcode.h"
#include "threevector.h"

namespace Smash {

  std::ostream& operator<<(std::ostream& os, DensityType dt);

  /** Get the factor that determines how much a particle contributes to the
   *  density type that is computed. E.g. positive pion contributes with
   *  factor 1 to total particle density and with factor 0 to baryon density.
   *  Proton contributes with factor 1 to baryon density, anti-proton - with
   *  factor -1 to baryon density, and so on.
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
   *  \f$ (2 \pi \sigma^2)^{3/2} \f$. Division over norm is splitted
   *  for efficiency: it is not nice to recalculate the same constant
   *  norm at every call.
   *
   * \param[in] r vector from the particle to the point of interest
   * \param[in] p particle 4-momentum to account for Lorentz contraction
   * \param[in] m particle mass, \f$ m = \sqrt{E^2 - p^2} \f$
   * \param[in] two_sigma_sqr \f$ 2 \sigma^2 \f$,
   *            \f$ \sigma \f$ - width of gaussian smearing
   * \param[in] r_cut_sqr radius, where gaussian is cut, squared
   * \param[in] compute_gradient option, true - compute gradient, false - no
   * \return smearing factor itself and optionally also its gradient
   */
  std::pair<double, ThreeVector> unnormalized_smearing_factor(
                       const ThreeVector &r, const FourVector &p,
                       const double m,
                       const double two_sigma_sqr, const double r_cut_sqr,
                       const bool compute_gradient = false);
  /**
   * Norm of the smearing function, \f$ (2 \pi \sigma^2)^{3/2}\f$
   *
   * \param[in] two_sigma_sqr \f$2 \sigma^2 \f$,
   *            \f$ \sigma \f$ - width of gaussian smearing
   */
  inline double smearing_factor_norm(const double two_sigma_sqr) {
    const double tmp = two_sigma_sqr * M_PI;
    return tmp * std::sqrt(tmp);
  }

  /**
   * Norm of the smearing factor gradient, \f$ (2 \pi \sigma^2)^{3/2} \cdot
   * (2 \sigma^2) \f$
   *
   * \param[in] two_sigma_sqr \f$2 \sigma^2 \f$,
   *            \f$ \sigma \f$ - width of gaussian smearing
   */
  inline double smearing_factor_grad_norm(const double two_sigma_sqr) {
    const double tmp = two_sigma_sqr * M_PI;
    return tmp * std::sqrt(tmp) * 0.5 * two_sigma_sqr;
  }

  /**
   * Gaussians used for smearing are cut at radius \f$r_{cut} = a \sigma \f$
   * for calculation speed-up. In the limit of \f$a \to \infty \f$ smearing
   * factor is normalized to 1:
   * \f[ \frac{4 \pi}{(2 \pi \sigma^2)^{3/2}}
   *     \int_0^{\infty} e^{-r^2/2 \sigma^2} r^2 dr = 1 \f]
   * However, for finite \f$ a\f$ integral is less than one:
   * \f[ g(a) \equiv \frac{4 \pi}{(2 \pi \sigma^2)^{3/2}}
   *    \int_0^{a \sigma} e^{-r^2/2 \sigma^2} r^2 dr =
   *    -\sqrt{\frac{2}{\pi}} a e^{-a^2/2} + Erf[a/\sqrt{2}]
   * \f] This \f$ g(a) \f$ is typically close to 1. For example,
   * for \f$r_{cut} = 3 \sigma \f$, and thus \f$ a=3 \f$, g(3) = 0.9707;
   * g(4) = 0.9987. The aim of this function is to compensate for this factor.
   *
   * \param[in] rcut_in_sigma \f$ a = r_{cut} / \sigma\f$
   * \return \f$ g(a) \f$
   */
  inline float smearing_factor_rcut_correction(const float rcut_in_sigma) {
    const float x = rcut_in_sigma / std::sqrt(2.0);
    return - 2.0 /std::sqrt(M_PI) * x * std::exp(-x*x) + std::erf(x);
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
   *  To avoid the problems with Eckart frame definition, densities for
   *  positive and negative charges, \f$\rho_+ \f$ and \f$ \rho_-\f$,
   *  are computed separately and result is \f$\rho_+ - \rho_-\f$.
   *
   * \param[in] r Arbitrary space point where 4-current is calculated
   * \param[in] plist List of all particles to be used in \f$j^{\mu}\f$
   *            calculation. If the distance between particle and calculation
   *            point r, \f$ |r-r_i| > r_{cut} \f$ then particle input
   *            to density will be ignored.
   *
   * Next three values are taken from ExperimentalParameters structure:
   * \param[in] par Set of parameters packed in one structure.
   *            From them the cutting radius r_cut \f$ r_{cut} / \sigma \f$,
   *            number of test-particles ntest and the gaussian width
   *            gs_sigma are needed.
   * \param[in] dens_type type of four-currect to be calculated:
   *            baryon, proton or neutron options are currently available
   * \param[in] compute_gradient true - compute gradient, false - no
   * \fpPrecision Density gradient is returned as double, because it is
   *   ThreeVector and ThreeVector currently comes only as double.
   *   Density itself is double for uniformity: if gradient is double,
   *   density should also be.
   */
  std::pair<double, ThreeVector> rho_eckart(const ThreeVector &r,
                const ParticleList &plist, const ExperimentParameters &par,
                DensityType dens_type, bool compute_gradient);
  /// convenience overload of the above
  std::pair<double, ThreeVector> rho_eckart(const ThreeVector &r,
                const Particles &plist, const ExperimentParameters &par,
                DensityType dens_type, bool compute_gradient);


/** A class for time-efficient (time-memory trade-off) calculation of density
 *  on the lattice. It holds two FourVectors - positive and negative
 *  summands of 4-current, and the density itself. Four-currents are
 *  additive by particles, density is not. That is why such structure is
 *  used. It is efficient to calculate additive jmu's in one loop over
 *  particles and then calculate density from jmu's. Splitting into
 *  positive and negative parts of jmu is necessary to avoid
 *  problems with the definition of Eckart rest frame.
 *
 *  Intended usage of the class:
 *  -# Add particles from some list using add_particle(...). This sets
 *     jmu_pos and jmu_neg
 *  -# Compute density from jmus using compute_density(...).
 *  -# Get jmus and density whenever necessary via density(),
 *     jmu_pos(), jmu_neg()
 */
class DensityOnLattice {
 public:
  /// Default constructor
  DensityOnLattice() : jmu_pos_(FourVector()),
                       jmu_neg_(FourVector()),
                       density_(0.0) {}
  /** Adds particle to 4-current: \f$j^{\mu} += p^{\mu}/p^0 \cdot factor \f$
   *  Factor can in principle be any scalar multiplier. Physics-wise it
   *  accounts for smearing on the lattice and for particle contribution
   *  to given density type (e.g. anti-proton contributes with factor -1
   *  to baryon density, proton - with factor 1).
   */
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
  /// Returns the density if it was previously computed
  double density() const { return density_; }
  /// Returns positive part of current four-vector
  FourVector jmu_pos() const { return jmu_pos_; }
  /// Returns negative part of current four-vector
  FourVector jmu_neg() const { return jmu_neg_; }

 private:
  /// Positive and negative parts of four-current
  FourVector jmu_pos_, jmu_neg_;
  /// The density itself
  double density_;
};

  /// Conveniency typedef for lattice of density
  typedef RectangularLattice<DensityOnLattice> DensityLattice;

  /** Calculates density on the lattice in an time-efficient way.
   *  \param lat pointer to the lattice
   *  \param update tells if called for update at printout or at timestep
   *  \param dens_type density type to be computed on the lattice
   *  \param par a structure containing testparticles number and gaussian
   *         smearing parameters.
   *  \param particles the particles vector
   */
  void update_density_lattice(DensityLattice* lat,
                            const LatticeUpdate update,
                            const DensityType dens_type,
                            const ExperimentParameters &par,
                            const Particles &particles);
}  // namespace Smash

#endif  // SRC_INCLUDE_DENSITY_H_
