/*
 *
 *    Copyright (c) 2013-2021
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#ifndef SRC_INCLUDE_SMASH_DENSITY_H_
#define SRC_INCLUDE_SMASH_DENSITY_H_

#include <iostream>
#include <tuple>
#include <typeinfo>
#include <utility>
#include <vector>

#include "energymomentumtensor.h"
#include "experimentparameters.h"
#include "forwarddeclarations.h"
#include "fourvector.h"
#include "lattice.h"
#include "particledata.h"
#include "particles.h"
#include "pdgcode.h"
#include "threevector.h"

namespace smash {
static constexpr int LDensity = LogArea::Density::id;

/**
 * Allows to choose which kind of density to calculate.
 * The baryon density is necessary for the Skyrme potential.
 * For the symmetry potential one needs to know the isospin density.
 */
enum class DensityType {
  None = 0,
  Hadron = 1,
  Baryon = 2,
  BaryonicIsospin = 3,
  Pion = 4,
  Isospin3_tot = 5,
  Charge = 6,
  Strangeness = 7,
};

/**
 * Create the output operator for the densities
 *
 * \param[out] os Output operator for the densities
 * \param[in] dt Type of density (e.g. baryon density)
 * \return An output operator for the densities
 */
std::ostream &operator<<(std::ostream &os, DensityType dt);

/**
 * Get the factor that determines how much a particle contributes to the
 * density type that is computed. E.g. positive pion contributes with
 * factor 1 to total particle density and with factor 0 to baryon density.
 * Proton contributes with factor 1 to baryon density, anti-proton - with
 * factor -1 to baryon density, and so on.
 *
 * \param[in] type type of the particle to be tested
 * \param[in] dens_type The density type
 * \return The corresponding factor (0 if the particle doesn't
 *         contribute at all).
 */
double density_factor(const ParticleType &type, DensityType dens_type);

/**
 * Norm of the Gaussian smearing function
 *
 * \param[in] two_sigma_sqr \f$2 \sigma^2 \f$ [fm\f$^2\f$],
 *            \f$ \sigma \f$ - width of gaussian smearing
 * \return \f$ (2 \pi \sigma^2)^{3/2}\f$ [fm\f$^3\f$]
 */
inline double smearing_factor_norm(const double two_sigma_sqr) {
  const double tmp = two_sigma_sqr * M_PI;
  return tmp * std::sqrt(tmp);
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
inline double smearing_factor_rcut_correction(const double rcut_in_sigma) {
  const double x = rcut_in_sigma / std::sqrt(2.0);
  return -2.0 / std::sqrt(M_PI) * x * std::exp(-x * x) + std::erf(x);
}

/**
 * A class to pre-calculate and store parameters relevant for density
 * calculation. It has to be initialized only once per SMASH run.
 */
class DensityParameters {
 public:
  /**
   * Constructor of DensityParameters.
   *
   * \param[in] par Struct containing the Gaussian smearing width \f$\sigma\f$,
   *            the cutoff factor \f$a\f$ where the cutoff radius
   *            \f$r_{\rm cut}=a\sigma\f$, the test-particle number, the number
   *            of ensembles, the mode of calculating the derivatives, the
   *            smearing mode, the central weight for Discrete smearing, the
   *            range (in units of lattice spacing) for Triangular smearing
   *            and the flag about using only participants or also spectators
   */
  DensityParameters(const ExperimentParameters &par)  // NOLINT
      : sig_(par.gaussian_sigma),
        r_cut_(par.gauss_cutoff_in_sigma * par.gaussian_sigma),
        ntest_(par.testparticles),
        nensembles_(par.n_ensembles),
        derivatives_(par.derivatives_mode),
        rho_derivatives_(par.rho_derivatives_mode),
        smearing_(par.smearing_mode),
        central_weight_(par.discrete_weight),
        triangular_range_(par.triangular_range),
        only_participants_(par.only_participants) {
    r_cut_sqr_ = r_cut_ * r_cut_;
    const double two_sig_sqr = 2 * sig_ * sig_;
    two_sig_sqr_inv_ = 1. / two_sig_sqr;
    const double norm = smearing_factor_norm(two_sig_sqr);
    const double corr_factor =
        smearing_factor_rcut_correction(par.gauss_cutoff_in_sigma);
    norm_factor_sf_ = 1. / (norm * ntest_ * nensembles_ * corr_factor);
  }
  /// \return Testparticle number
  int ntest() const { return ntest_; }
  /// \return Number of ensembles
  int nensembles() const { return nensembles_; }
  /// \return Mode of gradient calculation
  DerivativesMode derivatives() const { return derivatives_; }
  /// \return Mode of rest frame density derivatives (on or off)
  RestFrameDensityDerivativesMode rho_derivatives() const {
    return rho_derivatives_;
  }
  /// \return Smearing mode
  SmearingMode smearing() const { return smearing_; }
  /// \return Weight of the central cell in the discrete smearing
  double central_weight() const { return central_weight_; }
  /// \return Range of the triangular smearing, in units of lattice spacing
  double triangular_range() const { return triangular_range_; }
  /// \return Cut-off radius [fm]
  double r_cut() const { return r_cut_; }
  /// \return Squared cut-off radius [fm\f$^2\f$]
  double r_cut_sqr() const { return r_cut_sqr_; }
  /// \return \f$ (2 \sigma^2)^{-1} \f$ [fm\f$^{-2}\f$]
  double two_sig_sqr_inv() const { return two_sig_sqr_inv_; }
  /**
   * \return Normalization for smearing factor. Unnormalized smearing factor
   *         \f$ sf(\vec{r}) \f$ has to be multiplied by this to have
   *         \f$ \int d^3r \, sf(\vec{r}) = 1 \f$.
   */
  double norm_factor_sf() const { return norm_factor_sf_; }
  /// \return counting only participants (true) or also spectators (false)
  bool only_participants() const { return only_participants_; }

 private:
  /// Gaussian smearing width [fm]
  const double sig_;
  /// Cut-off radius [fm]
  const double r_cut_;
  /// Squared cut-off radius [fm\f$^2\f$]
  double r_cut_sqr_;
  /// \f$ (2 \sigma^2)^{-1} \f$ [fm\f$^{-2}\f$]
  double two_sig_sqr_inv_;
  /// Normalization for Gaussian smearing factor
  double norm_factor_sf_;
  /// Testparticle number
  const int ntest_;
  /// Number of ensembles
  const int nensembles_;
  /// Mode of calculating the gradients
  const DerivativesMode derivatives_;
  /// Whether to calculate the rest frame density derivatives
  const RestFrameDensityDerivativesMode rho_derivatives_;
  /// Mode of smearing
  const SmearingMode smearing_;
  /// Weight of the central cell in the discrete smearing
  const double central_weight_;
  /// Range of the triangular smearing
  const double triangular_range_;
  /// Flag to take into account only participants
  bool only_participants_;
};

/**
 * Implements gaussian smearing for any quantity.
 * Computes smearing factor taking Lorentz contraction into account.
 * Integral of unnormalized smearing factor over space should be
 *  \f$ (2 \pi \sigma^2)^{3/2} \f$. Division over norm is split
 *  for efficiency: it is not nice to recalculate the same constant
 *  norm at every call.
 *
 * \param[in] r vector from the particle to the point of interest [fm]
 * \param[in] p particle 4-momentum to account for Lorentz contraction [GeV]
 * \param[in] m_inv particle mass, \f$ (E^2 - p^2)^{-1/2} \f$ [GeV]
 * \param[in] dens_par object containing precomputed parameters for
 *            density calculation.
 * \param[in] compute_gradient option, true - compute gradient, false - no
 * \return (smearing factor, the gradient of the smearing factor or a zero
 *         three vector)
 */
std::pair<double, ThreeVector> unnormalized_smearing_factor(
    const ThreeVector &r, const FourVector &p, const double m_inv,
    const DensityParameters &dens_par, const bool compute_gradient = false);

/*!\Userguide
 * \anchor current_eckart
 * Calculates Eckart rest frame density and 4-current of a given density type
 * and optionally the gradient of the density in an arbitary frame (grad j0),
 * the curl of the 3-current, and the time, x, y, and z derivatives of the
 * 4-current.
 * \f[j^{\mu} = (\sqrt{2\pi} \sigma )^{-3} \sum_{i=1}^N C_i u^{\mu}_i
 * exp \left(- \frac{(\vec r -\vec r_i + \frac{\gamma_i^2}{1 + \gamma_i}
 * \vec \beta_i (\vec \beta_i, \vec r - \vec r_i))^2}{2\sigma^2} \right)\f]
 * \f[ \rho^{Eckart} = \sqrt{j^{\mu} j_{\mu}} \f]
 * Here \f$ C_i \f$ is a corresponding value of "charge". If baryon
 * current option is selected then \f$ C_i \f$ is 1 for baryons,
 * -1 for antibaryons and 0 otherwise. For proton/neutron
 * current \f$ C_i = 1\f$ for proton/neutron and 0 otherwise.
 *
 * To avoid the problems with Eckart frame definition, densities for
 * positive and negative charges, \f$\rho_+ \f$ and \f$ \rho_-\f$,
 * are computed separately and final density is \f$\rho_+ - \rho_-\f$.
 *
 * \param[in] r Arbitrary space point where 4-current is calculated [fm];
              ignored if smearing is false
 * \param[in] plist List of all particles to be used in \f$j^{\mu}\f$
 *            calculation. If smearing is false or if the distance
 *            between particle and calculation point r,
 *            \f$ |r-r_i| > r_{cut} \f$ then particle input
 *            to density will be ignored.
 *
 * Next four values are taken from ExperimentalParameters structure:
 *
 * \param[in] par Set of parameters packed in one structure.
 *            From them the cutting radius r_cut \f$ r_{cut} / \sigma \f$,
 *            number of test-particles ntest and the gaussian width
 *            gs_sigma are needed.
 * \param[in] dens_type type of four-currect to be calculated:
 *            baryon, proton or neutron options are currently available
 * \param[in] compute_gradient true - compute gradient, false - no
 * \param[in] smearing whether to use gaussian smearing or not. If false,
 *            this parameter will use ALL particles equally to calculate the
 *            current, and that as such it will not be normalized wrt volume.
 *            This should be true for any internal calculation of any quantity
 *            and only makes sense to turn off for output purposes in a box.
 * \return (rest frame density in the local Eckart frame [fm\f$^{-3}\f$],
 *          \f$ j^\mu \f$ as a 4-vector,
 *          \f$ \vec{\nabla}\cdot j^0 \f$ or a 0 3-vector,
 *          \f$ \vec{\nabla} \times \vec{j} \f$ or a 0 3-vector,
 *          \f$ \partial_t j^\mu \f$ or a 0 4-vector,
 *          \f$ \partial_x j^\mu \f$ or a 0 4-vector,
 *          \f$ \partial_y j^\mu \f$ or a 0 4-vector,
 *          \f$ \partial_z j^\mu \f$ or a 0 4-vector).
 */
std::tuple<double, FourVector, ThreeVector, ThreeVector, FourVector, FourVector,
           FourVector, FourVector>
current_eckart(const ThreeVector &r, const ParticleList &plist,
               const DensityParameters &par, DensityType dens_type,
               bool compute_gradient, bool smearing);
/// convenience overload of the above (ParticleList -> Particles)
std::tuple<double, FourVector, ThreeVector, ThreeVector, FourVector, FourVector,
           FourVector, FourVector>
current_eckart(const ThreeVector &r, const Particles &plist,
               const DensityParameters &par, DensityType dens_type,
               bool compute_gradient, bool smearing);

/**
 * A class for time-efficient (time-memory trade-off) calculation of density
 * on the lattice. It holds six FourVectors - positive and negative
 * summands of 4-current, and the time and spatial derivatives of the compound
 * current. These four-vectors are  additive by particles. It is efficient to
 * calculate additive \f$j^\mu\f$ and \f$\partial_\nu j^\mu \f$ in one loop over
 * particles and then calculate the Eckart density, the gradient of the density,
 * the curl, the time derivative of the current, and derivatives of the rest
 * frame density accordingly.
 * Splitting into  positive and negative parts of \f$j^\mu\f$ is necessary to
 * avoid problems with the definition of Eckart rest frame.
 *
 * Intended usage of the class:
 * -# Add particles from some list using add_particle(), setting jmu_pos and
 *    jmu_neg. Calculate derivatives using either add_particle_for_derivatives()
 *    (in case of Gaussian derivatives) or calculating finite difference
 *    derivatives; this sets djmu_dxnu. If needed, calculate rest frame density
 *    derivatives, setting drho_dxnu.
 * -# Get the net current via jmu_net().
 * -# Get the net rest frame density via rho().
 * -# Get the derivatives of the net current via djmu_dxnu()
 * -# Get the derivatives of the net rest frame baryon density via drho_dxnu()
 * -# Get \f$\vec{\nabla} j^0\f$ via grad_j0()
 * -# Get \f$\vec{\nabla} \times \vec{j}\f$ via curl_vecj()
 * -# Get \f$\partial_t \vec{j}\f$ via dvecj_dt()
 * -# Get \f$(\vec{\nabla} \rho) \times \vec{j}\f$ via grad_rho_cross_vecj()
 */
class DensityOnLattice {
 public:
  /// Default constructor
  DensityOnLattice()
      : jmu_pos_(FourVector()),
        jmu_neg_(FourVector()),
        djmu_dxnu_({FourVector(), FourVector(), FourVector(), FourVector()}),
        drho_dxnu_(FourVector()) {}

  /**
   * Adds particle to 4-current: \f$j^{\mu} += p^{\mu}/p^0 \cdot factor \f$.
   * Two private class members jmu_pos_ and jmu_neg_ indicating the 4-current
   * of the positively and negatively charged particles are updated by this
   * function.
   *
   * \param[in] part Particle would be added to the current density
   *            on the lattice.
   * \param[in] FactorTimesSf particle contribution to given density type (e.g.
   *            anti-proton contributes with factor -1 to baryon density,
   *            proton - with factor 1) times the smearing factor.
   */
  void add_particle(const ParticleData &part, double FactorTimesSf) {
    const FourVector part_four_velocity = FourVector(1.0, part.velocity());
    if (FactorTimesSf > 0.0) {
      jmu_pos_ += part_four_velocity * FactorTimesSf;
    } else {
      jmu_neg_ += part_four_velocity * FactorTimesSf;
    }
  }

  /**
   * Adds particle to the time and spatial derivatives of the 4-current.
   * An array of four private 4-vectors djmu_dxnu_ indicating the derivatives
   * of the compound current are updated by this function.
   *
   * \param[in] part Particle would be added to the current density
   *            on the lattice.
   * \param[in] factor particle contribution to given density type (e.g.
   *            anti-proton contributes with factor -1 to baryon density,
   *            proton - with factor 1).
   * \param[in] sf_grad Smearing factor of the gradients
   */
  void add_particle_for_derivatives(const ParticleData &part, double factor,
                                    ThreeVector sf_grad) {
    const FourVector PartFourVelocity = FourVector(1.0, part.velocity());
    for (int k = 1; k <= 3; k++) {
      djmu_dxnu_[k] += factor * PartFourVelocity * sf_grad[k - 1];
      djmu_dxnu_[0] -=
          factor * PartFourVelocity * sf_grad[k - 1] * part.velocity()[k - 1];
    }
  }

  /**
   * Compute the net Eckart density on the local lattice
   *
   * Note that the net Eckart density is calculated by taking the difference
   * between the Eckart density of the positively charged particles and that
   * of the negatively charged particles, which are, in general, defined in
   * different frames. So the net Eckart density is not the net density in the
   * Eckart local rest frame. However, this is the only way we can think of
   * to be applied to the case where the density current is space-like. And
   * fortunately, the net eckart densities are only used for calculating the
   * potentials which are valid only in the low-energy collisions where the
   * amount of the negatively charged particles are negligible. May be in the
   * future, the net Eckart density can be calculated in a smarter way.
   *
   * \param[in] norm_factor Normalization factor
   * \return Net Eckart density on the local lattice \f$\rho\f$ [fm\f$^{-3}\f$]
   */
  double rho(const double norm_factor = 1.0) {
    return (jmu_pos_.abs() - jmu_neg_.abs()) * norm_factor;
  }

  /**
   * Compute curl of the current on the local lattice
   *
   * \param[in] norm_factor Normalization factor
   * \return \f$\vec{\nabla}\times\vec{j}\f$ [fm \f$^{-4}\f$]
   */
  ThreeVector curl_vecj(const double norm_factor = 1.0) {
    ThreeVector curl_vec_j = ThreeVector();
    curl_vec_j.set_x1(djmu_dxnu_[2].x3() - djmu_dxnu_[3].x2());
    curl_vec_j.set_x2(djmu_dxnu_[3].x1() - djmu_dxnu_[1].x3());
    curl_vec_j.set_x3(djmu_dxnu_[1].x2() - djmu_dxnu_[2].x1());
    curl_vec_j *= norm_factor;
    return curl_vec_j;
  }

  /**
   * Compute gradient of the the zeroth component of the four-current j^mu
   * (that is of the computational frame density) on the local lattice
   *
   * \param[in] norm_factor Normalization factor
   * \return \f$\vec{\nabla} j^0\f$ [fm \f$^{-4}\f$]
   */
  ThreeVector grad_j0(const double norm_factor = 1.0) {
    ThreeVector j0_grad = ThreeVector();
    for (int i = 1; i < 4; i++) {
      j0_grad[i - 1] = djmu_dxnu_[i].x0() * norm_factor;
    }
    return j0_grad;
  }

  /**
   * Compute time derivative of the current density on the local lattice
   *
   * \param[in] norm_factor Normalization factor
   * \return \f$\partial_t \vec j\f$ [fm \f$^{-4}\f$]
   */
  ThreeVector dvecj_dt(const double norm_factor = 1.0) {
    return djmu_dxnu_[0].threevec() * norm_factor;
  }

  /**
   * \return Net current density
   *
   * There is a "+" operator in between, because the negative symbol
   * of the charge has already be included in FactorTimesSF.
   */
  FourVector jmu_net() const { return jmu_pos_ + jmu_neg_; }

  /**
   * Add to the positive density current.
   * \param[in] additional_jmu_B Value of positive density current to be added
   */
  void add_to_jmu_pos(FourVector additional_jmu_B) {
    jmu_pos_ += additional_jmu_B;
  }

  /**
   * Add to the negative density current.
   * \param[in] additional_jmu_B Value of negative density current to be added
   */
  void add_to_jmu_neg(FourVector additional_jmu_B) {
    jmu_neg_ += additional_jmu_B;
  }

  /**
   * Return the FourGradient of the rest frame density
   * \f$\partial_{\nu}\rho\f$
   * \return the FourGradient of the rest frame density
   *         \f$\partial_{\nu}\rho\f$
   */
  FourVector drho_dxnu() const { return drho_dxnu_; }

  /**
   * Return the FourGradient of the net baryon current
   * \f$\partial_{\nu} j^\mu\f$
   * \return the array of FourGradients of \f$\partial_{\nu} j^\mu\f$
   */
  std::array<FourVector, 4> djmu_dxnu() const { return djmu_dxnu_; }

  /**
   * Compute the  cross product of \f$\vec{\nabla}\rho\f$ and \f$j^\mu\f$
   * \return the cross product of \f$\vec{\nabla} \rho\f$ and \f$\vec{j}\f$
   */
  ThreeVector grad_rho_cross_vecj() const {
    const ThreeVector grad_rho = drho_dxnu_.threevec();
    const ThreeVector vecj = jmu_net().threevec();
    const ThreeVector Drho_cross_vecj = grad_rho.cross_product(vecj);

    return Drho_cross_vecj;
  }

  /**
   * Overwrite the time derivative of the current to zero.
   */
  void overwrite_djmu_dt_to_zero() {
    djmu_dxnu_[0] = FourVector(0.0, 0.0, 0.0, 0.0);
  }

  /**
   * Overwrite the time derivative of the rest frame density to zero.
   */
  void overwrite_drho_dt_to_zero() { drho_dxnu_[0] = 0.0; }

  /**
   * Overwrite the rest frame density derivatives to provided values.
   * \param[in] computed_drho_dxnu a FourGradient of the rest frame density rho
   */
  void overwrite_drho_dxnu(FourVector computed_drho_dxnu) {
    drho_dxnu_ = computed_drho_dxnu;
  }

  /**
   * Overwrite all density current derivatives to provided values.
   * \param[in] djmu_dt time derivative of the current FourVector jmu
   * \param[in] djmu_dx x derivative of the current FourVector jmu
   * \param[in] djmu_dy y derivative of the current FourVector jmu
   * \param[in] djmu_dz z derivative of the current FourVector jmu
   */
  void overwrite_djmu_dxnu(FourVector djmu_dt, FourVector djmu_dx,
                           FourVector djmu_dy, FourVector djmu_dz) {
    djmu_dxnu_[0] = djmu_dt;
    djmu_dxnu_[1] = djmu_dx;
    djmu_dxnu_[2] = djmu_dy;
    djmu_dxnu_[3] = djmu_dz;
  }

 private:
  /// Four-current density of the positively charged particle.
  FourVector jmu_pos_;
  /// Four-current density of the negatively charged particle.
  FourVector jmu_neg_;
  /// Four-gradient of the four-current density, \f$\partial_\nu j^\mu \f$
  std::array<FourVector, 4> djmu_dxnu_;
  /// Four-gradient of the rest frame density, \f$\partial_\nu \rho \f$
  FourVector drho_dxnu_;
};

/// Conveniency typedef for lattice of density
typedef RectangularLattice<DensityOnLattice> DensityLattice;

/**
 * Updates the contents on the lattice.
 *
 * \param[out] lat The lattice on which the content will be updated
 * \param[in] update tells if called for update at printout or at timestep
 * \param[in] dens_type density type to be computed on the lattice
 * \param[in] par a structure containing testparticles number and gaussian
 *            smearing parameters.
 * \param[in] ensembles the particles vector for each ensemble
 * \param[in] compute_gradient Whether to compute the gradients
 * \tparam T LatticeType
 */
template <typename T>
void update_lattice(RectangularLattice<T> *lat, const LatticeUpdate update,
                    const DensityType dens_type, const DensityParameters &par,
                    const std::vector<Particles> &ensembles,
                    const bool compute_gradient) {
  // Do not proceed if lattice does not exists/update not required
  if (lat == nullptr || lat->when_update() != update) {
    return;
  }

  lat->reset();
  // get the normalization factor for the covariant Gaussian smearing
  const double norm_factor_gaus = par.norm_factor_sf();
  // get the volume of the cell and weights for discrete smearing
  const double V_cell =
      (lat->cell_sizes())[0] * (lat->cell_sizes())[1] * (lat->cell_sizes())[2];
  // weights for coarse smearing
  const double big = par.central_weight();
  const double small = (1.0 - big) / 6.0;
  // get the radii for triangular smearing
  const std::array<double, 3> triangular_radius = {
      par.triangular_range() * (lat->cell_sizes())[0],
      par.triangular_range() * (lat->cell_sizes())[1],
      par.triangular_range() * (lat->cell_sizes())[2]};
  const double prefactor_triangular =
      1.0 /
      (par.ntest() * par.nensembles() * triangular_radius[0] *
       triangular_radius[0] * triangular_radius[1] * triangular_radius[1] *
       triangular_radius[2] * triangular_radius[2]);

  for (const Particles &particles : ensembles) {
    for (const ParticleData &part : particles) {
      if (par.only_participants()) {
        // if this conditions holds, the hadron is a spectator
        if (part.get_history().collisions_per_particle == 0) {
          continue;
        }
      }
      const double dens_factor = density_factor(part.type(), dens_type);
      if (std::abs(dens_factor) < really_small) {
        continue;
      }
      const FourVector p_mu = part.momentum();
      const ThreeVector pos = part.position().threevec();

      // act accordingly to which smearing is used
      if (par.smearing() == SmearingMode::CovariantGaussian) {
        const double m = p_mu.abs();
        if (unlikely(m < really_small)) {
          logg[LDensity].warn("Gaussian smearing is undefined for momentum ",
                              p_mu);
          continue;
        }
        const double m_inv = 1.0 / m;

        // unweighted contribution to density
        const double common_weight = dens_factor * norm_factor_gaus;
        lat->iterate_in_cube(
            pos, par.r_cut(), [&](T &node, int ix, int iy, int iz) {
              // find the weight for smearing
              const ThreeVector r = lat->cell_center(ix, iy, iz);
              const auto sf = unnormalized_smearing_factor(
                  pos - r, p_mu, m_inv, par, compute_gradient);
              node.add_particle(part, sf.first * common_weight);
              if (par.derivatives() == DerivativesMode::CovariantGaussian) {
                node.add_particle_for_derivatives(part, dens_factor,
                                                  sf.second * norm_factor_gaus);
              }
            });
      } else if (par.smearing() == SmearingMode::Discrete) {
        // unweighted contribution to density
        const double common_weight =
            dens_factor / (par.ntest() * par.nensembles() * V_cell);
        lat->iterate_nearest_neighbors(
            pos, [&](T &node, int iterated_index, int center_index) {
              node.add_particle(
                  part, common_weight *
                            // the contribution to density is weighted depending
                            // on what node it is added to
                            (iterated_index == center_index ? big : small));
            });
      } else if (par.smearing() == SmearingMode::Triangular) {
        // unweighted contribution to density
        const double common_weight = dens_factor * prefactor_triangular;
        lat->iterate_in_rectangle(
            pos, triangular_radius, [&](T &node, int ix, int iy, int iz) {
              // compute the position of the node
              const ThreeVector cell_center = lat->cell_center(ix, iy, iz);
              // compute smearing weight
              const double weight_x =
                  triangular_radius[0] - std::abs(cell_center[0] - pos[0]);
              const double weight_y =
                  triangular_radius[1] - std::abs(cell_center[1] - pos[1]);
              const double weight_z =
                  triangular_radius[2] - std::abs(cell_center[2] - pos[2]);
              // add the contribution to the node
              node.add_particle(part,
                                common_weight * weight_x * weight_y * weight_z);
            });
      }
    }  // end of for (const ParticleData &part : particles)
  }    // end of for (const Particles &particles : ensembles)
}

/**
 * Updates the contents on the lattice of DensityOnLattice type.
 *
 * \param[out] lat The lattice of DensityOnLattice type on which the content
 *             will be updated
 * \param[in] old_jmu Auxiliary lattice, filled with current values at t0,
 *            needed for calculating time derivatives
 * \param[in] new_jmu Auxiliary lattice,filled with current values at t0 + dt,
 *            needed for calculating time derivatives
 * \param[in] four_grad_lattice Auxiliary lattice for calculating the
 *            fourgradient of the current
 * \param[in] update Tells if called for update at printout or at timestep
 * \param[in] dens_type Density type to be computed on the lattice
 * \param[in] par a structure containing testparticles number and gaussian
 *            smearing parameters.
 * \param[in] ensembles The particles vector for each ensemble
 * \param[in] time_step Time step used in the simulation
 * \param[in] compute_gradient Whether to compute the gradients
 */
void update_lattice(
    RectangularLattice<DensityOnLattice> *lat,
    RectangularLattice<FourVector> *old_jmu,
    RectangularLattice<FourVector> *new_jmu,
    RectangularLattice<std::array<FourVector, 4>> *four_grad_lattice,
    const LatticeUpdate update, const DensityType dens_type,
    const DensityParameters &par, const std::vector<Particles> &ensembles,
    const double time_step, const bool compute_gradient);
}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_DENSITY_H_
