/*
 *    Copyright (c) 2016-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_GRANDCAN_THERMALIZER_H_
#define SRC_INCLUDE_GRANDCAN_THERMALIZER_H_

#include <memory>
#include <vector>

#include "angles.h"
#include "clock.h"
#include "configuration.h"
#include "density.h"
#include "distributions.h"
#include "forwarddeclarations.h"
#include "hadgas_eos.h"
#include "lattice.h"
#include "particledata.h"
#include "quantumnumbers.h"

namespace Smash {

/**
 * The ThermLatticeNode class is intended to compute thermodynamical quantities
 * in a cell given particles. It accumulates the upper row of the
 * energy-momentum tensor \f$ T^{\mu 0}\f$, baryon density nb and strangeness
 * densities in the computational frame. From these quantities it allows to
 * compute the local rest frame quantites: temperature T, chemical potentials
 * mub and mus, the velocity of the local rest frame with respect to
 * computational frame.
 *
 * An example of the intended use is:
 *
 * ThermLatticeNode node;
 * for (const ParticleData &particle: some_particles_list) {
 *   const double some_smearing_factor = ...;
 *   node.add_particle(particle, some_smearing_factor);
 * }
 * HadronGasEos eos = HadronGasEos(false);
 * node.compute_rest_frame_quantities(eos);
 *
 * Note that before calling the compute_rest_frame_quantities T, mu, p and e
 * and simply set to zero.
 */
class ThermLatticeNode {
 public:
  ThermLatticeNode();
  /// Add particle contribution to Tmu0, nb and ns
  void add_particle(const ParticleData& p, double factor);
  /// Compute T, mu, v given Tmu0, nb and ns
  void compute_rest_frame_quantities(HadronGasEos& eos);
  /// Simply set rest frame quantities. Useful for testing.
  void set_rest_frame_quantities(double T0, double mub0, double mus0,
                                 const ThreeVector v0);

  FourVector Tmu0() const { return Tmu0_; }
  double nb() const { return nb_; }
  double ns() const { return ns_; }
  double e() const { return e_; }
  double p() const { return p_; }
  ThreeVector v() const { return v_; }
  double T() const { return T_; }
  double mub() const { return mub_; }
  double mus() const { return mus_; }

 private:
  /// Four-momentum of the cell
  FourVector Tmu0_;
  /// Net baryon density of the cell in the computational frame
  double nb_;
  /// Net strangeness density of the cell in the computational frame
  double ns_;
  /// Energy density in the rest frame
  double e_;
  /// Pressure in the rest frame
  double p_;
  /// Velocity of the rest frame
  ThreeVector v_;
  /// Temperature
  double T_;
  /// Baryon chemical potential
  double mub_;
  /// Strangeness chemical potential
  double mus_;
};

std::ostream& operator<<(std::ostream& s, const ThermLatticeNode& node);

enum class HadronClass {
  /// All baryons
  Baryon = 0,
  /// All anti-baryons
  Antibaryon = 1,
  /// Mesons with strangeness S > 0
  PositiveSMeson = 2,
  /// Mesons with strangeness S < 0
  NegativeSMeson = 3,
  /// Non-strange mesons (S = 0) with electric cherge Q > 0
  PositiveQZeroSMeson = 4,
  /// Non-strange mesons (S = 0) with electric cherge Q < 0
  NegativeQZeroSMeson = 5,
  /// Neutral non-strange mesons
  ZeroQZeroSMeson = 6,
};

/** The GrandCanThermalizer class implements the following functionality:
 *  1. Create a lattice and find the local rest frame energy density in each
 *     cell from the particles.
 *  2. Remove particles from the cells, where the energy density is high enough.
 *     Save the energy, momentum and quantum numbers of the removed particles.
 *  3. Sample new particles instead of the removed ones according to the
 *     grand-canonical thermal distribution, but with an additional constraint:
 *     the energy, momentum and quantum numbers should be the same as those of
 *     the removed particles.
 *
 *  The step 3. is a challenging task, so several algorithms are implemented
 *  that try to fulfil the requirements. The algorithms are a trade-off between
 *  mathematical rigour and computational speed. All of them are shown
 *  to reproduce the mean values of multiplicities correctly. However, this
 *  is not the case for multiplicity fluctuations. For details see
 *  \iref{Oliinychenko:2016vkg}.
 */
class GrandCanThermalizer {
 public:
  /// Create the thermalizer: allocate the lattice
  GrandCanThermalizer(const std::array<double, 3> lat_sizes,
                      const std::array<int, 3> n_cells,
                      const std::array<double, 3> origin, bool periodicity,
                      double e_critical, double t_start, double delta_t,
                      ThermalizationAlgorithm algo);
  GrandCanThermalizer(Configuration& conf,
                      const std::array<double, 3> lat_sizes,
                      const std::array<double, 3> origin, bool periodicity)
      : GrandCanThermalizer(
            lat_sizes, conf.take({"Cell_Number"}), origin, periodicity,
            conf.take({"Critical_Edens"}), conf.take({"Start_Time"}),
            conf.take({"Timestep"}),
            conf.take({"Algorithm"}, ThermalizationAlgorithm::BiasedBF)) {}
  /// Check that the clock is close to n * period of thermalization
  bool is_time_to_thermalize(const Clock& clock) const {
    const double t = clock.current_time();
    const int n = static_cast<int>(std::floor((t - t_start_) / period_));
    return (t > t_start_ &&
            t < t_start_ + n * period_ + clock.timestep_duration());
  }
  /// Compute all the thermodynamical quantities on the lattice from particles.
  void update_lattice(const Particles& particles, const DensityParameters& par,
                      bool ignore_cells_under_treshold = true);
  /// Simply returns a vector uniformly sampled from the rectangular cell.
  ThreeVector uniform_in_cell() const;
  /** Changes energy and momenta of the particles in plist to match the
   *  required_total_momentum. The procedure is described in
   *  \iref{Oliinychenko:2016vkg}.
   */
  void renormalize_momenta(ParticleList& plist,
                           const FourVector required_total_momentum);

  // Functions for BF-sampling algorithm

  /**
   *  The sample_multinomial function samples integer numbers n_i distributed
   *  according to the multinomial distribution with sum N: \f$ p(n_1, n_2,
   *  \dots) = \prod a_i^{n_i} \times \frac{N!}{n_1!n_2! \dots} \f$ if \f$ \sum
   *  n_i = N \f$  and \f$ p = 0 \f$ otherwise.
   *
   * The array mult_sort_ contains real numbers \f$ a_i \f$. The numbers \f$
   * n_i \f$ are saved in the mult_int_ array. Only particles of class
   * particle_class are sampled, where particle_class is defined by the
   * get_class function.
   */
  void sample_multinomial(HadronClass particle_class, int N);
  /**
   * The total number of particles of sort type_index is defined by mult_int_
   * array. This function samples mult_int_[type_index] particles. It chooses
   * the cell to sample, picks up momentum and coordinate from the
   * corresponding distributions.
   */
  void sample_in_random_cell_BF_algo(ParticleList& plist, const double time,
                                     size_t type_index);
  /**
   * Samples particles to the sampled_list according to the BF algorithm.
   * Quantum numbers of the sampled particles are required to be as in
   * conserved_initial.
   */
  void thermalize_BF_algo(ParticleList& sampled_list,
                          QuantumNumbers& conserved_initial, double time,
                          int ntest);

  // Functions for mode-sampling algorithm

  /// Computes average number of particles in each cell.
  template <typename F>
  void compute_N_in_cells_mode_algo(F&& condition) {
    N_in_cells_.clear();
    N_total_in_cells_ = 0.0;
    for (auto cell_index : cells_to_sample_) {
      const ThermLatticeNode cell = (*lat_)[cell_index];
      const double gamma = 1.0 / std::sqrt(1.0 - cell.v().sqr());
      double N_tot = 0.0;
      for (ParticleTypePtr i : eos_typelist_) {
        if (condition(i->strangeness(), i->baryon_number(), i->charge())) {
          // N_i = n u^mu dsigma_mu = (isochronous hypersurface) n * V * gamma
          N_tot +=
              cell_volume_ * gamma * HadronGasEos::partial_density(
                                         *i, cell.T(), cell.mub(), cell.mus());
        }
      }
      N_in_cells_.push_back(N_tot);
      N_total_in_cells_ += N_tot;
    }
  }

  /** Samples one particle. The species, cell, momentum and coordinate
   *  are chosen from the corresponding distributions. The condition
   *  function limits the choice of possible species.
   *
   *  Condition is the function of the signature
   *  bool condition(int strangeness, int baryon_number, int charge);
   */
  template <typename F>
  ParticleData sample_in_random_cell_mode_algo(const double time,
                                               F&& condition) {
    // Choose random cell, probability = N_in_cell/N_total
    double r = Random::uniform(0.0, N_total_in_cells_);
    double partial_sum = 0.0;
    int index_only_thermalized = -1;
    while (partial_sum < r) {
      index_only_thermalized++;
      partial_sum += N_in_cells_[index_only_thermalized];
    }
    const int cell_index = cells_to_sample_[index_only_thermalized];
    const ThermLatticeNode cell = (*lat_)[cell_index];
    const ThreeVector cell_center = lat_->cell_center(cell_index);
    const double gamma = 1.0 / std::sqrt(1.0 - cell.v().sqr());
    const double N_in_cell = N_in_cells_[index_only_thermalized];

    // Which sort to sample - probability N_i/N_tot
    r = Random::uniform(0.0, N_in_cell);
    double N_sum = 0.0;
    ParticleTypePtr type_to_sample;
    for (ParticleTypePtr i : eos_typelist_) {
      if (!condition(i->strangeness(), i->baryon_number(), i->charge())) {
        continue;
      }
      N_sum += cell_volume_ * gamma * HadronGasEos::partial_density(
                                          *i, cell.T(), cell.mub(), cell.mus());
      if (N_sum >= r) {
        type_to_sample = i;
        break;
      }
    }

    ParticleData particle(*type_to_sample);
    // Note: it's pole mass for resonances!
    const double m = type_to_sample->mass();
    // Position
    particle.set_4position(FourVector(time, cell_center + uniform_in_cell()));
    // Momentum
    double momentum_radial = sample_momenta_from_thermal(cell.T(), m);
    Angles phitheta;
    phitheta.distribute_isotropically();
    particle.set_4momentum(m, phitheta.threevec() * momentum_radial);
    particle.boost_momentum(-cell.v());

    return particle;
  }

  /*
   * Samples particles to the sampled_list according to the mode algorithm.
   * Quantum numbers of the sampled particles are required to be as in
   * conserved_initial.
   */
  void thermalize_mode_algo(ParticleList& sampled_list,
                            QuantumNumbers& conserved_initial, double time);

  /// Main thermalize function, that chooses algorithm
  void thermalize(Particles& particles, double time, int ntest);

  void print_statistics(const Clock& clock) const;

  RectangularLattice<ThermLatticeNode>& lattice() const { return *lat_; }
  double e_crit() const { return e_crit_; }

 private:
  ParticleTypePtrList list_eos_particles() const {
    ParticleTypePtrList res;
    for (const ParticleType& ptype : ParticleType::list_all()) {
      if (HadronGasEos::is_eos_particle(ptype)) {
        res.push_back(&ptype);
      }
    }
    return res;
  }
  HadronClass get_class(size_t typelist_index) const {
    const int B = eos_typelist_[typelist_index]->baryon_number();
    const int S = eos_typelist_[typelist_index]->strangeness();
    const int ch = eos_typelist_[typelist_index]->charge();
    // clang-format off
    return (B > 0) ? HadronClass::Baryon :
           (B < 0) ? HadronClass::Antibaryon :
           (S > 0) ? HadronClass::PositiveSMeson :
           (S < 0) ? HadronClass::NegativeSMeson :
           (ch > 0) ? HadronClass::PositiveQZeroSMeson :
           (ch < 0) ? HadronClass::NegativeQZeroSMeson :
                      HadronClass::ZeroQZeroSMeson;
    // clang-format on
  }
  /// Returns multiplicity of the hadron class cl
  double mult_class(const HadronClass cl) const {
    return mult_classes_[static_cast<size_t>(cl)];
  }
  std::vector<double> N_in_cells_;
  std::vector<size_t> cells_to_sample_;
  HadronGasEos eos_ = HadronGasEos(true);
  std::unique_ptr<RectangularLattice<ThermLatticeNode>> lat_;
  const ParticleTypePtrList eos_typelist_;
  const size_t N_sorts_;
  std::vector<double> mult_sort_;
  std::vector<int> mult_int_;
  std::array<double, 7> mult_classes_;
  double N_total_in_cells_;
  double cell_volume_;
  const double e_crit_;
  const double t_start_;
  const double period_;
  const ThermalizationAlgorithm algorithm_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_GRANDCAN_THERMALIZER_H_
