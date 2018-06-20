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

namespace smash {

/**
 * The ThermLatticeNode class is intended to compute thermodynamical quantities
 * in a cell given a set of particles. It accumulates the upper row of the
 * energy-momentum tensor \f$ T^{\mu 0}\f$, net baryon density nb and net
 * strangeness densities in the computational frame. From these quantities
 * it allows to compute the local rest frame quantites: temperature T,
 * chemical potentials mub and mus, the velocity of the local rest frame with
 * respect to the computational frame.
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
 * are set to zero.
 */
class ThermLatticeNode {
 public:
  /**
   * Default constructor of thermal quantities on the lattice returning
   * thermodynamic quantities in computational frame
   * \return Tmu0_ Four vector \f$T^{\mu 0}\f$, nb_ Net baryon density at
   * this location, ns_ Net strangeness density, e_ Energy density,
   * v_ 3 vector velocity of local rest frame, T_ Temperature
   * mub_ Net baryon chemical potential, mus_ Net strangeness chemical potential
   */
  ThermLatticeNode();
  /**
   *  Add particle contribution to Tmu0, nb and ns
   *  May look like unused at first glance, but it is actually used
   *  by update_general_lattice, where the node type of the lattice
   *  is templated.
   */
  void add_particle(const ParticleData& p, double factor);
  /**
   * Temperature, chemical potentials and rest frame velocity are
   * calculated given the hadron gas equation of state object
   * \param[in] eos \see HadronGasEos based on Tmu0, nb and ns
   * \return Temperature T, net baryon chemical potential mub,
   * net strangeness chemical potential mus and the velocity of the
   * Landau rest frame, under assumption that energy-momentum tensor
   * has an ideal-fluid form. For more details and discussion see
   * \iref{Oliinychenko:2015lva}. The advantage of this rest frame
   * transformation is that it conserves energy and momentum, even
   * though the dissipative part of the energy-momentum tensor is neglected.
   */
  void compute_rest_frame_quantities(HadronGasEos& eos);
  /**
   * Set all the rest frame quantities to some values, this is useful
   * for testing.
   * \param[out] T0 Rest frame temperature
   * \param[out] mub0 Rest frame net baryon chemical potential
   * \param[out] mus0 Rest frame net strangeness chemical potential
   * \param[out] v0 Velocity of the rest frame
   */
  void set_rest_frame_quantities(double T0, double mub0, double mus0,
                                 const ThreeVector v0);
  /// Get Four-momentum flow of the cell
  FourVector Tmu0() const { return Tmu0_; }
  /// Get net baryon density of the cell in the computational frame
  double nb() const { return nb_; }
  /// Get net strangeness density of the cell in the computational frame
  double ns() const { return ns_; }
  /// Get energy density in the rest frame
  double e() const { return e_; }
  /// Get pressure in the rest frame
  double p() const { return p_; }
  /// Get 3-velocity of the rest frame
  ThreeVector v() const { return v_; }
  /// Get the temperature
  double T() const { return T_; }
  /// Get the net baryon chemical potential
  double mub() const { return mub_; }
  /// Get the net strangeness chemical potential
  double mus() const { return mus_; }

 private:
  /// Four-momentum flow of the cell
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
  /// Net baryon chemical potential
  double mub_;
  /// Net strangeness chemical potential
  double mus_;
};

/**
 * This operator writes all the thermodynamic quantities at a certain
 * position to the file out
 * \param[in] s location of the output
 * \param[in] node position on the lattice, where the output is generated
 */
std::ostream& operator<<(std::ostream& s, const ThermLatticeNode& node);

/**
 * Specifier to classify the different hadron species according to
 * their quantum numbers
 */
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

/*!\Userguide
 * \page input_forced_thermalization_ Forced Thermalization
 *
 * \key Cell_Number (list of 3 doubles, required, no default): \n
 * Number of cells in each direction (x,y,z).
 *
 * \key Critical_Edens (double, required, )
 * Critical energydensity above which forced thermalization is applied, in
 * GeV/fm^3.
 *
 * \key Start_Time (double, required, no default): \n
 * Time after which forced thermalization may be applied (in fm/c), if
 * energydensity is sufficiently high.
 *
 * \key Timestep (double, required, no default): \n
 * Timestep of thermalization, in fm/c.
 *
 * \key Algorithm (string, optional, default = "biased BF") \n
 * Algorithm applied to enforce thermalization. See \iref{Oliinychenko:2016vkg}
 * for more details.
 * \li \key "unbiased BF" - slowest, but theoretically most robust
 * \li \key "biased BF" - faster, but theoretically less robust
 * \li \key "mode sampling" - fastest, but least robust
 *
 */

/**
 * The GrandCanThermalizer class implements the following functionality:
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
  /**
   * Default constructor for the GranCanThermalizer to allocate the lattice
   * \param[in] lat_sizes Size of lattice in x,y and z-direction in fm.
   * \param[in] n_cells Number of cells in x, y and z-direction.
   * \param[in] origin Coordinates of the left, down, near corner of
   * the lattice in fm.
   * \param[in] periodicity Boolean to decide, if the lattice is periodically
   * extended to infinity or not
   * \param[in] e_critical Critical energy density above which the cells are
   * thermalized
   * \param[in] t_start Starting time of the simulation
   * \param[in] delta_t Timestep of the simulation
   * \param[in] algo Choice of algorithm for the canonical sampling
   */
  GrandCanThermalizer(const std::array<double, 3> lat_sizes,
                      const std::array<int, 3> n_cells,
                      const std::array<double, 3> origin, bool periodicity,
                      double e_critical, double t_start, double delta_t,
                      ThermalizationAlgorithm algo);
  /// \see GrandCanThermalizer Exactly the same but taking values from config
  GrandCanThermalizer(Configuration& conf,
                      const std::array<double, 3> lat_sizes,
                      const std::array<double, 3> origin, bool periodicity)
      : GrandCanThermalizer(
            lat_sizes, conf.take({"Cell_Number"}), origin, periodicity,
            conf.take({"Critical_Edens"}), conf.take({"Start_Time"}),
            conf.take({"Timestep"}),
            conf.take({"Algorithm"}, ThermalizationAlgorithm::BiasedBF)) {}
  /**
   * Check that the clock is close to n * period of thermalization, since
   * the thermalization only happens at these times
   * \param[in] clock Current system time
   */
  bool is_time_to_thermalize(const Clock& clock) const {
    const double t = clock.current_time();
    const int n = static_cast<int>(std::floor((t - t_start_) / period_));
    return (t > t_start_ &&
            t < t_start_ + n * period_ + clock.timestep_duration());
  }
  /**
   * Compute all the thermodynamical quantities on the lattice from particles.
   * \param[in] particles Current list of particles \see Particles
   * \param[in] par Parameters necessary for density determination
   * \see DensityParameters
   * \param[in] ignore_cells_under_threshold Boolean that is true by default
   */
  void update_lattice(const Particles& particles, const DensityParameters& par,
                      bool ignore_cells_under_threshold = true);
  /// \return 3 vector uniformly sampled from the rectangular cell.
  ThreeVector uniform_in_cell() const;
  /**
   * Changes energy and momenta of the particles in plist to match the
   *  required_total_momentum. The procedure is described in
   *  \iref{Oliinychenko:2016vkg}.
   * \param[in] plist List of particles \see ParticleList
   * \param[in] required_total_momentum The necessary total momentum of the cell
   */
  void renormalize_momenta(ParticleList& plist,
                           const FourVector required_total_momentum);

  // Functions for BF-sampling algorithm

  /**
   * The sample_multinomial function samples integer numbers n_i distributed
   * according to the multinomial distribution with sum N: \f$ p(n_1, n_2,
   * \dots) = \prod a_i^{n_i} \times \frac{N!}{n_1!n_2! \dots} \f$ if \f$ \sum
   * n_i = N \f$  and \f$ p = 0 \f$ otherwise.
   * \param[in] particle_class A certain group of hadron species \see
   * HadronClass \param[out] N Number of particles to be sampled
   */
  void sample_multinomial(HadronClass particle_class, int N);
  /**
   * The total number of particles of species type_index is defined by mult_int_
   * array that is returned by \see sample_multinomial.
   * This function samples mult_int_[type_index] particles. It chooses
   * randomly the cell to sample and picks up momentum and coordinate from the
   * corresponding distributions.
   * \param[out] plist \see ParticleList of newly produced particles
   * \param[in] time Current time in the simulation to become zero component of
   * sampled particles
   * \param[in] type_index Species that should be sampled
   */
  void sample_in_random_cell_BF_algo(ParticleList& plist, const double time,
                                     size_t type_index);
  /**
   * Samples particles according to the BF algorithm by making use of the
   * \see sample_in_random_cell_BF_algo.
   * Quantum numbers of the sampled particles are required to be equal to the
   * original particles in this region.
   * \param[in] conserved_initial The quantum numbers of the total ensemble of
   * of particles in the region to be thermalized
   * \param[in] time Current time of the simulation
   * \param[in] ntest Number of testparticles
   * \return Particle list with newly sampled particles according to
   * Becattini-Feroni algorithm
   */
  void thermalize_BF_algo(QuantumNumbers& conserved_initial, double time,
                          int ntest);

  // Functions for mode-sampling algorithm

  /**
   * Computes average number of particles in each cell for the mode algorithm.
   * \param[in] condition Specifies the current mode (1 to 7)
   */
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
          N_tot += cell_volume_ * gamma *
                   HadronGasEos::partial_density(*i, cell.T(), cell.mub(),
                                                 cell.mus());
        }
      }
      N_in_cells_.push_back(N_tot);
      N_total_in_cells_ += N_tot;
    }
  }

  /**
   * Samples one particle and the species, cell, momentum and coordinate
   * are chosen from the corresponding distributions. The condition
   * function limits the choice of possible species.
   *
   * Condition is a function of the signature of quantum number S, B and Q.
   * bool condition(int strangeness, int baryon_number, int charge);
   * \param[in] time Current time in simulation
   * \param[in] condition Specifies the actual mode (1 to 7)
   */
  template <typename F>
  ParticleData sample_in_random_cell_mode_algo(const double time,
                                               F&& condition) {
    // Choose random cell, probability = N_in_cell/N_total
    double r = random::uniform(0.0, N_total_in_cells_);
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
    r = random::uniform(0.0, N_in_cell);
    double N_sum = 0.0;
    ParticleTypePtr type_to_sample;
    for (ParticleTypePtr i : eos_typelist_) {
      if (!condition(i->strangeness(), i->baryon_number(), i->charge())) {
        continue;
      }
      N_sum +=
          cell_volume_ * gamma *
          HadronGasEos::partial_density(*i, cell.T(), cell.mub(), cell.mus());
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
    particle.set_formation_time(time);
    return particle;
  }

  /**
   * Samples particles to the according to the mode algorithm.
   * Quantum numbers of the sampled particles are required to be as in
   * conserved_initial.
   * \param[in] conserved_initial Quantum numbers of the original particles
   * in the region to be thermalized
   * \param[in] time Current time of the simulation
   */
  void thermalize_mode_algo(QuantumNumbers& conserved_initial, double time);
  /**
   * Main thermalize function, that chooses the algorithm to follow
   * (BF or mode sampling).
   * \param[out] particles List of sampled particles in thermalized region
   * \param[in] time Current time of the simulation
   * \param[in] ntest number of testparticles
   */
  void thermalize(const Particles& particles, double time, int ntest);

  /**
   * Generates standard output with information about the thermodynamic
   * properties of the lattice, the thermalized region and the volume to
   * be thermalized above the critical energy density
   * \param[in] clock Current time of the simulation
   */
  void print_statistics(const Clock& clock) const;
  /// Getter function for the lattice
  RectangularLattice<ThermLatticeNode>& lattice() const { return *lat_; }
  /// Get the critical energy density
  double e_crit() const { return e_crit_; }
  /// List of particles to be removed from the simulation
  ParticleList particles_to_remove() const { return to_remove_; }
  /// List of newly created particles to be inserted in the simulation
  ParticleList particles_to_insert() const { return sampled_list_; }

 private:
  /**
   * Extracts the particles in the hadron gas equation of state from the
   * complete list of particle types in SMASH
   */
  ParticleTypePtrList list_eos_particles() const {
    ParticleTypePtrList res;
    for (const ParticleType& ptype : ParticleType::list_all()) {
      if (HadronGasEos::is_eos_particle(ptype)) {
        res.push_back(&ptype);
      }
    }
    return res;
  }
  /**
   * Defines the class of hadrons by quantum numbers
   * \param[in] typelist_index Index for a certain quantum number
   */
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
  /// \param[out] cl Multiplicity of the hadron class
  double mult_class(const HadronClass cl) const {
    return mult_classes_[static_cast<size_t>(cl)];
  }
  /// Number of particles to be sampled in one cell
  std::vector<double> N_in_cells_;
  /// Cells above critical energy density
  std::vector<size_t> cells_to_sample_;
  /// Hadron gas equation of state
  HadronGasEos eos_ = HadronGasEos(true);
  /// The lattice on which the thermodynamic quantities are calculated
  std::unique_ptr<RectangularLattice<ThermLatticeNode>> lat_;
  /// Particles to be removed after this thermalization step
  ParticleList to_remove_;
  /// Newly generated particles by thermalizer
  ParticleList sampled_list_;
  /**
   * List of particle types from which equation of state is computed.
   * Most particles are included, but not all of them.
   * For example, photons and leptons are not included. Heavy hadrons, that
   * can originate from Pythia, but do not interact in SMASH are not included.
   * The latter are excluded to avoid violations of charm and bottomness
   * conservation, when HadronGasEoS is used for forced thermalization.
   */
  const ParticleTypePtrList eos_typelist_;
  /// Number of different species to be sampled
  const size_t N_sorts_;
  /// Real number multiplicity for each particle type
  std::vector<double> mult_sort_;
  /// Integer multiplicity for each particle type
  std::vector<int> mult_int_;
  /**
   * The different hadron species according to the enum defined
   * in \see HadronClass
   */
  std::array<double, 7> mult_classes_;
  /// Total number of particles over all cells in thermalization region
  double N_total_in_cells_;
  /**
   * Volume of a single cell, necessary to convert thermal densities to
   * actual particle numbers
   */
  double cell_volume_;
  /// Critical energy density above which cells are thermalized
  const double e_crit_;
  /// Starting time of the simulation
  const double t_start_;
  /// Defines periodicity of the lattice in fm
  const double period_;
  /// Algorithm to choose for sampling of particles obeying conservation laws
  const ThermalizationAlgorithm algorithm_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_GRANDCAN_THERMALIZER_H_
