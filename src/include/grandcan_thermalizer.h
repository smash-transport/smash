/*
 *    Copyright (c) 2016-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_GRANDCAN_THERMALIZER_H_
#define SRC_INCLUDE_GRANDCAN_THERMALIZER_H_

#include "clock.h"
#include "configuration.h"
#include "density.h"
#include "forwarddeclarations.h"
#include "hadgas_eos.h"
#include "lattice.h"
#include "particledata.h"
#include "quantumnumbers.h"

namespace Smash {

class ThermLatticeNode {
 public:
  ThermLatticeNode();
  void add_particle(const ParticleData& p, double factor);
  void compute_rest_frame_quantities(HadronGasEos& eos);

  FourVector Tmu0() const { return Tmu0_; }
  double nb() const { return nb_; }
  double ns() const { return ns_; }
  double e()  const { return e_; }
  double p()  const { return p_; }
  ThreeVector v()  const { return v_; }
  double T()   const { return T_; }
  double mub() const { return mub_; }
  double mus() const { return mus_; }
 private:
  FourVector Tmu0_;
  double nb_;
  double ns_;
  double e_;
  double p_;
  ThreeVector v_;
  double T_;
  double mub_;
  double mus_;
};

std::ostream &operator<<(std::ostream &s, const ThermLatticeNode &node);

class GrandCanThermalizer {
 public:
  /// Create the thermalizer: allocate the lattice
  GrandCanThermalizer(const std::array<float, 3> lat_sizes,
                      const std::array<int, 3> n_cells,
                      const std::array<float, 3> origin,
                      bool periodicity,
                      float e_critical,
                      float t_start,
                      float delta_t,
                      ThermalizationAlgorithm algo);
  GrandCanThermalizer(Configuration& conf,
                      const std::array<float, 3> lat_sizes,
                      const std::array<float, 3> origin,
                      bool periodicity) :
    GrandCanThermalizer(lat_sizes,
                        conf.take({"Cell_Number"}),
                        origin,
                        periodicity,
                        conf.take({"Critical_Edens"}),
                        conf.take({"Start_Time"}),
                        conf.take({"Timestep"}),
                        conf.take({"Algorithm"},
                            ThermalizationAlgorithm::BiasedBF)) {};
  bool is_time_to_thermalize(const Clock& clock) const {
    const float t = clock.current_time();
    const int n = static_cast<int>(std::floor((t - t_start_)/period_));
    return (t > t_start_ &&
            t < t_start_ + n*period_ + clock.timestep_duration());
  }

  void update_lattice(const Particles& particles,
                      const DensityParameters& par,
                      bool ignore_cells_under_treshold = true);
  ThreeVector uniform_in_cell() const;
  void sample_multinomial(int particle_class, int N);
  void renormalize_momenta(ParticleList& plist,
           const FourVector required_total_momentum);

  // Functions for BF-sampling algorithm
  void sample_in_random_cell_BF_algo(ParticleList& plist,
                                     const double time,
                                     size_t type_index);
  void thermalize_BF_algo(Particles& particles,
                          double time, int ntest);

  // Functions for mode-sampling algorithm
  void compute_N_in_cells_mode_algo(
                        std::function<bool(int, int, int)> condition);
  ParticleData sample_in_random_cell_mode_algo(const double time,
                        std::function<bool(int, int, int)> condition);
  void thermalize_mode_algo(Particles& particles, double time);

  void thermalize(Particles& particles,
                  double time, int ntest) {
    if (algorithm_ == ThermalizationAlgorithm::BiasedBF ||
        algorithm_ == ThermalizationAlgorithm::UnbiasedBF) {
      thermalize_BF_algo(particles, time, ntest);
    } else if (algorithm_ == ThermalizationAlgorithm::ModeSampling) {
      thermalize_mode_algo(particles, time);
    } else {
      throw std::invalid_argument("This thermalization algorithm is"
                                  " not yet implemented");
    }
  }


  void print_statistics(const Clock& clock) const;

  RectangularLattice<ThermLatticeNode>& lattice() const {
    return *lat_;
  }
  float e_crit() const { return e_crit_; }

 private:
  ParticleTypePtrList list_eos_particles() const {
    ParticleTypePtrList res;
    for (const ParticleType &ptype : ParticleType::list_all()) {
      if (HadronGasEos::is_eos_particle(ptype)) {
        res.push_back(&ptype);
      }
    }
    return res;
  }
  int get_class(size_t typelist_index) const {
    const int B = eos_typelist_[typelist_index]->baryon_number();
    const int S = eos_typelist_[typelist_index]->strangeness();
    const int ch = eos_typelist_[typelist_index]->charge();
    return (B > 0) ? 0 :
           (B < 0) ? 1 :
           (S > 0) ? 2 :
           (S < 0) ? 3 :
           (ch > 0) ? 4 :
           (ch < 0) ? 5 : 6;
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
  float cell_volume_;
  const float e_crit_;
  const float t_start_;
  const float period_;
  const ThermalizationAlgorithm algorithm_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_GRANDCAN_THERMALIZER_H_
