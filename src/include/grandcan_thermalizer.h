/*
 *    Copyright (c) 2016-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_GRANDCAN_THERMALIZER_H_
#define SRC_INCLUDE_GRANDCAN_THERMALIZER_H_

#include "hadgas_eos.h"

#include "clock.h"
#include "configuration.h"
#include "density.h"
#include "particledata.h"
#include "lattice.h"
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
  GrandCanThermalizer(const std::array<float, 3> cell_sizes,
                      const std::array<int, 3> n_cells,
                      float e_critical,
                      float t_start,
                      float delta_t);
  GrandCanThermalizer(Configuration& conf) :
    GrandCanThermalizer(conf.take({"Lattice_Spacing"}),
                        conf.take({"Cell_Number"}),
                        conf.take({"Critical_Edens"}),
                        conf.take({"Start_Time"}),
                        conf.take({"Timestep"})) {};
  bool is_time_to_thermalize(const Clock& clock) const {
    const float t = clock.current_time();
    const int n = static_cast<int>(std::floor((t - t_start_)/period_));
    return (t > t_start_ &&
            t < t_start_ + n*period_ + clock.timestep_duration());
  }
  void update_lattice(const Particles& particles, const DensityParameters& par,
                      bool ignore_cells_under_treshold = true);
  ThreeVector uniform_in_cell() const;
  void sample_in_random_cell(ParticleList& plist, const double time);
  void thermalize(Particles& particles, double time);
  void print_statistics() const;

  RectangularLattice<ThermLatticeNode>& lattice() const {
    return *lat_;
  }
  float e_crit() const { return e_crit_; }

 private:
  std::vector<size_t> cells_to_sample_;
  HadronGasEos eos_ = HadronGasEos(true);
  std::unique_ptr<RectangularLattice<ThermLatticeNode>> lat_;
  float cell_volume_;
  const float e_crit_;
  const float t_start_;
  const float period_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_GRANDCAN_THERMALIZER_H_
