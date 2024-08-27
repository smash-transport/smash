/*
 *    Copyright (c) 2014-2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#include "smash/alphaclusterednucleus.h"

#include <cmath>
#include <map>
#include <stdexcept>

#include "smash/configuration.h"
#include "smash/constants.h"
#include "smash/fourvector.h"
#include "smash/numerics.h"
#include "smash/random.h"
#include "smash/threevector.h"

namespace smash {

AlphaClusteredNucleus::AlphaClusteredNucleus(
    const std::map<PdgCode, int> &particle_list, int nTest)
    : Nucleus(particle_list, nTest) {}

AlphaClusteredNucleus::AlphaClusteredNucleus(Configuration &config, int nTest,
                                             bool auto_alphaclustering)
    : Nucleus(config, nTest) {
  int A = Nucleus::number_of_particles();
  int Z = Nucleus::number_of_protons();
  if (A == 16 && Z == 8) {
    // Set the Woods-Saxon parameters to those of Helium
    set_diffusiveness(0.322);
    set_nuclear_radius(1.676);
  } else {
    throw std::domain_error(
        "Alpha-Clustering is only implemented for oxygen nuclei. Please, check "
        "the 'Alpha_Clustered' section in your input file.");
  }
  if (!auto_alphaclustering) {
    tetrahedron_sidelength_ = config.take({"Alpha_Clustered", "Sidelength"});
  }
  scale_tetrahedron_vertex_positions(tetrahedron_sidelength_);
  if (config.has_value({"Orientation"})) {
    Configuration sub_conf = config.extract_sub_configuration({"Orientation"});
    set_orientation_from_config(sub_conf);
  }
}

ThreeVector AlphaClusteredNucleus::distribute_nucleon() {
  ThreeVector alpha_clustering_shift =
      tetrahedron_vertex_positions_[tetrahedron_vertex_index_ %
                                    tetrahedron_vertex_positions_.size()];
  tetrahedron_vertex_index_ += 1;
  // Get the solid angle of the nucleon.
  Angles dir;
  dir.distribute_isotropically();
  // diffusiveness_ zero or negative? Use hard sphere.
  if (almost_equal(Nucleus::get_diffusiveness(), 0.)) {
    return dir.threevec() * Nucleus::get_nuclear_radius() *
               std::cbrt(random::canonical()) +
           alpha_clustering_shift;
  }
  if (almost_equal(Nucleus::get_nuclear_radius(), 0.)) {
    return smash::ThreeVector() + alpha_clustering_shift;
  }
  double radius_scaled =
      Nucleus::get_nuclear_radius() / Nucleus::get_diffusiveness();
  double prob_range1 = 1.0;
  double prob_range2 = 3. / radius_scaled;
  double prob_range3 = 2. * prob_range2 / radius_scaled;
  double prob_range4 = 1. * prob_range3 / radius_scaled;
  double ranges234 = prob_range2 + prob_range3 + prob_range4;
  double t;
  /// \li Decide which branch \f$\tilde p^{({\rm I - IV})}\f$ to go into
  do {
    double which_range = random::uniform(-prob_range1, ranges234);
    if (which_range < 0.0) {
      t = radius_scaled * (std::cbrt(random::canonical()) - 1.);
    } else {
      t = -std::log(random::canonical());
      if (which_range >= prob_range2) {
        t -= std::log(random::canonical());
        if (which_range >= prob_range2 + prob_range3) {
          t -= std::log(random::canonical());
        }
      }
    }
    /**
     * \li Generate \f$t\f$ from the distribution in the respective
     * branches
     * \li \a Reject that number with a probability
     * \f$1-(1+\exp(-|t|))^{-1}\f$ (the efficiency of this should be
     * \f$\gg \frac{1}{2}\f$)
     */
  } while (random::canonical() > 1. / (1. + std::exp(-std::abs(t))));
  /// \li Shift and rescale \f$t\f$ to \f$r = d\cdot t + r_0\f$
  double position_scaled = t + radius_scaled;
  double position = position_scaled * Nucleus::get_diffusiveness();
  return dir.threevec() * position + alpha_clustering_shift;
}

void AlphaClusteredNucleus::scale_tetrahedron_vertex_positions(
    double sidelength) {
  for (auto &&i : tetrahedron_vertex_positions_) {
    i = i * std::sqrt(6) / 4 * sidelength;
  }
}

}  // namespace smash
