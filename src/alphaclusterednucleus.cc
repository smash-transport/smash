/*
 *    Copyright (c) 2024
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
#include "smash/input_keys.h"
#include "smash/numerics.h"
#include "smash/random.h"
#include "smash/threevector.h"

namespace smash {

AlphaClusteredNucleus::AlphaClusteredNucleus(
    const std::map<PdgCode, int> &particle_list, int n_test)
    : Nucleus(particle_list, n_test) {}

AlphaClusteredNucleus::AlphaClusteredNucleus(Configuration &config, int n_test,
                                             bool automatic)
    : Nucleus(config, n_test) {
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
  const bool is_projectile = is_about_projectile(config);

  if (!automatic) {
    const auto &side_length_key = [&is_projectile]() {
      return is_projectile
                 ? InputKeys::modi_collider_projectile_alphaClustered_sideLength
                 : InputKeys::modi_collider_target_alphaClustered_sideLength;
    }();
    tetrahedron_side_length_ = config.take(side_length_key);
  }
  scale_tetrahedron_vertex_positions(tetrahedron_side_length_);
  const auto &orientation_section = [&is_projectile]() {
    return is_projectile ? InputSections::m_c_p_orientation
                         : InputSections::m_c_t_orientation;
  }();
  if (config.has_section(orientation_section)) {
    Configuration sub_conf = config.extract_sub_configuration({"Orientation"});
    set_orientation_from_config(sub_conf);
  }
}

ThreeVector AlphaClusteredNucleus::distribute_nucleon() {
  ThreeVector alpha_clustering_shift =
      tetrahedron_vertex_positions_[tetrahedron_vertex_index_ %
                                    tetrahedron_vertex_positions_.size()];
  tetrahedron_vertex_index_ += 1;

  return Nucleus::distribute_nucleon() + alpha_clustering_shift;
}

void AlphaClusteredNucleus::scale_tetrahedron_vertex_positions(
    double side_length) {
  for (auto &position : tetrahedron_vertex_positions_) {
    position = position * std::sqrt(6) / 4 * side_length;
  }
}

}  // namespace smash
