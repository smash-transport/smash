/*
 *    Copyright (c) 2024-2025
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
        "Alpha-clustering is only implemented for oxygen nuclei. Please, check "
        "the 'Alpha_Clustered' section in your input file.");
  }

  /* NOTE: Here the Nucleus parent class has been initialised and the Projectile
   * or Target subsection could be already be entirely extracted. Therefore
   * unconditionally call here the is_about_projectile function might fail.
   * This has to be called where appropriate and repeating it is not a problem.
   */
  if (!automatic && has_projectile_or_target(config)) {
    const bool is_projectile = is_about_projectile(config);
    const auto &side_length_key = [&is_projectile]() {
      return is_projectile
                 ? InputKeys::modi_collider_projectile_alphaClustered_sideLength
                 : InputKeys::modi_collider_target_alphaClustered_sideLength;
    }();
    tetrahedron_side_length_ = config.take(side_length_key);
  }
  scale_tetrahedron_vertex_positions(tetrahedron_side_length_);
  if (has_projectile_or_target(config)) {
    const bool is_projectile = is_about_projectile(config);
    const auto &orientation_section = [&is_projectile]() {
      return is_projectile ? InputSections::m_c_p_orientation
                           : InputSections::m_c_t_orientation;
    }();
    if (config.has_section(orientation_section)) {
      Configuration sub_conf =
          config.extract_complete_sub_configuration(orientation_section);
      set_orientation_from_config(sub_conf);
    }
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
    position = position / position.abs() * std::sqrt(6) / 4 * side_length;
  }
}

}  // namespace smash
