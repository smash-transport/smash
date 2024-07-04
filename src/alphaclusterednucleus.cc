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
  const auto &orientation_section = [&config]() {
    return is_configuration_about_projectile(config)
               ? InputSections::m_c_p_orientation
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
    double sidelength) {
  for (auto &&i : tetrahedron_vertex_positions_) {
    i = i * std::sqrt(6) / 4 * sidelength;
  }
}

}  // namespace smash
