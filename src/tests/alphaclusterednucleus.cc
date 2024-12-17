/*
 *
 *    Copyright (c) 2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/alphaclusterednucleus.h"

#include <map>
#include <vector>

#include "smash/constants.h"
#include "smash/fourvector.h"
#include "smash/nucleus.h"
#include "smash/particledata.h"
#include "smash/pdgcode.h"
#include "smash/pow.h"

namespace particles_txt {
#include <particles.txt.h>
}

using namespace smash;

static Configuration get_oxygen_configuration() {
  return Configuration{R"(
    Modi:
      Collider:
        Projectile:
          Particles: {2212: 8, 2112: 8}
  )"};
}

TEST(init_particle_types) {
  ParticleType::create_type_list(particles_txt::data);
}

// Check if each cluster is composed of 2 neutrons and 2 protons
TEST(correct_clusters) {
  auto config = get_oxygen_configuration();
  AlphaClusteredNucleus alpha_clustered_oxygen(config, 1, true);
  alpha_clustered_oxygen.set_diffusiveness(0.0);
  alpha_clustered_oxygen.set_nuclear_radius(0.0001);
  alpha_clustered_oxygen.arrange_nucleons();
  // This for loop goes through the first 4 nucleons in the list because they
  // are in the 4 different clusters
  for (auto i = alpha_clustered_oxygen.begin();
       i != alpha_clustered_oxygen.end() - 12; i++) {
    int proton_count = 0;
    int neutron_count = 0;
    for (auto j = alpha_clustered_oxygen.begin();
         j != alpha_clustered_oxygen.end(); j++) {
      if (((i->position().threevec() - j->position().threevec()).abs() <
           0.01) &&
          (j->is_proton())) {
        proton_count += 1;
      } else if (((i->position().threevec() - j->position().threevec()).abs() <
                  0.01) &&
                 (j->is_neutron())) {
        neutron_count += 1;
      }
    }
    COMPARE(neutron_count, 2);
    COMPARE(proton_count, 2);
  }
}

/* NOTE: Since arranging nucleons is not deterministic, here both diffusiveness
 * and the nuclear radius are set to 0, such that the position of the Helium
 * nucleons is the same.
 */
TEST(tetrahedron_scaling) {
  auto config = Configuration{R"(
    Modi:
      Collider:
        Projectile:
          Particles: {2212: 8, 2112: 8}
          Alpha_Clustered:
            Side_Length: 10
  )"};
  AlphaClusteredNucleus alpha_clustered_oxygen(config, 1, false);
  alpha_clustered_oxygen.set_diffusiveness(0.0);
  alpha_clustered_oxygen.set_nuclear_radius(0.0);
  alpha_clustered_oxygen.arrange_nucleons();
  for (auto it = alpha_clustered_oxygen.begin();
       it != alpha_clustered_oxygen.end(); it++) {
    COMPARE_RELATIVE_ERROR(it->position().threevec().abs() * 4 / std::sqrt(6),
                           10.,
                           1.e-12);  // Checks if the side length is 10
  }
}

TEST(scaling_tetrahedron_twice_does_not_change_it) {
  auto config = get_oxygen_configuration();
  AlphaClusteredNucleus alpha_clustered_oxygen(config, 1, true);
  alpha_clustered_oxygen.set_diffusiveness(0.0);
  alpha_clustered_oxygen.set_nuclear_radius(0.0);
  alpha_clustered_oxygen.arrange_nucleons();
  std::vector<ThreeVector> initial_positions{};
  std::transform(alpha_clustered_oxygen.cbegin(), alpha_clustered_oxygen.cend(),
                 std::back_inserter(initial_positions),
                 [](auto in) { return in.position().threevec(); });
  alpha_clustered_oxygen.scale_tetrahedron_vertex_positions(3.42);
  alpha_clustered_oxygen.arrange_nucleons();
  std::vector<ThreeVector> final_positions{};
  std::transform(alpha_clustered_oxygen.cbegin(), alpha_clustered_oxygen.cend(),
                 std::back_inserter(final_positions),
                 [](auto in) { return in.position().threevec(); });
  for (auto it1 = initial_positions.begin(), it2 = final_positions.begin();
       it1 != initial_positions.end(); it1++, it2++) {
    COMPARE_RELATIVE_ERROR(it1->abs(), it2->abs(), 1.e-12);
  }
}
