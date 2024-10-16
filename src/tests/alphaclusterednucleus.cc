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

// Check if the tetrahedron has the side length that is specified
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
  alpha_clustered_oxygen.set_nuclear_radius(0.01);
  alpha_clustered_oxygen.arrange_nucleons();
  for (auto i = alpha_clustered_oxygen.begin();
       i != alpha_clustered_oxygen.end(); i++) {
    COMPARE_RELATIVE_ERROR(i->position().threevec().abs(),
                           std::sqrt(6) / 4 * 10,
                           0.05);  // Checks if the side length is 10
  }
}
