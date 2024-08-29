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

std::map<PdgCode, int> oxygen_list = {{0x2212, 8},   // protons
                                      {0x2112, 8}};  // neutrons

TEST(init_particle_types) {
  ParticleType::create_type_list(particles_txt::data);
}

// Check if each cluster is composed of 2 neutrons and 2 protons
TEST(correct_clusters) {
  AlphaClusteredNucleus alphaclusteredoxygen(oxygen_list, 1);
  alphaclusteredoxygen.set_diffusiveness(0.0);
  alphaclusteredoxygen.set_nuclear_radius(0.0001);
  alphaclusteredoxygen.arrange_nucleons();
  // This for loop goes through the first 4 nucleons in the list because they
  // are in the 4 different clusters
  for (auto i = alphaclusteredoxygen.begin();
       i != alphaclusteredoxygen.end() - 12; i++) {
    int proton_count = 0;
    int neutron_count = 0;
    for (auto j = alphaclusteredoxygen.begin(); j != alphaclusteredoxygen.end();
         j++) {
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

// Check if the tetrahedron has the sidelength that is specified
TEST(tetrahedron_scaling) {
  AlphaClusteredNucleus alphaclusteredoxygen(oxygen_list, 1);
  alphaclusteredoxygen.set_diffusiveness(0.0);
  alphaclusteredoxygen.set_nuclear_radius(0.01);
  alphaclusteredoxygen.scale_tetrahedron_vertex_positions(
      5);  // Sets the sidelength to 5
  alphaclusteredoxygen.arrange_nucleons();
  for (auto i = alphaclusteredoxygen.begin(); i != alphaclusteredoxygen.end();
       i++) {
    COMPARE_RELATIVE_ERROR(i->position().threevec().abs(), std::sqrt(6) / 4 * 5,
                           0.1);  // Checks if the Sidelength is 5
  }
}