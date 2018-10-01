/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "../include/smash/constants.h"
#include "../include/smash/deformednucleus.h"
#include "../include/smash/fourvector.h"
#include "../include/smash/nucleus.h"
#include "../include/smash/particledata.h"
#include "../include/smash/pdgcode.h"

#include <map>
#include <vector>

namespace particles_txt {
#include <particles.txt.h>
}

using namespace smash;

std::map<PdgCode, int> small_list = {{0x2212, 1}};
std::map<PdgCode, int> large_list = {{0x2212, 92}, {0x2112, 146}};

TEST(init_particle_types) {
  ParticleType::create_type_list(particles_txt::data);
}

TEST(rotate_phi) {
  DeformedNucleus dnucleus(small_list, 1);
  // Plan is to rotate the (0, 1, 0, 1) vector by pi/2.
  // Rotation by pi/2 means (0, 1, 0, 1) -> (0, 0, 1, 1)
  dnucleus.set_azimuthal_angle(M_PI / 2);
  FourVector expectation = FourVector(0., 0., 1., 1.);
  for (auto i = dnucleus.begin(); i != dnucleus.end(); i++) {
    i->set_4position(FourVector(0., 1., 0., 1.));
  }
  dnucleus.rotate();
  FourVector actual;
  for (auto i = dnucleus.begin(); i != dnucleus.end(); i++) {
    actual = i->position();
  }
  COMPARE_ABSOLUTE_ERROR(actual.x0(), expectation.x0(), 1e-7);
  COMPARE_ABSOLUTE_ERROR(actual.x1(), expectation.x1(), 1e-7);
  COMPARE_ABSOLUTE_ERROR(actual.x2(), expectation.x2(), 1e-7);
  COMPARE_ABSOLUTE_ERROR(actual.x3(), expectation.x3(), 1e-7);
}

// TEST(rotate_theta) {
//   // same as rotate_phi
// }

// TEST(rotate_both) {
//   // same as rotate_phi
// }

// TEST(dws_rcut) {
//   // similar to nucleus test
// }

// TEST(dws_thetacut) {
//   // similar to nucleus test
// }
