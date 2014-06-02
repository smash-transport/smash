/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "../include/particledata.h"
#include "../include/pdgcode.h"

using namespace Smash;

TEST(init_from_pdgcode) {
  PdgCode pdg = 0x211;
  ParticleData p{pdg};

  COMPARE(p.id(), -1);
  COMPARE(p.pdgcode(), pdg);
  COMPARE(p.id_process(), -1);
  COMPARE(p.collision_time(), 0.);
}
