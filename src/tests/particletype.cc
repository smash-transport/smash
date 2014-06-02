/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "../include/particletype.h"

using namespace Smash;

TEST(assign_default) {
  ParticleType A;
  COMPARE(A.name(), "unknown");
  COMPARE(A.mass(), -1.f);
  COMPARE(A.width(), -1.f);
  COMPARE(A.pdgcode(), PdgCode::invalid());
  COMPARE(A.isospin(), 0);
  COMPARE(A.charge(), 0);
}

TEST(assign) {
  PdgCode smashon("9003234");
  // there is a double for mass and a float for width. This is
  // intentional.
  ParticleType A("smashon", 3.243, 0.234f, smashon);
  COMPARE(A.name(), "smashon");
  COMPARE(A.mass(), 3.243f);
  COMPARE(A.width(), 0.234f);
  COMPARE(A.pdgcode(), smashon);
  COMPARE(A.isospin(), smashon.isospin_total());
  COMPARE(A.charge(), smashon.charge());
  COMPARE(A.spin(), smashon.spin());
}
