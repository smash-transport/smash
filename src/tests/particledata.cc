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

TEST(init_default) {
  ParticleData p;
  COMPARE(p.id(), -1);
  COMPARE(p.pdgcode(), PdgCode::invalid());
  COMPARE(p.id_process(), -1);
  COMPARE(p.collision_time(), 0.0);
  COMPARE(p.momentum().x0(), 0.0);
  COMPARE(p.momentum().x1(), 0.0);
  COMPARE(p.momentum().x2(), 0.0);
  COMPARE(p.momentum().x3(), 0.0);
  COMPARE(p.position().x0(), 0.0);
  COMPARE(p.position().x1(), 0.0);
  COMPARE(p.position().x2(), 0.0);
  COMPARE(p.position().x3(), 0.0);
}

TEST(init_from_pdgcode) {
  PdgCode pdg = 0x211;
  ParticleData p{pdg};

  COMPARE(p.id(), -1);
  COMPARE(p.pdgcode(), pdg);
  COMPARE(p.id_process(), -1);
  COMPARE(p.collision_time(), 0.0);
  COMPARE(p.momentum().x0(), 0.0);
  COMPARE(p.momentum().x1(), 0.0);
  COMPARE(p.momentum().x2(), 0.0);
  COMPARE(p.momentum().x3(), 0.0);
  COMPARE(p.position().x0(), 0.0);
  COMPARE(p.position().x1(), 0.0);
  COMPARE(p.position().x2(), 0.0);
  COMPARE(p.position().x3(), 0.0);
}

TEST(init_id) {
  ParticleData p(2);
  COMPARE(p.id(), 2);
  COMPARE(p.pdgcode(), PdgCode::invalid());
  COMPARE(p.id_process(), -1);
  COMPARE(p.collision_time(), 0.0);
  COMPARE(p.momentum().x0(), 0.0);
  COMPARE(p.momentum().x1(), 0.0);
  COMPARE(p.momentum().x2(), 0.0);
  COMPARE(p.momentum().x3(), 0.0);
  COMPARE(p.position().x0(), 0.0);
  COMPARE(p.position().x1(), 0.0);
  COMPARE(p.position().x2(), 0.0);
  COMPARE(p.position().x3(), 0.0);
}

TEST(set_get) {
  ParticleData p;
  p.set_id(4);
  COMPARE(p.id(), 4);
  PdgCode smashon = 0x12345;
  p.set_pdgcode(smashon);
  COMPARE(p.pdgcode(), smashon);
  COMPARE(p.is_hadron(), smashon.is_hadron());
  p.set_id_process(5);
  COMPARE(p.id_process(), 5);
  p.set_collision_time(1.234);
  COMPARE(p.collision_time(), 1.234);
  p.set_collision(2.345);
  COMPARE(p.collision_time(), 2.345);
  p.set_collision_past(6);
  COMPARE(p.collision_time(), 0.0);
  COMPARE(p.id_process(), 6);
  FourVector m(1.0, 1.2, 1.4, 1.6);
  p.set_momentum(m);
  COMPARE(p.momentum().x0(), 1.0);
  COMPARE(p.momentum().x1(), 1.2);
  COMPARE(p.momentum().x2(), 1.4);
  COMPARE(p.momentum().x3(), 1.6);
  ThreeVector M(1.1, 1.3, 1.5);
  p.set_momentum(1.0, M);
  COMPARE(p.momentum().x0(), sqrt(1.0 + M.sqr()));
  COMPARE(p.momentum().x1(), 1.1);
  COMPARE(p.momentum().x2(), 1.3);
  COMPARE(p.momentum().x3(), 1.5);
}

TEST(set_get2) {
  ParticleData p;
  FourVector X(3.5, 3.6, 3.7, 2345.3);
  p.set_position(X);
  COMPARE(p.position().x0(), 3.5);
  COMPARE(p.position().x1(), 3.6);
  COMPARE(p.position().x2(), 3.7);
  COMPARE(p.position().x3(), 2345.3);
  ThreeVector M(2.1, 2.3, 2.5);
  p.set_momentum(2.0, M.x1(), M.x2(), M.x3());
  COMPARE(p.momentum().x0(), sqrt(4.0 + M.sqr()));
  COMPARE(p.momentum().x1(), 2.1);
  COMPARE(p.momentum().x2(), 2.3);
  COMPARE(p.momentum().x3(), 2.5);
  ThreeVector v = p.velocity();
  COMPARE(v.x1(), 2.1/sqrt(4.0 + M.sqr()));
  COMPARE(v.x2(), 2.3/sqrt(4.0 + M.sqr()));
  COMPARE(v.x3(), 2.5/sqrt(4.0 + M.sqr()));
  p.boost(v);
  COMPARE_RELATIVE_ERROR(p.momentum().x0(), 2.0, 1e-15);
  COMPARE_ABSOLUTE_ERROR(p.momentum().x1(), 0.0, 1e-15);
  COMPARE_ABSOLUTE_ERROR(p.momentum().x2(), 0.0, 1e-15);
  COMPARE_ABSOLUTE_ERROR(p.momentum().x3(), 0.0, 1e-15);
}

TEST(comparisons) {
  ParticleData p(1);
  ParticleData q(2);
  ParticleData r(1);
  VERIFY(!(p == q));
  VERIFY(  p == r );
  VERIFY(  p == 1 );
  VERIFY(  p <  2 );
  VERIFY(  p <  q );
}
