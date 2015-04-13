/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "setup.h"
#include "../include/pdgcode.h"

using namespace Smash;

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
      "smashon 0.123 1.2 661\n"
      "pi0 0.1350 -1.0 111\n"
      "pi+ 0.1396 -1.0 211\n"
      "rho0 0.7755 0.149 113\n"
      "rho+ 0.7755 0.149 213\n"
      "eta 0.5479 1.0e-6 221\n"
      "omega 0.7827 0.0085 223\n"
      "p 0.9383 -1.0 2212\n"
      "n 0.9396 -1.0 2112\n"
      "Delta++ 1.232 0.117 2224\n"
      "Delta+ 1.232 0.117 2214\n"
      "Delta0 1.232 0.117 2114\n"
      "Delta- 1.232 0.117 1114\n");
}

TEST(create_particledata_piplus) {
  PdgCode pdg = 0x211;
  ParticleData p{ParticleType::find(pdg)};

  COMPARE(p.id(), -1);
  COMPARE(p.pdgcode(), pdg);
  COMPARE(p.id_process(), -1);
  COMPARE(p.momentum().x0(), 0.0);
  COMPARE(p.momentum().x1(), 0.0);
  COMPARE(p.momentum().x2(), 0.0);
  COMPARE(p.momentum().x3(), 0.0);
  COMPARE(p.position().x0(), 0.0);
  COMPARE(p.position().x1(), 0.0);
  COMPARE(p.position().x2(), 0.0);
  COMPARE(p.position().x3(), 0.0);

  p.set_id(2);
  COMPARE(p.id(), 2);
}

TEST(set_get) {
  PdgCode smashon = 0x661;
  ParticleData p = Test::smashon();
  p.set_id(4);
  COMPARE(p.id(), 4);
  COMPARE(p.pdgcode(), smashon);
  COMPARE(p.is_hadron(), smashon.is_hadron());
  p.set_id_process(5);
  COMPARE(p.id_process(), 5);
  p.set_id_process(6);
  COMPARE(p.id_process(), 6);
  FourVector m(1.0, 1.2, 1.4, 1.6);
  p.set_4momentum(m);
  COMPARE(p.momentum(), FourVector(1.0, 1.2, 1.4, 1.6));
  ThreeVector M(1.1, 1.3, 1.5);
  p.set_4momentum(1.0, M);
  COMPARE(p.momentum(), FourVector(sqrt(1.0 + M.sqr()), 1.1, 1.3, 1.5));
}

TEST(set_get2) {
  ParticleData p = Test::smashon(Test::Position{3.5, 3.6, 3.7, 2345.3});
  COMPARE(p.position().x0(), 3.5);
  COMPARE(p.position().x1(), 3.6);
  COMPARE(p.position().x2(), 3.7);
  COMPARE(p.position().x3(), 2345.3);
  ThreeVector M(2.1, 2.3, 2.5);
  p.set_4momentum(2.0, M.x1(), M.x2(), M.x3());
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
  ParticleData p = Test::smashon(1);
  ParticleData q = Test::smashon(2);
  ParticleData r = Test::smashon(1);
  VERIFY(!(p == q));
  VERIFY(  p == r );
  VERIFY(  p == 1 );
  VERIFY(  p <  2 );
  VERIFY(  p <  q );
}
