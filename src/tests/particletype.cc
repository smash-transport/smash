/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "../include/particles.h"
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

TEST_CATCH(load_from_incorrect_string, Particles::LoadFailure) {
  ParticleType::create_type_list("Hallo Welt! (wave)");
}

TEST_CATCH(load_one_particle_with_incorrect_newline, Particles::LoadFailure) {
  const std::string parts("pi0 0.1350\n-1.0 111");
  ParticleType::create_type_list(parts);
}

TEST(create_type_list) {
  ParticleType::create_type_list(
      "test0 0.1350 -1.0 -331\n"
      "# Hello\n"
      "  # Will you ignore me? #### sldfkjsdf\n"
      "\t\t  \t # yes?\n"
      "test1 0.1350  -2.0 -441 # This is pi0. Swell.\n"
      "\t\n\t  test2  0.1350 \t -3.0 -551\n "
      "pi0 0.1350 -1.0 111\n"
      "pi+ 0.1396 -1.0 211\n"
      "pi- 0.1396 -1.0 -211\n"
      "rho0 0.7755 0.149 113\n"
      "rho+ 0.7755 0.149 213\n"
      "rho- 0.7755 0.149 -213\n"
      "eta 0.5479 1.0e-6 221\n"
      "omega 0.7827 0.0085 223\n"
      "p 0.9383 -1.0 2212\n"
      "pbar 0.9383 -1.0 -2212\n"
      "n 0.9396 -1.0 2112\n"
      "nbar 0.9396 -1.0 -2112\n"
      "Delta++ 1.232 0.117 2224\n"
      "Delta+ 1.232 0.117 2214\n"
      "Delta0 1.232 0.117 2114\n"
      "Delta- 1.232 0.117 1114\n"
      "Deltabar++ 1.232 0.117 -2224\n"
      "Deltabar+ 1.232 0.117 -2214\n"
      "Deltabar0 1.232 0.117 -2114\n"
      "Deltabar- 1.232 0.117 -1114\n");

  COMPARE(ParticleType::list_all().size(), 23u);

  ParticleType type = ParticleType::find(-0x331);
  COMPARE(type.mass(), 0.135f);
  COMPARE(type.width(), -1.f);
  COMPARE(type.pdgcode(), PdgCode(-0x331));
  type = ParticleType::find(-0x441);
  COMPARE(type.mass(), 0.135f);
  COMPARE(type.width(), -2.f);
  COMPARE(type.pdgcode(), PdgCode(-0x441));
  type = ParticleType::find(-0x551);
  COMPARE(type.mass(), 0.135f);
  COMPARE(type.width(), -3.f);
  COMPARE(type.pdgcode(), PdgCode(-0x551));

  type = ParticleType::find(-0x1114);
  COMPARE(type.mass(), 1.232f);
  COMPARE(type.width(), .117f);
  COMPARE(type.pdgcode().dump(), 0x80001114);
  COMPARE(type.isospin(), 3);
  COMPARE(type.charge(), 1);
  COMPARE(type.spin(), 3);

  type = ParticleType::find(0x2112);
  COMPARE(type.mass(), .9396f);
  COMPARE(type.width(), -1.f);
  COMPARE(type.pdgcode().dump(), 0x2112u);
  COMPARE(type.isospin(), 1);
  COMPARE(type.charge(), 0);
  COMPARE(type.spin(), 1);
}

TEST(exists) {
  VERIFY(ParticleType::exists(0x111));
  VERIFY(!ParticleType::exists(-0x661));
}

