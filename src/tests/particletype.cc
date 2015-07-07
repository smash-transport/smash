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

TEST(assign) {
  PdgCode smashon("9003234");
  // there is a double for mass and a float for width. This is
  // intentional.
  ParticleType A("smashon", 3.243, 0.234f, smashon);
#ifdef NDEBUG
  COMPARE(A.name(), std::string{});
#else
  COMPARE(A.name(), "smashon");
#endif
  COMPARE(A.mass(), 3.243f);
  COMPARE(A.width_at_pole(), 0.234f);
  COMPARE(A.pdgcode(), smashon);
  COMPARE(A.isospin(), smashon.isospin_total());
  COMPARE(A.charge(), smashon.charge());
  COMPARE(A.spin(), smashon.spin());
}

TEST_CATCH(load_from_incorrect_string, ParticleType::LoadFailure) {
  ParticleType::create_type_list("Hallo Welt! (wave)");
}

TEST_CATCH(load_one_particle_with_incorrect_newline, ParticleType::LoadFailure) {
  const std::string parts("pi0 0.1350\n-1.0 111");
  ParticleType::create_type_list(parts);
}

TEST_CATCH(load_duplicate_particle, ParticleType::LoadFailure) {
  ParticleType::create_type_list("π⁺ 0.138 0.0 211\nπ⁺ 0.138 0.0 211\n");
}

TEST(create_type_list) {
  ParticleType::create_type_list(
      "test0 0.1350 -1.0 661\n"
      "# Hello\n"
      "  # Will you ignore me? #### sldfkjsdf\n"
      "\t\t  \t # yes?\n"
      "test1 0.1350  -2.0 663 # This is pi0. Swell.\n"
      "\t\n\t  test2  0.1350 \t -3.0 665\n "
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

  COMPARE(ParticleType::list_all().size(), 23u);

  ParticleTypePtr type = &ParticleType::find(0x661);
  COMPARE(type->mass(), 0.135f);
  COMPARE(type->width_at_pole(), -1.f);
  COMPARE(type->pdgcode(), PdgCode(0x661));
  type = &ParticleType::find(0x663);
  COMPARE(type->mass(), 0.135f);
  COMPARE(type->width_at_pole(), -2.f);
  COMPARE(type->pdgcode(), PdgCode(0x663));
  type = &ParticleType::find(0x665);
  COMPARE(type->mass(), 0.135f);
  COMPARE(type->width_at_pole(), -3.f);
  COMPARE(type->pdgcode(), PdgCode(0x665));

  type = &ParticleType::find(-0x1114);
  COMPARE(type->mass(), 1.232f);
  COMPARE(type->width_at_pole(), .117f);
  COMPARE(type->pdgcode().dump(), 0x80001114);
  COMPARE(type->isospin(), 3);
  COMPARE(type->charge(), 1);
  COMPARE(type->spin(), 3u);

  type = &ParticleType::find(0x2112);
  COMPARE(type->mass(), .9396f);
  COMPARE(type->width_at_pole(), -1.f);
  COMPARE(type->pdgcode().dump(), 0x2112u);
  COMPARE(type->isospin(), 1);
  COMPARE(type->charge(), 0);
  COMPARE(type->spin(), 1u);
}

TEST(list_all_iteration) {
  std::size_t count = 0;
  for (const auto &type : ParticleType::list_all()) {
    const PdgCode pdg = type.pdgcode();
    const ParticleType &type2 = ParticleType::find(pdg);
    COMPARE(&type, &type2);
    ++count;
  }
  COMPARE(count, ParticleType::list_all().size());
}

TEST(exists) {
  VERIFY(ParticleType::exists(0x111));
  VERIFY(!ParticleType::exists(0x667));
}

