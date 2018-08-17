/*
 *
 *    Copyright (c) 2014-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "setup.h"

#include "../include/smash/decaymodes.h"
#include "../include/smash/isoparticletype.h"
#include "../include/smash/particletype.h"

using namespace smash;

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
      "π               0.138   7.7e-9      111     211\n"
      "σ               0.800   0.400   9000221\n"
      "ρ               0.776   0.149       113     213\n"
      "ω               0.783   8.49e-3     223\n"
      "N       0.938 0         2112    2212\n"
      "Δ       1.232 0.117    1114    2114    2214    2224\n"
      "Λ        1.116 0         3122\n"
      "Λ(1520)  1.520 0.0156    3124\n"
      "Λ(1690)  1.690 0.0600   13124\n"
      "Σ       1.189 0        3112    3212    3222\n"
      "e⁻ 0.000511 0 11\n");
}

TEST_CATCH(load_decaymodes_missing_pdg,
           IsoParticleType::ParticleNotFoundFailure) {
  const std::string decays_input("unknown_particle \n");
  DecayModes::load_decaymodes(decays_input);
}

TEST_CATCH(load_decaymodes_no_decays, DecayModes::MissingDecays) {
  const std::string decays_input("ρ  # rho\n");
  DecayModes::load_decaymodes(decays_input);
}

TEST_CATCH(load_decaymodes_incorrect_start,
           IsoParticleType::ParticleNotFoundFailure) {
  const std::string decays_input("ρ⁺  # rho+\n");
  DecayModes::load_decaymodes(decays_input);
}

TEST_CATCH(load_decaymodes_duplicate, DecayModes::LoadFailure) {
  const std::string decays_input(
      "σ \n"
      "1.  0  π π\n"
      "\n"
      "σ \n"
      "1.  0  π π\n");
  DecayModes::load_decaymodes(decays_input);
}

const double tolerance = 2.0e-7;

const std::string decays_input(
    " ρ\t# rho\n"
    "\n"
    "0.99\t1\tπ π\t# pi pi \n"
    "0.01  1  e⁻ e⁺\n"
    "\n"
    "\n"
    "ω      # omega\n"
    "1. 0 π ρ   # pi rho\n"
    "\n"
    "Δ\n"
    "1.  1  N π\n"
    "\n"
    "Λ(1520)\n"
    "1.0 1 Λ π π\n"
    "\n"
    "Λ(1690)\n"
    "1.0 1 Σ π π\n"
    "\n"
    "σ \n"
    "1.  0  π π\n");

TEST(load_decay_modes) {
  DecayModes::load_decaymodes(decays_input);

  UnitTest::setFuzzyness<double>(2);
  // check that the decays of the rho and omega are generated correctly
  {
    const auto &rho_0 = ParticleType::find(0x113).decay_modes();
    VERIFY(!rho_0.is_empty());
    const auto &modelist = rho_0.decay_mode_list();
    COMPARE(modelist.size(), 2u);
    COMPARE_ABSOLUTE_ERROR(modelist[0]->weight(), 0.99, tolerance);
    COMPARE(modelist[0]->particle_number(), 2u);
    COMPARE(modelist[0]->particle_types()[0]->pdgcode(), 0x211);
    COMPARE(modelist[0]->particle_types()[1]->pdgcode(), -0x211);
    COMPARE(modelist[1]->weight(), 0.01);
    COMPARE(modelist[1]->particle_number(), 2u);
    COMPARE(modelist[1]->particle_types()[0]->pdgcode(), 0x11);
    COMPARE(modelist[1]->particle_types()[1]->pdgcode(), -0x11);
  }
  {
    const auto &rhoplus = ParticleType::find(0x213).decay_modes();
    VERIFY(!rhoplus.is_empty());
    const auto &modelist = rhoplus.decay_mode_list();
    COMPARE(modelist.size(), 1u);
    COMPARE(modelist[0]->weight(), 1.);
    COMPARE(modelist[0]->particle_number(), 2u);
    COMPARE(modelist[0]->particle_types()[0]->pdgcode(), 0x111);
    COMPARE(modelist[0]->particle_types()[1]->pdgcode(), 0x211);
  }
  {
    const auto &rhominus = ParticleType::find(-0x213).decay_modes();
    VERIFY(!rhominus.is_empty());
    const auto &modelist = rhominus.decay_mode_list();
    COMPARE(modelist.size(), 1u);
    COMPARE(modelist[0]->weight(), 1.);
    COMPARE(modelist[0]->particle_number(), 2u);
    COMPARE(modelist[0]->particle_types()[0]->pdgcode(), 0x111);
    COMPARE(modelist[0]->particle_types()[1]->pdgcode(), -0x211);
  }
  {
    const auto &omega = ParticleType::find(0x223).decay_modes();
    VERIFY(!omega.is_empty());
    const auto &modelist = omega.decay_mode_list();
    COMPARE(modelist.size(), 3u);
    COMPARE_ABSOLUTE_ERROR(modelist[0]->weight(), 1. / 3., tolerance);
    COMPARE(modelist[0]->particle_number(), 2u);
    COMPARE(modelist[0]->particle_types()[0]->pdgcode(), 0x111);
    COMPARE(modelist[0]->particle_types()[1]->pdgcode(), 0x113);
    COMPARE_ABSOLUTE_ERROR(modelist[1]->weight(), 1. / 3., tolerance);
    COMPARE(modelist[1]->particle_number(), 2u);
    COMPARE(modelist[1]->particle_types()[0]->pdgcode(), 0x211);
    COMPARE(modelist[1]->particle_types()[1]->pdgcode(), -0x213);
    COMPARE_ABSOLUTE_ERROR(modelist[2]->weight(), 1. / 3., tolerance);
    COMPARE(modelist[2]->particle_number(), 2u);
    COMPARE(modelist[2]->particle_types()[0]->pdgcode(), -0x211);
    COMPARE(modelist[2]->particle_types()[1]->pdgcode(), 0x213);
  }
  // check that the decays of the anti-Delta multiplet are generated correctly
  {
    // anti-Delta--
    const auto &Delta = ParticleType::find(-0x2224).decay_modes();
    VERIFY(!Delta.is_empty());
    const auto &modelist = Delta.decay_mode_list();
    COMPARE(modelist.size(), 1u);
    COMPARE(modelist[0]->weight(), 1.);
    COMPARE(modelist[0]->particle_number(), 2u);
    COMPARE(modelist[0]->particle_types()[0]->pdgcode(), -0x2212);
    COMPARE(modelist[0]->particle_types()[1]->pdgcode(), -0x211);
  }
  {
    // anti-Delta-
    const auto &Delta = ParticleType::find(-0x2214).decay_modes();
    VERIFY(!Delta.is_empty());
    const auto &modelist = Delta.decay_mode_list();
    COMPARE(modelist.size(), 2u);
    COMPARE_ABSOLUTE_ERROR(modelist[0]->weight(), 1. / 3., tolerance);
    COMPARE(modelist[0]->particle_number(), 2u);
    COMPARE(modelist[0]->particle_types()[0]->pdgcode(), -0x2112);
    COMPARE(modelist[0]->particle_types()[1]->pdgcode(), -0x211);
    FUZZY_COMPARE(modelist[1]->weight(), 2. / 3.);
    COMPARE(modelist[1]->particle_number(), 2u);
    COMPARE(modelist[1]->particle_types()[0]->pdgcode(), -0x2212);
    COMPARE(modelist[1]->particle_types()[1]->pdgcode(), 0x111);
  }
  {
    // anti-Delta0
    const auto &Delta = ParticleType::find(-0x2114).decay_modes();
    VERIFY(!Delta.is_empty());
    const auto &modelist = Delta.decay_mode_list();
    COMPARE(modelist.size(), 2u);
    FUZZY_COMPARE(modelist[0]->weight(), 2. / 3.);
    COMPARE(modelist[0]->particle_number(), 2u);
    COMPARE(modelist[0]->particle_types()[0]->pdgcode(), -0x2112);
    COMPARE(modelist[0]->particle_types()[1]->pdgcode(), 0x111);
    COMPARE_ABSOLUTE_ERROR(modelist[1]->weight(), 1. / 3., tolerance);
    COMPARE(modelist[1]->particle_number(), 2u);
    COMPARE(modelist[1]->particle_types()[0]->pdgcode(), -0x2212);
    COMPARE(modelist[1]->particle_types()[1]->pdgcode(), 0x211);
  }
  {
    // anti-Delta+
    const auto &Delta = ParticleType::find(-0x1114).decay_modes();
    VERIFY(!Delta.is_empty());
    const auto &modelist = Delta.decay_mode_list();
    COMPARE(modelist.size(), 1u);
    COMPARE(modelist[0]->weight(), 1.);
    COMPARE(modelist[0]->particle_number(), 2u);
    COMPARE(modelist[0]->particle_types()[0]->pdgcode(), -0x2112);
    COMPARE(modelist[0]->particle_types()[1]->pdgcode(), 0x211);
  }
}

TEST(load_decaymodes_3body) {
  DecayModes::load_decaymodes(decays_input);
  {
    const auto &antiLambda = ParticleType::find(-0x3124).decay_modes();
    VERIFY(!antiLambda.is_empty());
    const auto &modelist = antiLambda.decay_mode_list();
    COMPARE(modelist.size(), 2u);
    COMPARE_ABSOLUTE_ERROR(modelist[0]->weight(), 1. / 3., tolerance);
    COMPARE_ABSOLUTE_ERROR(modelist[1]->weight(), 2. / 3., tolerance);
    for (int i = 0; i < 2; i++) {
      COMPARE(modelist[i]->particle_types()[2]->pdgcode(), -0x3122);
    }
  }
  {
    const auto &antiLambda = ParticleType::find(-0x13124).decay_modes();
    VERIFY(!antiLambda.is_empty());
    const auto &modelist = antiLambda.decay_mode_list();
    COMPARE(modelist.size(), 3u);
    for (int i = 0; i < 3; i++) {
      COMPARE_ABSOLUTE_ERROR(modelist[i]->weight(), 1. / 3., tolerance);
      int charge = 0;
      // 3 neutral particles are forbidden by isospin.
      bool all_charges_are_zero = true;
      for (const auto &p : modelist[i]->particle_types()) {
        charge += p->charge();
        all_charges_are_zero &= (p->charge() == 0);
      }
      VERIFY(!all_charges_are_zero);
      COMPARE(charge, 0);
    }
  }
  {
    const auto &omega = ParticleType::find(0x223).decay_modes();
    VERIFY(!omega.is_empty());
    const auto &modelist = omega.decay_mode_list();
    COMPARE(modelist.size(), 3u);
    COMPARE_ABSOLUTE_ERROR(modelist[0]->weight(), 1. / 3., tolerance);
    VERIFY(modelist[0]->particle_types()[0]->pdgcode().is_pion());
  }
}

TEST_CATCH(add_no_particles, DecayModes::InvalidDecay) {
  DecayModes m;
  m.add_mode(&ParticleType::find(0x113), 1., 0, {});
}

TEST_CATCH(add_one_particle, DecayModes::InvalidDecay) {
  DecayModes m;
  m.add_mode(&ParticleType::find(0x113), 1., 0, {&ParticleType::find(0x211)});
}

TEST(add_two_particles) {
  DecayModes m;
  VERIFY(m.is_empty());
  m.add_mode(&ParticleType::find(0x113), 1., 0,
             {&ParticleType::find(0x211), &ParticleType::find(-0x211)});
  VERIFY(!m.is_empty());
}
