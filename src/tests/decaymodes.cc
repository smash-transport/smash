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

#include "../include/decaymodes.h"
#include "../include/isoparticletype.h"
#include "../include/particletype.h"

using namespace Smash;

TEST(init_particle_types) {
  Test::create_actual_particletypes();
}

TEST_CATCH(load_decaymodes_missing_pdg, IsoParticleType::ParticleNotFoundFailure) {
  const std::string decays_input(
      "unknown_particle \n"
      );
  DecayModes::load_decaymodes(decays_input);
}

TEST_CATCH(load_decaymodes_no_decays, DecayModes::MissingDecays) {
  const std::string decays_input(
      "ρ  # rho\n"
      );
  DecayModes::load_decaymodes(decays_input);
}

TEST_CATCH(load_decaymodes_incorrect_start, IsoParticleType::ParticleNotFoundFailure) {
  const std::string decays_input(
      "ρ⁺  # rho+\n"
      );
  DecayModes::load_decaymodes(decays_input);
}

TEST_CATCH(load_decaymodes_duplicate, DecayModes::LoadFailure) {
  const std::string decays_input(
      "σ \n"
      "1.  0  π π\n"
      "\n"
      "σ \n"
      "1.  0  π π\n"
      );
  DecayModes::load_decaymodes(decays_input);
}


const float tolerance = 1.0e-7;

TEST(load_decay_modes) {
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
      );
  DecayModes::load_decaymodes(decays_input);

  // check that the decays of the rho and omega are generated correctly
  {
    const auto &rho_0 = ParticleType::find(0x113).decay_modes();
    VERIFY(!rho_0.is_empty());
    const auto &modelist = rho_0.decay_mode_list();
    COMPARE(modelist.size(), 3u);
    COMPARE_ABSOLUTE_ERROR(modelist[0]->weight(), 0.495f, tolerance);
    COMPARE(modelist[0]->particle_number(), 2u);
    COMPARE(modelist[0]->particle_types()[0]->pdgcode(),  0x211);
    COMPARE(modelist[0]->particle_types()[1]->pdgcode(), -0x211);
    COMPARE_ABSOLUTE_ERROR(modelist[1]->weight(), 0.495f, tolerance);
    COMPARE(modelist[1]->particle_number(), 2u);
    COMPARE(modelist[1]->particle_types()[0]->pdgcode(), -0x211);
    COMPARE(modelist[1]->particle_types()[1]->pdgcode(),  0x211);
    COMPARE(modelist[2]->weight(), 0.01f);
    COMPARE(modelist[2]->particle_number(), 2u);
    COMPARE(modelist[2]->particle_types()[0]->pdgcode(),  0x11);
    COMPARE(modelist[2]->particle_types()[1]->pdgcode(), -0x11);
  }
  {
    const auto &rhoplus = ParticleType::find(0x213).decay_modes();
    VERIFY(!rhoplus.is_empty());
    const auto &modelist = rhoplus.decay_mode_list();
    COMPARE(modelist.size(), 2u);
    COMPARE(modelist[0]->weight(), 0.5);
    COMPARE(modelist[0]->particle_number(), 2u);
    COMPARE(modelist[0]->particle_types()[0]->pdgcode(), 0x111);
    COMPARE(modelist[0]->particle_types()[1]->pdgcode(), 0x211);
    COMPARE(modelist[1]->weight(), 0.5);
    COMPARE(modelist[1]->particle_number(), 2u);
    COMPARE(modelist[1]->particle_types()[0]->pdgcode(), 0x211);
    COMPARE(modelist[1]->particle_types()[1]->pdgcode(), 0x111);
  }
  {
    const auto &rhominus = ParticleType::find(-0x213).decay_modes();
    VERIFY(!rhominus.is_empty());
    const auto &modelist = rhominus.decay_mode_list();
    COMPARE(modelist.size(), 2u);
    COMPARE(modelist[0]->weight(), 0.5);
    COMPARE(modelist[0]->particle_number(), 2u);
    COMPARE(modelist[0]->particle_types()[0]->pdgcode(),  0x111);
    COMPARE(modelist[0]->particle_types()[1]->pdgcode(), -0x211);
    COMPARE(modelist[1]->weight(), 0.5);
    COMPARE(modelist[1]->particle_number(), 2u);
    COMPARE(modelist[1]->particle_types()[0]->pdgcode(), -0x211);
    COMPARE(modelist[1]->particle_types()[1]->pdgcode(),  0x111);
  }
  {
    const auto &omega = ParticleType::find(0x223).decay_modes();
    VERIFY(!omega.is_empty());
    const auto &modelist = omega.decay_mode_list();
    COMPARE(modelist.size(), 3u);
    COMPARE_ABSOLUTE_ERROR(modelist[0]->weight(), 1.f/3.f, tolerance);
    COMPARE(modelist[0]->particle_number(), 2u);
    COMPARE(modelist[0]->particle_types()[0]->pdgcode(), 0x111);
    COMPARE(modelist[0]->particle_types()[1]->pdgcode(), 0x113);
    COMPARE_ABSOLUTE_ERROR(modelist[1]->weight(), 1.f/3.f, tolerance);
    COMPARE(modelist[1]->particle_number(), 2u);
    COMPARE(modelist[1]->particle_types()[0]->pdgcode(),  0x211);
    COMPARE(modelist[1]->particle_types()[1]->pdgcode(), -0x213);
    COMPARE_ABSOLUTE_ERROR(modelist[2]->weight(), 1.f/3.f, tolerance);
    COMPARE(modelist[2]->particle_number(), 2u);
    COMPARE(modelist[2]->particle_types()[0]->pdgcode(), -0x211);
    COMPARE(modelist[2]->particle_types()[1]->pdgcode(),  0x213);
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
    COMPARE(modelist[0]->particle_types()[1]->pdgcode(),  -0x211);
  }
  {
    // anti-Delta-
    const auto &Delta = ParticleType::find(-0x2214).decay_modes();
    VERIFY(!Delta.is_empty());
    const auto &modelist = Delta.decay_mode_list();
    COMPARE(modelist.size(), 2u);
    COMPARE_ABSOLUTE_ERROR(modelist[0]->weight(), 1.f/3.f, tolerance);
    COMPARE(modelist[0]->particle_number(), 2u);
    COMPARE(modelist[0]->particle_types()[0]->pdgcode(), -0x2112);
    COMPARE(modelist[0]->particle_types()[1]->pdgcode(),  -0x211);
    COMPARE(modelist[1]->weight(), 2.f/3.f);
    COMPARE(modelist[1]->particle_number(), 2u);
    COMPARE(modelist[1]->particle_types()[0]->pdgcode(), -0x2212);
    COMPARE(modelist[1]->particle_types()[1]->pdgcode(),   0x111);
  }
  {
    // anti-Delta0
    const auto &Delta = ParticleType::find(-0x2114).decay_modes();
    VERIFY(!Delta.is_empty());
    const auto &modelist = Delta.decay_mode_list();
    COMPARE(modelist.size(), 2u);
    COMPARE(modelist[0]->weight(), 2.f/3.f);
    COMPARE(modelist[0]->particle_number(), 2u);
    COMPARE(modelist[0]->particle_types()[0]->pdgcode(), -0x2112);
    COMPARE(modelist[0]->particle_types()[1]->pdgcode(),   0x111);
    COMPARE_ABSOLUTE_ERROR(modelist[1]->weight(), 1.f/3.f, tolerance);
    COMPARE(modelist[1]->particle_number(), 2u);
    COMPARE(modelist[1]->particle_types()[0]->pdgcode(), -0x2212);
    COMPARE(modelist[1]->particle_types()[1]->pdgcode(),   0x211);
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
    COMPARE(modelist[0]->particle_types()[1]->pdgcode(),   0x211);
  }
}

TEST(load_decaymodes_3body) {
  const std::string decays_input(
      "Λ(1520)\n"
      "1.0 1 Λ π π\n"
      "\n"
      "Λ(1690)\n"
      "1.0 1 Σ π π\n"
      );
  DecayModes::load_decaymodes(decays_input);
  {
    const auto &Lambda = ParticleType::find(-0x3124).decay_modes();
    VERIFY(!Lambda.is_empty());
    const auto &modelist = Lambda.decay_mode_list();
    COMPARE(modelist.size(), 3u);
    for (int i = 0; i < 3; i++) {
      COMPARE_ABSOLUTE_ERROR(modelist[i]->weight(), 1.f/3.f, tolerance);
      COMPARE(modelist[i]->particle_types()[0]->pdgcode(), -0x3122);
    }
  }
  {
    const auto &Lambda = ParticleType::find(-0x13124).decay_modes();
    VERIFY(!Lambda.is_empty());
    const auto &modelist = Lambda.decay_mode_list();
    COMPARE(modelist.size(), 6u);
    for (int i = 0; i < 6; i++) {
      COMPARE_ABSOLUTE_ERROR(modelist[i]->weight(), 1.f/6.f, tolerance);
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
}


TEST_CATCH(add_no_particles, DecayModes::InvalidDecay) {
  DecayModes m;
  m.add_mode(1.f, 0, {});
}

TEST_CATCH(add_one_particle, DecayModes::InvalidDecay) {
  DecayModes m;
  m.add_mode(1.f, 0, {&ParticleType::find(0x211)});
}

TEST(add_two_particles) {
  DecayModes m;
  VERIFY(m.is_empty());
  m.add_mode(1.f, 0,
             {&ParticleType::find(0x211), &ParticleType::find(-0x211)});
  VERIFY(!m.is_empty());
}
