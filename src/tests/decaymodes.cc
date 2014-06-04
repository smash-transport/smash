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
#include "../include/decaymodes.h"

using namespace Smash;

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
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
}

TEST_CATCH(load_decaymodes_missing_pdg, DecayModes::ReferencedParticleNotFound) {
  const std::string decays_input(
      "-441 \n"
      );
  DecayModes::load_decaymodes(decays_input);
}

TEST_CATCH(load_decaymodes_no_decays, DecayModes::MissingDecays) {
  const std::string decays_input(
      "113 # rho0\n"
      );
  DecayModes::load_decaymodes(decays_input);
}

TEST_CATCH(load_decaymodes_incorrect_start, PdgCode::InvalidPdgCode) {
  const std::string decays_input(
      "113. # rho0\n"
      );
  DecayModes::load_decaymodes(decays_input);
}

TEST(load_decaymodes_two_channels) {
  const std::string decays_input(
      " 113\t# rho0\n"
      "\n"
      " 1.0\t1\t211 -211\t# pi+ pi- \n"
      " \n"
      "\n"
      "223      # omega\n"
      "0.33 0 111 113   # pi0 rho0\n"
      "\n"
      "0.33 0 211 -213  # pi+ rho-\n"
      "0.33 0 -211 213  # pi- rho+\n"
      );
  DecayModes::load_decaymodes(decays_input);

  {
    const auto &rho0 = DecayModes::find(0x113);
    VERIFY(!rho0.is_empty());
    const auto &modelist = rho0.decay_mode_list();
    COMPARE(modelist.size(), 1u);
    COMPARE(modelist[0].weight(), 1.);
    COMPARE(modelist[0].pdg_list().size(), 2u);
    COMPARE(modelist[0].pdg_list()[0].dump(), 0x211u);
    COMPARE(modelist[0].pdg_list()[1].dump(), 0x80000211u);
  }
  {
    const auto &omega = DecayModes::find(0x223);
    VERIFY(!omega.is_empty());
    const auto &modelist = omega.decay_mode_list();
    COMPARE(modelist.size(), 3u);
    FUZZY_COMPARE(float(modelist[0].weight()), 1.f/3.f);
    FUZZY_COMPARE(float(modelist[1].weight()), 1.f/3.f);
    FUZZY_COMPARE(float(modelist[2].weight()), 1.f/3.f);
    COMPARE(modelist[0].pdg_list().size(), 2u);
    COMPARE(modelist[0].pdg_list()[0].dump(), 0x111u);
    COMPARE(modelist[0].pdg_list()[1].dump(), 0x113u);
    COMPARE(modelist[1].pdg_list().size(), 2u);
    COMPARE(modelist[1].pdg_list()[0].dump(), 0x211u);
    COMPARE(modelist[1].pdg_list()[1].dump(), 0x80000213u);
    COMPARE(modelist[2].pdg_list().size(), 2u);
    COMPARE(modelist[2].pdg_list()[0].dump(), 0x80000211u);
    COMPARE(modelist[2].pdg_list()[1].dump(), 0x213u);
  }
}

TEST_CATCH(add_no_particles, DecayModes::InvalidDecay) {
  DecayModes m;
  m.add_mode(1.f, 0, {});
}

TEST_CATCH(add_one_particle, DecayModes::InvalidDecay) {
  DecayModes m;
  m.add_mode(1.f, 0, {0});
}

TEST(add_two_particles) {
  DecayModes m;
  VERIFY(m.is_empty());
  m.add_mode(1.f, 0, {0, 1});
  VERIFY(!m.is_empty());
}
