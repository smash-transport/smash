/*
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "unittest.h"  // This include has to be first

#include "histogram.h"
#include "setup.h"

#include "../include/smash/decayaction.h"
#include "../include/smash/decaymodes.h"
#include "../include/smash/formfactors.h"
#include "../include/smash/kinematics.h"
#include "../include/smash/particletype.h"

using namespace smash;

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
      "π⁰ 0.138 0.0 - 111\n"
      "π⁺ 0.138 0.0 - 211\n"
      "ρ⁰ 0.776 0.149 - 113\n"
      "ρ⁺ 0.776 0.149 - 213\n"
      "ω 0.783 0.0085 - 223\n"
      "e⁻ 0.000511 0.0 + 11\n");
}

TEST(init_decay_modes) {
  DecayModes::load_decaymodes(
      "ρ \n"
      "1.      1  π π   \n"
      "4.72e-5 1  e⁻ e⁺ \n"
      "\n"
      "ω \n"
      "1.  1  π ρ \n");
}

TEST(omega_decay) {
  // set up omega decay action
  const ParticleType &type_omega = ParticleType::find(0x223);
  ParticleData omega{type_omega};
  omega.set_4momentum(0.782,                     // pole mass
                      ThreeVector(0., 0., 0.));  // at rest
  const auto act = make_unique<DecayAction>(omega, 0.);
  const auto srts = omega.effective_mass();
  act->add_decays(type_omega.get_partial_widths_hadronic(srts));

  const double dm = 0.001;       // bin size
  Histogram1d hist_charged(dm);  // histogram for charged rhos
  Histogram1d hist_neutral(dm);  // histogram for neutral rhos

  // sample the final state
  const int N_samples = 1E6;
  printf("sampling ...\n");
  for (int i = 0; i < N_samples; i++) {
    act->generate_final_state();
    const ParticleList &fs = act->outgoing_particles();
    if (fs.size() != 2) {
      std::cout << "unexpected FS size: " << fs.size() << "\n";
    }
    const ParticleData *rho = nullptr;
    if (!fs[1].type().pdgcode().is_pion()) {
      rho = &fs[1];
    } else if (!fs[0].type().pdgcode().is_pion()) {
      rho = &fs[0];
    }
    double m = rho->effective_mass();
    if (rho->type().charge() == 0) {
      hist_neutral.add(m);
    } else {
      hist_charged.add(m);
    }
  }

  // test with the analytical function
  const ParticleType &type_rho_zero = ParticleType::find(0x113);  // rho0
  const ParticleType &type_rho_plus = ParticleType::find(0x213);  // rho+
  const ParticleType &type_pi = ParticleType::find(0x111);        // pi0
  const auto mass_stable = type_pi.mass();

  printf("testing ρ⁰ distribution ...\n");
  hist_neutral.test([&](double m) {
    double pcm = pCM(srts, mass_stable, m);
    return type_rho_zero.spectral_function(m) * pcm *
           blatt_weisskopf_sqr(pcm, 1);
  }
                    //,"masses_rho_neutral.dat"
  );

  printf("testing ρ⁺ distribution ...\n");
  hist_charged.test([&](double m) {
    double pcm = pCM(srts, mass_stable, m);
    return type_rho_plus.spectral_function(m) * pcm *
           blatt_weisskopf_sqr(pcm, 1);
  }
                    //,"masses_rho_charged.dat"
  );
}
