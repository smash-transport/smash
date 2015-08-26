/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "unittest.h"
#include "setup.h"
#include "histogram.h"

#include "../include/particletype.h"
#include "../include/decaymodes.h"
#include "../include/decayaction.h"

using namespace Smash;

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
      "π⁰ 0.138 0.0 111\n"
      "π⁺ 0.138 0.0 211\n"
      "ρ⁰ 0.776 0.149 113\n"
      "ρ⁺ 0.776 0.149 213\n"
      "ω 0.783 0.0085 223\n"
      "e⁻ 0.000511 0.0 11\n");
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
  const auto act = make_unique<DecayAction>(omega, 0.f);
  act->add_decays(type_omega.get_partial_widths_hadronic(omega.effective_mass()));

  const float dm = 0.001;  // bin size
  Histogram1d hist(dm);    // histogram

  // sample the final state threlve gazillion times
  const int N_samples = 1E7;
  for (int i = 0; i < N_samples; i++) {
    act->generate_final_state();
    const ParticleList &fs = act->outgoing_particles();
    if (fs.size() != 2) {
      std::cout << "unexpected FS size: " << fs.size() << "\n";
    }
    float m = 0.;
    if (!fs[1].type().pdgcode().is_pion()) {
      m = fs[1].effective_mass();
    } else if (!fs[0].type().pdgcode().is_pion()) {
      m = fs[0].effective_mass();
    }
    hist.add(m);
  }

  // print out the histogram
  hist.print_to_file("mass_sampling.dat");
}
