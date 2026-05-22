/*
 *
 *    Copyright (c) 2026
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include <fstream>

#include "histogram.h"
#include "setup.h"
#include "smash/particletype.h"

using namespace smash;

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
      "A 3.0 0.3 + 50661\n"
      "B 0.6 0.3 + 10661\n"
      "C 1.0 0.6 + 20661\n"
      "D 0.1 0.0 + 30661\n");

  const std::string decays_input(
      "A       \n"
      "0.8 0 B B C \n"
      "0.2 0 C D   \n"
      "B       \n"
      "0.5 0  D D  \n"
      "0.5 0  D D D\n"
      "C       \n"
      "0.3 0  B D  \n"
      "0.4 0  D D  \n"
      "0.3 0  D D D\n");
  DecayModes::load_decaymodes(decays_input);
  ParticleType::check_consistency();
}

TEST(phasespace_manybody) {
  // The test only passes if there is enough energy for the phase space
  constexpr double sqrts = 10;

  const ParticleType &A = ParticleType::find(0x50661);
  const ParticleType &B = ParticleType::find(0x10661);
  const ParticleType &C = ParticleType::find(0x20661);
  ParticleTypePtrList types = {&A, &B, &C};
  std::vector<FourVector> sampled_momenta(types.size());

  const int N_sample = 100000;
  const double dm_hist = 0.05;
  Histogram1d hist_A(dm_hist);
  Histogram1d hist_B(dm_hist);
  Histogram1d hist_C(dm_hist);

  std::ostringstream buffer;
  for (int i = 0; i < N_sample; ++i) {
    smash::detail::sample_manybody_phasespace_impl(sqrts, types,
                                                   sampled_momenta);
    for (const auto &p : sampled_momenta) {
      buffer << p.abs() << " ";
    }
    buffer << "\n";

    hist_A.add(sampled_momenta[0].abs());
    hist_B.add(sampled_momenta[1].abs());
    hist_C.add(sampled_momenta[2].abs());
  }

  // Uncomment for printout to file
  /*
    std::ofstream out("manybody_histogram.dat");
    out << buffer.str();
    std::ofstream analytic("manybody_analytic.dat");
    for (double m = 0; m < 6; m += 0.02) {
      analytic << m
               << " " << A.full_spectral_function(m)
               << " " << B.full_spectral_function(m)
               << " " << C.full_spectral_function(m) << std::endl;
    }
  */

  hist_A.test([&](double m) { return A.full_spectral_function(m); });
  hist_B.test([&](double m) { return B.full_spectral_function(m); });
  hist_C.test([&](double m) { return C.full_spectral_function(m); });
}
