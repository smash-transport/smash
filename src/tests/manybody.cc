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
  constexpr double sqrts = 5;
  constexpr int samples = 1000;

  const ParticleType &A = ParticleType::find(0x50661);
  const ParticleType &B = ParticleType::find(0x10661);
  const ParticleType &C = ParticleType::find(0x20661);
  ParticleTypePtrList types = {&A, &B, &C};
  std::vector<FourVector> sampled_momenta(types.size());


  Action::sample_manybody_phasespace_impl(sqrts, types, sampled_momenta);
  // Remove the comment for printout
  /*
  std::ofstream out("manybody_histogram.dat");
  for (int i = 0; i < samples; ++i) {
    Action::sample_manybody_phasespace_impl(sqrts, types, sampled_momenta);
    for (const auto &p : sampled_momenta) {
      out << p.abs() << " ";
    }
    out << "\n";
  }

  std::ofstream analytic("manybody_analytic.dat");
  for (double m = 0; m < 6; m += 0.02) {
    analytic << m
             << " " << A.spectral_function(m)
             << " " << B.spectral_function(m)
             << " " << C.spectral_function(m) << std::endl;
  }
  */
}
