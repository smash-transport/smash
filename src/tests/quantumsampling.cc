/*
 *
 *    Copyright (c) 2020-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "histogram.h"
#include "setup.h"

#include "../include/smash/chemicalpotential.h"
#include "../include/smash/constants.h"
#include "../include/smash/hadgas_eos.h"
#include "../include/smash/quantumsampling.h"

using namespace smash;

TEST(create_particles_table) {
  Test::create_actual_particletypes();
  Test::create_actual_decaymodes();
}

TEST(chemical_potential_solver) {
  const double g = 2.0,
    m = 0.938,
    T = 0.015,
    stat = 1.0,
    precision = 1.e-8,
    mu_expected = 1.1,
    n_computed_by_mathematica = 0.847003688570, // [fm^-3]
    hbarc3 = hbarc * hbarc * hbarc;
  const double mu = effective_chemical_potential(g, m,
						 n_computed_by_mathematica * hbarc3,  // convert density to GeV^3
						 T, stat, precision);
  COMPARE_ABSOLUTE_ERROR(mu, mu_expected, precision) <<
    "Expected: " << mu_expected << ", got: " << mu;
}

TEST(sample_fermi_distribution) {
  PdgCode proton(0x2212);
  const int number_of_protons = 50;
  std::map<PdgCode, int> init_mult = {{proton, number_of_protons}};
  const double V = 50.0,  // fm^3
    T = 0.01;  // GeV
  QuantumSampling sampler(init_mult, V, T); 

  const double dmomentum = 0.01;  // GeV
  Histogram1d hist(dmomentum);
  constexpr int N_TEST = 1E6;  // number of samples
  hist.populate(N_TEST, [&]() { return sampler.sample(proton); });
  const double g = proton.spin_degeneracy(),
    m = ParticleType::find(proton).mass(),
    n = number_of_protons / V * hbarc * hbarc * hbarc,
    stat = 1.0,
    mu = effective_chemical_potential(g, m, n, T, stat, 1.e-6);
  hist.test(
	    [&](double p) {
	      return p * p / (std::exp((std::sqrt(p*p + m*m) - mu) / T) + 1.0);
	    }, "quantum_sampling.dat");
}
