/*
 *
 *    Copyright (c) 2020,2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

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
  // precision (default 1.e-8) is used in solving for mu reproducing n;
  // comparison_precision is used to see if the value of mu agrees with the one
  // obtained from Mathematica; with precision = 1.e-8, comparison_precision at
  // 1.e-6 is a good enough test across different points of the phase diagram;
  // reproducing n is tested in the effective_chemical_potential itself.
  const double precision = 1.e-8, hbarc3 = hbarc * hbarc * hbarc;
  const double comparison_precision = 1.e-6;
  // parameters of a hypothetical particle, subject to either
  // Fermi, Bose, or Boltzmann distribution
  const double g = 4.0, mass = 0.938;

  // initialize the ChemicalPotentialSolver object
  ChemicalPotentialSolver mu_solver;

  // 1) Fermi distribution
  // tested T = 0.001 [GeV]
  // tested nB = 0.0125 * 0.160fm^-3
  const double stat1 = 1.0, n1 = 0.002, T1 = 0.001,
               mu_computed_by_Mathematica_1 = 0.939464245179786;
  const double mu1 = mu_solver.effective_chemical_potential(
      g, mass, n1 * hbarc3,  // convert density to GeV^3
      T1, stat1, precision);
  COMPARE_ABSOLUTE_ERROR(mu1, mu_computed_by_Mathematica_1,
                         comparison_precision)
      << "Expected: " << mu_computed_by_Mathematica_1 << ", got: " << mu1;

  // 2) Bose distribution
  // tested T = 0.001 [GeV]
  // tested nB = 0.0125 * 0.160fm^-3
  const double stat2 = -1.0, n2 = 0.002, T2 = 0.001,
               mu_computed_by_Mathematica_2 = 0.937976538402567;
  const double mu2 = mu_solver.effective_chemical_potential(
      g, mass, n2 * hbarc3,  // convert density to GeV^3
      T2, stat2, precision);
  COMPARE_ABSOLUTE_ERROR(mu2, mu_computed_by_Mathematica_2,
                         comparison_precision)
      << "Expected: " << mu_computed_by_Mathematica_2 << ", got: " << mu2;

  // 3) Fermi distribution
  // tested T = 0.001 [GeV]
  // tested nB = 1 * 0.160fm^-3
  const double stat3 = 1.0, n3 = 0.160, T3 = 0.001,
               mu_computed_by_Mathematica_3 = 0.974159174500194;
  const double mu3 = mu_solver.effective_chemical_potential(
      g, mass, n3 * hbarc3,  // convert density to GeV^3
      T3, stat3, precision);
  COMPARE_ABSOLUTE_ERROR(mu3, mu_computed_by_Mathematica_3,
                         comparison_precision)
      << "Expected: " << mu_computed_by_Mathematica_3 << ", got: " << mu3;

  // 4) Fermi distribution
  // tested T = 0.001 [GeV]
  // tested nB = 5 * 0.160fm^-3
  const double stat4 = 1.0, n4 = 0.8, T4 = 0.001,
               mu_computed_by_Mathematica_4 = 1.04025839770074;
  const double mu4 = mu_solver.effective_chemical_potential(
      g, mass, n4 * hbarc3,  // convert density to GeV^3
      T4, stat4, precision);
  COMPARE_ABSOLUTE_ERROR(mu4, mu_computed_by_Mathematica_4,
                         comparison_precision)
      << "Expected: " << mu_computed_by_Mathematica_4 << ", got: " << mu4;

  // 5) Fermi distribution
  // tested T = 0.016 [GeV]
  // tested nB = 1.0 * 0.160fm^-3
  const double stat5 = 1.0, n5 = 0.160, T5 = 0.016,
               mu_computed_by_Mathematica_5 = 0.966539703510592;
  const double mu5 = mu_solver.effective_chemical_potential(
      g, mass, n5 * hbarc3,  // convert density to GeV^3
      T5, stat5, precision);
  COMPARE_ABSOLUTE_ERROR(mu5, mu_computed_by_Mathematica_5,
                         comparison_precision)
      << "Expected: " << mu_computed_by_Mathematica_5 << ", got: " << mu5;

  // 6) Bose distribution
  // tested T = 0.016 [GeV]
  // tested nB = 1.0 * 0.160fm^-3
  const double stat6 = -1.0, n6 = 0.160, T6 = 0.016,
               mu_computed_by_Mathematica_6 = 0.937999333245918;
  const double mu6 = mu_solver.effective_chemical_potential(
      g, mass, n6 * hbarc3,  // convert density to GeV^3
      T6, stat6, precision);
  COMPARE_ABSOLUTE_ERROR(mu6, mu_computed_by_Mathematica_6,
                         comparison_precision)
      << "Expected: " << mu_computed_by_Mathematica_6 << ", got: " << mu6;

  // 7) Fermi distribution
  // tested T = 0.150 [GeV]
  // tested nB = 2.0 * 0.160fm^-3
  const double stat7 = 1.0, n7 = 0.320, T7 = 0.150,
               mu_computed_by_Mathematica_7 = 0.648402961261576;
  const double mu7 = mu_solver.effective_chemical_potential(
      g, mass, n7 * hbarc3,  // convert density to GeV^3
      T7, stat7, precision);
  COMPARE_ABSOLUTE_ERROR(mu7, mu_computed_by_Mathematica_7,
                         comparison_precision)
      << "Expected: " << mu_computed_by_Mathematica_7 << ", got: " << mu5;

  // 8) Bose distribution
  // tested T = 0.150 [GeV]
  // tested nB = 2.0 * 0.160fm^-3
  const double stat8 = -1.0, n8 = 0.32, T8 = 0.150,
               mu_computed_by_Mathematica_8 = 0.635499008702213;
  const double mu8 = mu_solver.effective_chemical_potential(
      g, mass, n8 * hbarc3,  // convert density to GeV^3
      T8, stat8, precision);
  COMPARE_ABSOLUTE_ERROR(mu8, mu_computed_by_Mathematica_8,
                         comparison_precision)
      << "Expected: " << mu_computed_by_Mathematica_8 << ", got: " << mu8;

  // 9) Boltzmann distribution
  // tested T = 0.150 [GeV]
  // tested nB = 2.0 * 0.160fm^-3
  const double stat9 = 0.0, n9 = 0.32, T9 = 0.150,
               mu_computed_by_Mathematica_9 = 0.642000594923801;
  const double mu9 = mu_solver.effective_chemical_potential(
      g, mass, n9 * hbarc3,  // convert density to GeV^3
      T9, stat9, precision);
  COMPARE_ABSOLUTE_ERROR(mu9, mu_computed_by_Mathematica_9,
                         comparison_precision)
      << "Expected: " << mu_computed_by_Mathematica_9 << ", got: " << mu9;

  // 10) Fermi distribution
  // tested T = 0.200 [GeV]
  // tested nB = 0.1 * 0.160fm^-3
  const double stat10 = 1.0, n10 = 0.016, T10 = 0.200,
               mu_computed_by_Mathematica_10 = -0.158561260539330;
  const double mu10 = mu_solver.effective_chemical_potential(
      g, mass, n10 * hbarc3,  // convert density to GeV^3
      T10, stat10, precision);
  COMPARE_ABSOLUTE_ERROR(mu10, mu_computed_by_Mathematica_10,
                         comparison_precision)
      << "Expected: " << mu_computed_by_Mathematica_10 << ", got: " << mu10;
}

TEST(sample_fermi_distribution) {
  PdgCode proton(0x2212);
  const int number_of_protons = 50;
  std::map<PdgCode, int> init_mult = {{proton, number_of_protons}};
  const double V = 50.0,  // fm^3
      T = 0.01;           // GeV
  // initialize the ChemicalPotentialSolver object and sampler
  ChemicalPotentialSolver mu_solver;
  QuantumSampling sampler(init_mult, V, T);

  const double dmomentum = 0.01;  // GeV
  Histogram1d hist(dmomentum);
  constexpr int N_TEST = 1E6;  // number of samples
  hist.populate(N_TEST, [&]() { return sampler.sample(proton); });
  const double g = proton.spin_degeneracy(),
               m = ParticleType::find(proton).mass(),
               n = number_of_protons / V * hbarc * hbarc * hbarc, stat = 1.0,
               mu = mu_solver.effective_chemical_potential(g, m, n, T, stat,
                                                           1.e-6);
  hist.test(
      [&](double p) {
        return p * p / (std::exp((std::sqrt(p * p + m * m) - mu) / T) + 1.0);
      },
      "quantum_sampling.dat");
}
