/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "setup.h"
#include "../include/hadgas_eos.h"
#include "../include/constants.h"

using namespace Smash;

TEST(td_simple_gas) {
  ParticleType::create_type_list(
    "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
    "π⁰ 0.138 7.7e-9     111\n"
    "K⁰ 0.494 0.0        311\n"
    "N⁺ 0.938 0.0       2212\n");
  // Note that antiparticles are also created!
  const double T   = 0.1;
  const double mub = 0.8;
  const double mus = 0.1;
  COMPARE_ABSOLUTE_ERROR(HadgasEos::hadgas_net_baryon_density(T, mub, mus), 0.144397953, 1.e-6);
  COMPARE_ABSOLUTE_ERROR(HadgasEos::hadgas_net_strange_density(T, mub, mus), 0.002152896, 1.e-6);
  COMPARE_ABSOLUTE_ERROR(HadgasEos::hadgas_density(T, mub, mus), 0.15637968, 1.e-6);
  COMPARE_ABSOLUTE_ERROR(HadgasEos::hadgas_pressure(T, mub, mus), 0.015637968, 1.e-6);
  COMPARE_ABSOLUTE_ERROR(HadgasEos::hadgas_energy_density(T, mub, mus), 0.16493153, 1.e-6);
}

TEST(mu_zero_net_strangeness) {
  const double mub = 0.6;
  const double T   = 0.05;
  const double mus = HadgasEos::mus_net_strangeness0(T, mub);
  const double ns  = HadgasEos::hadgas_net_strange_density(T, mub, mus);
  VERIFY(std::abs(ns) < 1.e-4) << ns << ", mus = " << mus;
}

TEST(solve_EoS_substitute) {
  const double mub = 0.2;
  const double mus = 0.0;
  const double T = 0.30;
  const double e  = HadgasEos::hadgas_energy_density(T, mub, mus);
  const double nb = HadgasEos::hadgas_net_baryon_density(T, mub, mus);
  const double ns = HadgasEos::hadgas_net_strange_density(T, mub, mus);
  std::unique_ptr<HadgasEos> eos = make_unique<HadgasEos>();
  const std::array<double, 3> sol = eos->solve_hadgas_eos(e, nb, ns);
  COMPARE_ABSOLUTE_ERROR(sol[0], T, 1.e-4);
  COMPARE_ABSOLUTE_ERROR(sol[1], mub, 1.e-4);
  COMPARE_ABSOLUTE_ERROR(sol[2], mus, 1.e-4);
}
