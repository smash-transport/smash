/*
 *
 *    Copyright (c) 2016-2020,2022,2024-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/hadgas_eos.h"

#include "setup.h"
#include "smash/constants.h"

using namespace smash;

TEST(create_particles_table) {
  Test::create_actual_particletypes();
  Test::create_actual_decaymodes();
}

TEST(td_simple_gas) {
  const double T = 0.1;
  const double mub = 0.8;
  const double mus = 0.1;
  const double muq = 0.05;
  // The EOS (Equation of State) is sensitive to resonance parameters.
  // If decaymodes.txt or particles.txt are changed, these tests might fail.
  // The current values are from running SMASH with the parameters in the input
  // folder.
  const double net_baryon_density =
      HadronGasEos::net_baryon_density(T, mub, mus, muq);
  const double net_strange_density =
      HadronGasEos::net_strange_density(T, mub, mus, muq);
  const double net_charge_density =
      HadronGasEos::net_charge_density(T, mub, mus, muq);
  const double density = HadronGasEos::density(T, mub, mus, muq);
  const double pressure = HadronGasEos::pressure(T, mub, mus, muq);
  const double energy_density = HadronGasEos::energy_density(T, mub, mus, muq);
  std::stringstream expected_values;
  expected_values << std::setprecision(10);
  expected_values
      << "\nExpected values:\n net_baryon_density net_strange_density "
         "net_charge_density density pressure energy_density\n";
  expected_values << net_baryon_density << " " << net_strange_density << " "
                  << net_charge_density << " " << density << " " << pressure
                  << " " << energy_density << " " << std::endl;

  COMPARE_ABSOLUTE_ERROR(net_baryon_density, 0.5832753159, 1.e-6)
      .on_failure(expected_values.str());
  COMPARE_ABSOLUTE_ERROR(net_strange_density, -0.03686360221, 1.e-6)
      .on_failure(expected_values.str());
  COMPARE_ABSOLUTE_ERROR(net_charge_density, 0.4236041139, 1.e-6)
      .on_failure(expected_values.str());
  COMPARE_ABSOLUTE_ERROR(density, 0.6232158477, 1.e-6)
      .on_failure(expected_values.str());
  COMPARE_ABSOLUTE_ERROR(pressure, 0.06232158477, 1.e-6)
      .on_failure(expected_values.str());
  COMPARE_ABSOLUTE_ERROR(energy_density, 0.7325009276, 1.e-6)
      .on_failure(expected_values.str());
}

TEST(mu_zero_net_strangeness) {
  const double mub = 0.6;
  const double muq = 0.1;
  const double T = 0.05;
  const double mus = HadronGasEos::mus_net_strangeness0(T, mub, muq);
  const double ns = HadronGasEos::net_strange_density(T, mub, mus, muq);
  VERIFY(std::abs(ns) < 1.e-4) << ns << ", mus = " << mus;
}

TEST(solve_EoS_substitute) {
  const double mub = 0.2;
  const double muq = 0.1;
  const double mus = 0.0;
  const double T = 0.30;
  const double e = HadronGasEos::energy_density(T, mub, mus, muq);
  const double nb = HadronGasEos::net_baryon_density(T, mub, mus, muq);
  const double ns = HadronGasEos::net_strange_density(T, mub, mus, muq);
  const double nq = HadronGasEos::net_charge_density(T, mub, mus, muq);
  HadronGasEos eos = HadronGasEos(false, false);
  const std::array<double, 4> sol = eos.solve_eos(e, nb, ns, nq);
  COMPARE_ABSOLUTE_ERROR(sol[0], T, 1.e-4);
  COMPARE_ABSOLUTE_ERROR(sol[1], mub, 1.e-4);
  COMPARE_ABSOLUTE_ERROR(sol[2], mus, 1.e-4);
  COMPARE_ABSOLUTE_ERROR(sol[3], muq, 1.e-4);
}

TEST(EoS_table) {
  // make a small table of EoS
  HadronGasEos eos = HadronGasEos(false, false);
  EosTable table = EosTable(0.1, 0.05, 0.05, 5, 5, 5);
  table.compile_table(eos, "small_test_table_eos.dat");
  EosTable::table_element x;
  const double my_e = 0.39, my_nb = 0.09, my_nq = 0.06;
  table.get(x, my_e, my_nb, my_nq);
  // check if tabulated values are the right solutions
  COMPARE_ABSOLUTE_ERROR(HadronGasEos::energy_density(x.T, x.mub, x.mus, x.muq),
                         my_e, 1.e-2);
  COMPARE_ABSOLUTE_ERROR(
      HadronGasEos::net_baryon_density(x.T, x.mub, x.mus, x.muq), my_nb, 1.e-2);
  COMPARE_ABSOLUTE_ERROR(
      HadronGasEos::net_charge_density(x.T, x.mub, x.mus, x.muq), my_nq, 1.e-2);
  remove("small_test_table_eos.dat");
}

/*
TEST(make_test_table) {
  // To switch on these tests, comment out the previous ones.
  // Otherwise the particle table confilct emerges.
  Test::create_actual_particletypes();
  const double ns = 0.0;
  const double nb = 0.3;
  HadronGasEos eos = HadronGasEos(false, false);
  for (int ie = 0; ie < 1000; ie++) {
    const double e = nb + 0.001 + 0.001 * ie;
    const std::array<double, 3> sol = eos.solve_eos(e, nb, ns);
    std::cout << e << " " << HadronGasEos::pressure(sol[0], sol[1], sol[2]) <<
std::endl;
  }
}

TEST(make_test_table2) {
  const double mub = 0.0;
  const double mus = 0.0;
  for (int it = 0; it < 1000; it++) {
    const double T = 0.070 + 0.001 * it;
    const double e  = HadronGasEos::energy_density(T, mub, mus);
    const double p  = HadronGasEos::pressure(T, mub, mus);
    std::cout << T << " " << (e-3*p)/(T*T*T*T) << std::endl;
  }
}*/
