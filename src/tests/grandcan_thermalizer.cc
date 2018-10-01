/*
 *
 *    Copyright (c) 2016-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "unittest.h"  // This include has to be first

#include "setup.h"

#include "../include/smash/boxmodus.h"
#include "../include/smash/grandcan_thermalizer.h"
#include "../include/smash/logging.h"
#include "../include/smash/thermalizationaction.h"

using namespace smash;

TEST(init) { set_default_loglevel(einhard::INFO); }

TEST(create_part_list) {
  // Note that antiparticles are also created!
  /* At least one baryon and one strange particle have to be in the list.
     Otherwise nb or ns is always 0, which makes degenerate Jacobian in
     the system to be solved.
   */
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
      "N⁺ 0.938 0.0 + 2212\n"
      "K⁰ 0.494 0.0 -  311\n");
}

static BoxModus create_box_for_tests() {
  auto conf = Test::configuration();
  const int N = 30;
  const double L = 10.;
  const double T_init = 0.2;
  conf["Modus"] = "Box";
  conf["Modi"]["Box"]["Init_Multiplicities"]["2212"] = N;
  conf["Modi"]["Box"]["Init_Multiplicities"]["311"] = N;
  conf["Modi"]["Box"]["Length"] = L;
  conf["Modi"]["Box"]["Temperature"] = T_init;
  conf["Modi"]["Box"]["Initial_Condition"] = "thermal momenta";
  conf["Modi"]["Box"]["Start_Time"] = 0.0;
  const ExperimentParameters par = smash::Test::default_parameters();
  return BoxModus(conf["Modi"], par);
}

TEST(rest_frame_transformation) {
  //  1. Generate a box of particles
  //  2. Construct ThermLatticeNode from these particles
  //  3. Go to the rest frame
  //  4. Check Tmu0 = (e+p)umu unu - p gmunu
  //  5. Check that EoS is satisfied
  Particles P;
  const ExperimentParameters par = smash::Test::default_parameters();
  BoxModus b = create_box_for_tests();
  b.initial_conditions(&P, par);

  HadronGasEos eos = HadronGasEos(false);
  ThermLatticeNode node = ThermLatticeNode();
  const ThreeVector v_boost(0.1, 0.2, 0.8);
  const double L = b.length();
  for (auto &part : P) {
    part.boost(v_boost);
    node.add_particle(part, std::sqrt(1.0 - v_boost.sqr()) / (L * L * L));
  }
  node.compute_rest_frame_quantities(eos);

  // Tmu0 should satisfy ideal hydro form
  const double ep_gamm_sqr = (node.e() + node.p()) / (1.0 - node.v().sqr());
  COMPARE_ABSOLUTE_ERROR(node.Tmu0().x0(), ep_gamm_sqr - node.p(), 1.e-8);
  COMPARE_ABSOLUTE_ERROR(node.Tmu0().x1(), ep_gamm_sqr * node.v().x1(), 1.e-8);
  COMPARE_ABSOLUTE_ERROR(node.Tmu0().x2(), ep_gamm_sqr * node.v().x2(), 1.e-8);
  COMPARE_ABSOLUTE_ERROR(node.Tmu0().x3(), ep_gamm_sqr * node.v().x3(), 1.e-8);

  // EoS should be satisfied
  const double T = node.T();
  const double mub = node.mub();
  const double mus = node.mus();
  const double gamma = 1.0 / std::sqrt(1.0 - node.v().sqr());
  const double tolerance = 5.e-4;
  COMPARE_ABSOLUTE_ERROR(node.p(), eos.pressure(T, mub, mus), tolerance);
  COMPARE_ABSOLUTE_ERROR(node.e(), eos.energy_density(T, mub, mus), tolerance);
  COMPARE_ABSOLUTE_ERROR(node.nb(), eos.net_baryon_density(T, mub, mus) * gamma,
                         tolerance);
  COMPARE_ABSOLUTE_ERROR(
      node.ns(), eos.net_strange_density(T, mub, mus) * gamma, tolerance);
}

TEST(thermalization_action) {
  Particles P;
  BoxModus b = create_box_for_tests();
  const ExperimentParameters par = smash::Test::default_parameters();
  b.initial_conditions(&P, par);

  Configuration th_conf = Test::configuration();
  std::vector<int> cell_n = {1, 1, 1};
  th_conf["Cell_Number"] = cell_n;
  th_conf["Critical_Edens"] = 0.01;
  th_conf["Algorithm"] = "biased BF";
  th_conf["Start_Time"] = 0.0;
  th_conf["Timestep"] = 1.0;
  auto thermalizer = b.create_grandcan_thermalizer(th_conf);

  const DensityParameters dens_par = DensityParameters(par);
  std::cout << "Updating lattice" << std::endl;
  thermalizer->update_thermalizer_lattice(P, dens_par, true);
  std::cout << "Thermalizing" << std::endl;
  thermalizer->thermalize(P, 0.0, par.testparticles);
  ThermalizationAction th_act(*thermalizer, 0.0);
  // If all this did not crash - the test is passed.
}
