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
#include "../include/boxmodus.h"
#include "../include/grandcan_thermalizer.h"

using namespace Smash;

TEST(create_part_list) {
  // Note that antiparticles are also created!
  /* At least one baryon and one strange particle have to be in the list.
     Otherwise nb or ns is always 0, which makes degenerate Jacobian in
     the system to be solved.
   */
  ParticleType::create_type_list(
    "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
    "N⁺ 0.938 0.0       2212\n"
    "K⁰ 0.494 0.0        311\n");
}

TEST(rest_frame_transformation) {
  //  1. Generate a box of particles
  //  2. Construct ThermLatticeNode from these particles
  //  3. Go to the rest frame
  //  4. Check Tmu0 = (e+p)umu unu - p gmunu
  //  5. Check that EoS is satisfied
  auto conf = Test::configuration();
  const int N = 300;
  const float L = 10.0f;
  const float T_init = 0.2;
  const ThreeVector v_boost(0.1, 0.2, 0.8);
  conf["Modus"] = "Box";
  conf.take({"Modi", "Box", "Init_Multiplicities"});
  conf["Modi"]["Box"]["Init_Multiplicities"]["2212"] = N;
  conf["Modi"]["Box"]["Init_Multiplicities"]["311"] = N;
  conf["Modi"]["Box"]["Length"] = L;
  conf["Modi"]["Box"]["Temperature"] = T_init;
  const ExperimentParameters par = Smash::Test::default_parameters();
  std::unique_ptr<BoxModus> b = make_unique<BoxModus>(conf["Modi"], par);

  Particles P;
  b->initial_conditions(&P, par);

  HadronGasEos eos = HadronGasEos(false);
  ThermLatticeNode node = ThermLatticeNode();
  for (auto &part : P) {
    part.boost(v_boost);
    node.add_particle(part, std::sqrt(1.0 - v_boost.sqr())/(L*L*L));
  }
  node.compute_rest_frame_quantities(eos);

  // Tmu0 should satisfy ideal hydro form
  const double ep_gamm_sqr = (node.e()+node.p())/(1.0 - node.v().sqr());
  COMPARE_ABSOLUTE_ERROR(node.Tmu0().x0(), ep_gamm_sqr - node.p(), 1.e-8);
  COMPARE_ABSOLUTE_ERROR(node.Tmu0().x1(), ep_gamm_sqr*node.v().x1(), 1.e-8);
  COMPARE_ABSOLUTE_ERROR(node.Tmu0().x2(), ep_gamm_sqr*node.v().x2(), 1.e-8);
  COMPARE_ABSOLUTE_ERROR(node.Tmu0().x3(), ep_gamm_sqr*node.v().x3(), 1.e-8);

  // EoS should be satisfied
  const double T = node.T();
  const double mub = node.mub();
  const double mus = node.mus();
  const double gamma = 1.0/std::sqrt(1.0 - node.v().sqr());
  const double tolerance = 5.e-4;
  COMPARE_ABSOLUTE_ERROR(node.p(), eos.pressure(T, mub, mus), tolerance);
  COMPARE_ABSOLUTE_ERROR(node.e(), eos.energy_density(T, mub, mus), tolerance);
  COMPARE_ABSOLUTE_ERROR(node.nb(), eos.net_baryon_density(T, mub, mus)*gamma, tolerance);
  COMPARE_ABSOLUTE_ERROR(node.ns(), eos.net_strange_density(T, mub, mus)*gamma, tolerance);
}


