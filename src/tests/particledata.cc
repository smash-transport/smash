/*
 *
 *    Copyright (c) 2014-2018,2020,2022-2023
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "setup.h"
#include "smash/pdgcode.h"

using namespace smash;

TEST_CATCH(wrong_pion_mass, std::runtime_error) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
      "π⁰ 0.1350 0      -  111\n");
}

TEST_CATCH(wrong_omega_mass, std::runtime_error) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
      "ω  0.7827 0.0085 -  223\n");
}

TEST_CATCH(wrong_neutron_mass, std::runtime_error) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
      "N⁰ 0.9381  0      + 2112\n");
}

TEST_CATCH(wrong_delta_mass, std::runtime_error) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
      "Δ  1.2319  0.117  + 2224 2214 2114 1114\n");
}

TEST_CATCH(wrong_kaon_mass, std::runtime_error) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
      "K   0.493 0 - 311 321\n");
}

TEST_CATCH(wrong_deuteron_mass, std::runtime_error) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
      "d  1.8755 0 +  1000010020\n");
}

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
      "σ  0.123  1.2    +  661\n"
      "π⁰ 0.1380 0      -  111\n"
      "π⁺ 0.1380 0      -  211\n"
      "ρ  0.7755 0.149  -  113 213\n"
      "η  0.5479 1.3e-6 -  221\n"
      "ω  0.7830 0.0085 -  223\n"
      "N⁺ 0.938  0      + 2212\n"
      "N⁰ 0.938  0      + 2112\n"
      "Δ  1.232  0.117  + 2224 2214 2114 1114\n");
}

TEST(create_particledata_piplus) {
  PdgCode pdg = 0x211;
  ParticleData p{ParticleType::find(pdg)};

  COMPARE(p.id(), -1);
  COMPARE(p.pdgcode(), pdg);
  COMPARE(p.id_process(), 0u);
  COMPARE(p.momentum().x0(), 0.0);
  COMPARE(p.momentum().x1(), 0.0);
  COMPARE(p.momentum().x2(), 0.0);
  COMPARE(p.momentum().x3(), 0.0);
  COMPARE(p.position().x0(), 0.0);
  COMPARE(p.position().x1(), 0.0);
  COMPARE(p.position().x2(), 0.0);
  COMPARE(p.position().x3(), 0.0);

  p.set_id(2);
  COMPARE(p.id(), 2);
}

TEST(set_get) {
  PdgCode smashon = 0x661;
  ParticleData p = Test::smashon();
  p.set_id(4);
  COMPARE(p.id(), 4);
  COMPARE(p.pdgcode(), smashon);
  COMPARE(p.is_hadron(), smashon.is_hadron());
  p.set_history(3, 5, ProcessType::None, 1.2, ParticleList{});
  COMPARE(p.id_process(), 5u);
  COMPARE(p.get_history().collisions_per_particle, 3);
  COMPARE(p.get_history().time_last_collision, 1.2);
  p.set_history(4, 6, ProcessType::None, 2.5, ParticleList{});
  COMPARE(p.id_process(), 6u);
  COMPARE(p.get_history().collisions_per_particle, 4);
  COMPARE(p.get_history().time_last_collision, 2.5);
  FourVector m(1.0, 1.2, 1.4, 1.6);
  p.set_4momentum(m);
  COMPARE(p.momentum(), FourVector(1.0, 1.2, 1.4, 1.6));
  ThreeVector M(1.1, 1.3, 1.5);
  p.set_4momentum(1.0, M);
  COMPARE(p.momentum(), FourVector(sqrt(1.0 + M.sqr()), 1.1, 1.3, 1.5));
}

TEST(set_get2) {
  ParticleData p = Test::smashon(Test::Position{3.5, 3.6, 3.7, 2345.3});
  COMPARE(p.position().x0(), 3.5);
  COMPARE(p.position().x1(), 3.6);
  COMPARE(p.position().x2(), 3.7);
  COMPARE(p.position().x3(), 2345.3);
  ThreeVector M(2.1, 2.3, 2.5);
  p.set_4momentum(2.0, M.x1(), M.x2(), M.x3());
  COMPARE(p.momentum().x0(), sqrt(4.0 + M.sqr()));
  COMPARE(p.momentum().x1(), 2.1);
  COMPARE(p.momentum().x2(), 2.3);
  COMPARE(p.momentum().x3(), 2.5);
  ThreeVector v = p.velocity();
  COMPARE(v.x1(), 2.1 / sqrt(4.0 + M.sqr()));
  COMPARE(v.x2(), 2.3 / sqrt(4.0 + M.sqr()));
  COMPARE(v.x3(), 2.5 / sqrt(4.0 + M.sqr()));
  COMPARE_RELATIVE_ERROR(v.abs(), 0.8941469381, 1e-10);
  COMPARE_RELATIVE_ERROR(p.inverse_gamma(), 0.4477736628, 1e-10);
  p.boost(v);
  COMPARE_RELATIVE_ERROR(p.momentum().x0(), 2.0, 1e-15);
  COMPARE_ABSOLUTE_ERROR(p.momentum().x1(), 0.0, 1e-15);
  COMPARE_ABSOLUTE_ERROR(p.momentum().x2(), 0.0, 1e-15);
  COMPARE_ABSOLUTE_ERROR(p.momentum().x3(), 0.0, 1e-15);
  COMPARE_ABSOLUTE_ERROR(p.velocity().abs(), 0., 1e-15);
  COMPARE_RELATIVE_ERROR(p.inverse_gamma(), 1., 1e-15);
}

TEST(comparisons) {
  ParticleData p = Test::smashon(1);
  ParticleData q = Test::smashon(2);
  ParticleData r = Test::smashon(1);
  VERIFY(!(p == q));
  VERIFY(p == r);
  VERIFY(p == 1);
  VERIFY(p < 2);
  VERIFY(p < q);
}

TEST(translation) {
  ParticleData p = Test::smashon(Test::Position{0, 0, 0, 0});
  COMPARE(p.position(), FourVector(0, 0, 0, 0));
  COMPARE(p.translated({1, 0, 0}).position(), FourVector(0, 1, 0, 0));
  COMPARE(p.translated({1, 2, 0}).position(), FourVector(0, 1, 2, 0));
  COMPARE(p.translated({1, 2, 3}).position(), FourVector(0, 1, 2, 3));
}

TEST(parity) {
  const auto p = Parity::Pos;
  const auto n = Parity::Neg;
  COMPARE(-p, n);
  COMPARE(-n, p);
  COMPARE(n * p, n);
  COMPARE(p * n, n);
  COMPARE(p * p, p);
  COMPARE(n * n, p);
}

TEST(spin_projection) {
  for(int i=0; i<1000; i++){
    // Pi+
    ParticleData p1{ParticleType::find(smash::pdg::pi_p)};
    VERIFY(p1.spin_projection() == 0);
    // Proton
    ParticleData p2{ParticleType::find(smash::pdg::p)};
    auto s2 = p2.spin_projection();
    VERIFY(s2 == 1 || s2 == -1);
    // Delta-
    ParticleData p3{ParticleType::find(0x1114)};
    auto s3 = p3.spin_projection();
    VERIFY(s3 == -3 || s3 == -1 || s3 == 1 || s3 == 3);
  }
}

TEST(spin_flip) {
  // S=0
  ParticleData p1{ParticleType::find(smash::pdg::pi_p)};
  // S=3/2
  ParticleData p2{ParticleType::find(0x1114)};
  ParticleData p3{ParticleType::find(0x1114)};

  p2.set_spin_projection(-3);
  p3.set_spin_projection(1);

  // flip all spin projections
  p1.flip_spin_projection();
  p2.flip_spin_projection();
  p3.flip_spin_projection();

  VERIFY(
    p1.spin_projection() == 0 && 
    p2.spin_projection() == 3 &&
    p3.spin_projection() == -1  
  );
}