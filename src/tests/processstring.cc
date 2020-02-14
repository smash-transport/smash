/*
 *
 *    Copyright (c) 2015-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "histogram.h"
#include "setup.h"

#include "../include/smash/angles.h"
#include "../include/smash/random.h"
#include "../include/smash/scatteraction.h"
#include "Pythia8/Pythia.h"

#include <iostream>

using namespace smash;
using smash::Test::Momentum;
using smash::Test::Position;

// Test needs to be here in order for ParticleType objects to work
TEST(init_particle_types) {
  Test::create_actual_particletypes();
  Test::create_actual_decaymodes();
  ParticleType::check_consistency();
}

template <typename Chi, typename Analytical>
void test_distribution(int n_test, double dx, Chi get_chi,
                       Analytical get_analyt) {
  Histogram1d hist(dx);
  // sample distribution and populate histogram
  hist.populate(n_test, get_chi);
  // test with the analytical function
  hist.test(get_analyt);
}

using namespace smash;

/**
 * Compare sampled values of the LUND function to analytical ones via:
 * \f$ f(z) = \frac{1}{z} (1 - z)^a \exp{ \left(- \frac{b m_T^2}{z} \right) }
 * Using the same framework as the tests in random.cc do.
 */
TEST(string_zlund) {
  test_distribution(1e7, 0.0001,
                    []() { return StringProcess::sample_zLund(1, 1, 1); },
                    [](double x) { return 1 / x * (1. - x) * exp(-1. / x); });
}

TEST(string_incoming_lightcone_momenta) {
  std::unique_ptr<StringProcess> sp =
      make_unique<StringProcess>(1.0, 1.0, .0, 0.001, .0, .0, 1., 1., .0, .0,
                                 .5, .0, .0, .0, .0, true, 1. / 3., true, 0.);

  ParticleData a{ParticleType::find(0x2212)};
  a.set_4momentum(0.938, 0., 0., 1.);
  FourVector p_a = a.momentum();

  ParticleData b{ParticleType::find(0x2212)};
  b.set_4momentum(0.938, 0., 0., -1.);
  FourVector p_b = b.momentum();

  sp->init({a, b}, 0.);

  double Ap = sp->getPPosA();
  double An = sp->getPNegA();
  double Bp = sp->getPPosB();
  double Bn = sp->getPnegB();

  // longitudinal direction is +z so p± is (E±pz) / sqrt(2)
  FUZZY_COMPARE(Ap, (p_a.x0() + p_a.x3()) * M_SQRT1_2);
  FUZZY_COMPARE(An, (p_a.x0() - p_a.x3()) * M_SQRT1_2);
  FUZZY_COMPARE(Bp, (p_b.x0() + p_b.x3()) * M_SQRT1_2);
  FUZZY_COMPARE(Bn, (p_b.x0() - p_b.x3()) * M_SQRT1_2);
}

TEST(string_lightcone_final_two) {
  double a = .0;
  double b = .0;
  double c = .0;
  double d = .0;
  std::unique_ptr<StringProcess> sp = make_unique<StringProcess>(
      1.0, 1.0, 0.5, 0.001, 1.0, 2.5, 0.217, 0.081, 0.7, 0.7, 0.25, 0.68, 0.98,
      0.25, 1.0, true, 1. / 3., true, 0.2);

  // returns false because mTsqr_string < 0.
  VERIFY(sp->make_lightcone_final_two(false, -1., 1., .0, .0, a, b, c, d) ==
         false);
  // returns false because mTrn_string < mTrn_had_forward + mTrn_had_backward
  VERIFY(sp->make_lightcone_final_two(false, .0, .0, 1., 1., a, b, c, d) ==
         false);
  // returns false because lambda_sqr == 0.
  VERIFY(sp->make_lightcone_final_two(false, 1., 1. / 2., -0.337106, 0.837106,
                                      a, b, c, d) == false);
}

TEST(string_find_leading) {
  int nq1 = 1;
  int nq2 = -1;
  ParticleData a{ParticleType::find(-0x2212)};  // anti proton
  a.set_id(0);
  ParticleData b{ParticleType::find(0x2212)};  // proton
  b.set_id(1);
  // Create outgoing particle list
  ParticleList outgoing = {a, b};
  std::pair<int, int> indices = StringProcess::find_leading(nq1, nq2, outgoing);
  // Indices should switch because of the signs of nq1 and nq2
  COMPARE(indices.first, 1);
  COMPARE(indices.second, 0);
}

TEST(string_scaling_factor) {
  int nquark = 1;
  // pion+ particle data
  ParticleData data_pip{ParticleType::find(0x211)};
  double suppression_factor = 0.5;
  StringProcess::assign_scaling_factor(nquark, data_pip, suppression_factor);
  COMPARE_ABSOLUTE_ERROR(data_pip.initial_xsec_scaling_factor(), 0.25, 1e-7);

  nquark = 2;
  // pion+ particle data
  ParticleData data_proton{ParticleType::find(0x2212)};
  suppression_factor = 3. / 2.;
  StringProcess::assign_scaling_factor(nquark, data_proton, suppression_factor);
  COMPARE_ABSOLUTE_ERROR(data_proton.initial_xsec_scaling_factor(), 1.0, 1e-7);
}

TEST(string_all_scaling_factors) {
  int baryon_string = 0;
  ParticleData a{ParticleType::find(-0x2212)};  // anti proton
  a.set_id(0);
  a.set_4momentum(0.938, 0., 0., 1.);
  ParticleData b{ParticleType::find(0x2212)};  // proton
  b.set_id(1);
  b.set_4momentum(0.938, 0., 0., 0.5);
  ParticleData c{ParticleType::find(0x2212)};  // proton
  c.set_id(2);
  c.set_4momentum(0.938, 0., 0., 0.75);
  // Create outgoing particle list
  ParticleList outgoing = {a, b, c};
  // Zero velocity so no sorting happens
  ThreeVector velvec{0., 0., 0.};
  double suppression_factor = 3.0;
  StringProcess::assign_all_scaling_factors(baryon_string, outgoing, velvec,
                                            suppression_factor);
  COMPARE_ABSOLUTE_ERROR(outgoing[0].initial_xsec_scaling_factor(), 1.0, 1e-7);
  COMPARE_ABSOLUTE_ERROR(outgoing[1].initial_xsec_scaling_factor(), .0, 1e-7);
  COMPARE_ABSOLUTE_ERROR(outgoing[2].initial_xsec_scaling_factor(), 1.0, 1e-7);
}

TEST(string_orthonormal_basis) {
  ThreeVector evec_polar = ThreeVector(0., 1., 0.);
  std::array<ThreeVector, 3> evec_basis;
  StringProcess::make_orthonormal_basis(evec_polar, evec_basis);

  VERIFY(std::abs(evec_basis[0].x1()) < really_small);
  VERIFY(std::abs(evec_basis[0].x2() - 1.) < really_small);
  VERIFY(std::abs(evec_basis[0].x3()) < really_small);

  VERIFY(std::abs(evec_basis[1].x1()) < really_small);
  VERIFY(std::abs(evec_basis[1].x2()) < really_small);
  VERIFY(std::abs(evec_basis[1].x3() + 1.) < really_small);

  VERIFY(std::abs(evec_basis[2].x1() + 1.) < really_small);
  VERIFY(std::abs(evec_basis[2].x2()) < really_small);
  VERIFY(std::abs(evec_basis[2].x3()) < really_small);
}

TEST(string_find_excess_constituent) {
  std::array<int, 5> excess_quark;
  std::array<int, 5> excess_antiq;

  PdgCode pdg_piplus = PdgCode(0x211);
  PdgCode pdg_Kplus = PdgCode(0x321);
  StringProcess::find_excess_constituent(pdg_Kplus, pdg_piplus, excess_quark,
                                         excess_antiq);
  VERIFY(excess_quark[0] == 0);
  VERIFY(excess_quark[1] == 0);
  VERIFY(excess_quark[2] == 0);
  VERIFY(excess_quark[3] == 0);
  VERIFY(excess_quark[4] == 0);
  VERIFY(excess_antiq[0] == -1);
  VERIFY(excess_antiq[1] == 0);
  VERIFY(excess_antiq[2] == 1);
  VERIFY(excess_antiq[3] == 0);
  VERIFY(excess_antiq[4] == 0);

  PdgCode pdg_neutron = PdgCode(0x2112);
  PdgCode pdg_Omega = PdgCode(0x3334);
  StringProcess::find_excess_constituent(pdg_Omega, pdg_neutron, excess_quark,
                                         excess_antiq);
  VERIFY(excess_quark[0] == -2);
  VERIFY(excess_quark[1] == -1);
  VERIFY(excess_quark[2] == 3);
  VERIFY(excess_quark[3] == 0);
  VERIFY(excess_quark[4] == 0);
  VERIFY(excess_antiq[0] == 0);
  VERIFY(excess_antiq[1] == 0);
  VERIFY(excess_antiq[2] == 0);
  VERIFY(excess_antiq[3] == 0);
  VERIFY(excess_antiq[4] == 0);

  PdgCode pdg_anti_neutron = PdgCode(-0x2112);
  PdgCode pdg_anti_Xi0 = PdgCode(-0x3322);
  StringProcess::find_excess_constituent(pdg_anti_Xi0, pdg_anti_neutron,
                                         excess_quark, excess_antiq);
  VERIFY(excess_quark[0] == 0);
  VERIFY(excess_quark[1] == 0);
  VERIFY(excess_quark[2] == 0);
  VERIFY(excess_quark[3] == 0);
  VERIFY(excess_quark[4] == 0);
  VERIFY(excess_antiq[0] == -2);
  VERIFY(excess_antiq[1] == 0);
  VERIFY(excess_antiq[2] == 2);
  VERIFY(excess_antiq[3] == 0);
  VERIFY(excess_antiq[4] == 0);
}

TEST(string_quarks_from_diquark) {
  int id_diquark;
  int id1, id2, deg_spin;

  id1 = 0;
  id2 = 0;
  deg_spin = 0;
  // ud-diquark
  id_diquark = 2103;
  StringProcess::quarks_from_diquark(id_diquark, id1, id2, deg_spin);
  VERIFY(id1 == 2);
  VERIFY(id2 == 1);
  VERIFY(deg_spin == 3);

  id1 = 0;
  id2 = 0;
  deg_spin = 0;
  // ud-antidiquark
  id_diquark = -2101;
  StringProcess::quarks_from_diquark(id_diquark, id1, id2, deg_spin);
  VERIFY(id1 == -2);
  VERIFY(id2 == -1);
  VERIFY(deg_spin == 1);
}

TEST(string_diquark_from_quarks) {
  // ud-diquark
  int id1 = 1;
  int id2 = 2;
  int id_diquark = StringProcess::diquark_from_quarks(id1, id2);
  VERIFY(id_diquark == 2101 || id_diquark == 2103);
  // uu-diquark
  id1 = 2;
  id_diquark = StringProcess::diquark_from_quarks(id1, id2);
  VERIFY(id_diquark == 2203);
}

TEST(string_make_string_ends) {
  int id1, id2;

  // decompose neutron in d, ud-diquark or u, dd-diquark
  PdgCode pdg_neutron = PdgCode(0x2112);
  StringProcess::make_string_ends(pdg_neutron, id1, id2, 1. / 3.);
  VERIFY((id1 == 2 && (id2 == 1103)) || (id1 == 1 && (id2 == 2103 || 2101)));

  // decompose anti-neutron in dbar, ud-antidiquark or ubar, dd-antidiquark
  PdgCode pdg_anti_neutron = PdgCode(-0x2112);
  StringProcess::make_string_ends(pdg_anti_neutron, id1, id2, 1. / 3.);
  VERIFY((id2 == -2 && (id1 == -1103)) ||
         (id2 == -1 && (id1 == -2103 || -2101)));

  // decompose pion0 into u, ubar or d, dbar
  PdgCode pdg_pizero = PdgCode(0x111);
  StringProcess::make_string_ends(pdg_pizero, id1, id2, 1. / 3.);
  VERIFY((id1 == 1 && id2 == -1) || (id1 == 2 && id2 == -2));

  // decompose pion+ into u, dbar
  PdgCode pdg_piplus = PdgCode(0x211);
  StringProcess::make_string_ends(pdg_piplus, id1, id2, 1. / 3.);
  VERIFY(id1 == 2 && id2 == -1);

  // decompose pion- into d, ubar
  PdgCode pdg_piminus = PdgCode(-0x211);
  StringProcess::make_string_ends(pdg_piminus, id1, id2, 1. / 3.);
  VERIFY(id1 == 1 && id2 == -2);

  // decompose proton into u, ud-diquark or d, uu-diquark
  PdgCode pdg_proton = PdgCode(0x2212);
  StringProcess::make_string_ends(pdg_proton, id1, id2, 1. / 3.);
  VERIFY((id1 == 1 && id2 == 2203) || (id1 == 2 && (id2 == 2101 || 2103)));

  // decompose anti-proton ubar, ud-antidiquark or dbar, uu-antidiquark
  PdgCode pdg_antip = PdgCode(-0x2212);
  StringProcess::make_string_ends(pdg_antip, id1, id2, 1. / 3.);
  VERIFY((id2 == -1 && id1 == -2203) || (id2 == -2 && (id1 == -2101 || -2103)));
}

TEST(string_set_Vec4) {
  // make arbitrary lightlike 4-vector with random direction
  Angles angle_random = Angles(0., 0.);
  angle_random.distribute_isotropically();
  const double energy = 10.;
  const ThreeVector mom = energy * angle_random.threevec();
  Pythia8::Vec4 vector = Pythia8::Vec4(0., 0., 0., 0.);
  // set Pythia8::Vec4
  vector = StringProcess::set_Vec4(energy, mom);
  // check if Pythia8::Vec4 is same with 4-vector from energy and mom
  const double energy_scale = 0.5 * (vector.e() + energy);
  VERIFY(std::abs(vector.e() - energy) < really_small * energy_scale);
  VERIFY(std::abs(vector.px() - mom.x1()) < really_small * energy_scale);
  VERIFY(std::abs(vector.py() - mom.x2()) < really_small * energy_scale);
  VERIFY(std::abs(vector.pz() - mom.x3()) < really_small * energy_scale);
}

TEST(pdg_map_for_pythia) {
  int pdgid_mapped = 0;

  // pi+ is mapped onto pi+
  PdgCode pdg_piplus = PdgCode(0x211);
  pdgid_mapped = StringProcess::pdg_map_for_pythia(pdg_piplus);
  VERIFY(pdgid_mapped == 211);

  // pi0 is mapped onto pi+
  PdgCode pdg_pi0 = PdgCode(0x111);
  pdgid_mapped = StringProcess::pdg_map_for_pythia(pdg_pi0);
  VERIFY(pdgid_mapped == 211);

  // pi- is mapped onto pi-
  PdgCode pdg_piminus = PdgCode(-0x211);
  pdgid_mapped = StringProcess::pdg_map_for_pythia(pdg_piminus);
  VERIFY(pdgid_mapped == -211);

  // proton is mapped onto proton
  PdgCode pdg_proton = PdgCode(0x2212);
  pdgid_mapped = StringProcess::pdg_map_for_pythia(pdg_proton);
  VERIFY(pdgid_mapped == 2212);

  // neutron is mapped onto neutron
  PdgCode pdg_neutron = PdgCode(0x2112);
  pdgid_mapped = StringProcess::pdg_map_for_pythia(pdg_neutron);
  VERIFY(pdgid_mapped == 2112);

  // antiproton is mapped onto antiproton
  PdgCode pdg_antip = PdgCode(-0x2212);
  pdgid_mapped = StringProcess::pdg_map_for_pythia(pdg_antip);
  VERIFY(pdgid_mapped == -2212);

  // K+ is mapped onto pi+
  PdgCode pdg_Kplus = PdgCode(0x321);
  pdgid_mapped = StringProcess::pdg_map_for_pythia(pdg_Kplus);
  VERIFY(pdgid_mapped == 211);

  // K- is mapped onto pi-
  PdgCode pdg_Kminus = PdgCode(-0x321);
  pdgid_mapped = StringProcess::pdg_map_for_pythia(pdg_Kminus);
  VERIFY(pdgid_mapped == -211);

  // Lambda is mapped onto neutron
  PdgCode pdg_Lambda = PdgCode(0x3122);
  pdgid_mapped = StringProcess::pdg_map_for_pythia(pdg_Lambda);
  VERIFY(pdgid_mapped == 2112);

  // anti-Lambda is mapped onto anti-neutron
  PdgCode pdg_antiL = PdgCode(-0x3122);
  pdgid_mapped = StringProcess::pdg_map_for_pythia(pdg_antiL);
  VERIFY(pdgid_mapped == -2112);
}

TEST(string_scaling_factors) {
  ParticleData a{ParticleType::find(0x2212)};
  ParticleData b{ParticleType::find(0x2212)};
  ParticleList incoming{a, b};
  ParticleData c{ParticleType::find(-0x2212)};  // anti proton
  ParticleData d{ParticleType::find(0x2212)};   // proton
  ParticleData e{ParticleType::find(0x111)};    // pi0
  ParticleData f{ParticleType::find(0x111)};    // pi0
  c.set_id(0);
  d.set_id(1);
  e.set_id(2);
  f.set_id(3);
  c.set_4momentum(0.938, 0., 0., 1.);
  d.set_4momentum(0.938, 0., 0., 0.5);
  e.set_4momentum(0.138, 0., 0., -0.5);
  f.set_4momentum(0.138, 0., 0., -1.);
  ParticleList outgoing = {e, d, c, f};  // here in random order
  constexpr double coherence_factor = 0.7;
  ThreeVector evec_coll = ThreeVector(0., 0., 1.);
  int baryon_string =
      incoming[random::uniform_int(0, 1)].type().baryon_number();
  StringProcess::assign_all_scaling_factors(baryon_string, outgoing, evec_coll,
                                            coherence_factor);
  // outgoing list is now assumed to be sorted by z-velocity (so c,d,e,f)
  VERIFY(outgoing[0] == c);
  VERIFY(outgoing[1] == d);
  VERIFY(outgoing[2] == e);
  VERIFY(outgoing[3] == f);
  // Since the string is baryonic,
  // the most forward proton has to carry the diquark,
  // which leads to a scaling factor of 0.7*2/3 and the most backward pion (f)
  // gets the other quark and a scaling factor of 0.7*1/2
  COMPARE(outgoing[0].initial_xsec_scaling_factor(), 0.);
  COMPARE(outgoing[1].initial_xsec_scaling_factor(),
          coherence_factor * 2. / 3.);
  COMPARE(outgoing[2].initial_xsec_scaling_factor(), 0.);
  COMPARE(outgoing[3].initial_xsec_scaling_factor(), coherence_factor / 2.0);

  incoming = {e, f};  // Mesonic string
  e.set_4momentum(0.138, {0., 0., 1.0});
  f.set_4momentum(0.138, {0., 0., 0.5});
  c.set_4momentum(0.938, {0., 0., -0.5});
  d.set_4momentum(0.938, {0., 0., -1.0});
  outgoing = {f, c, d, e};  // again in random order
  // Since it is a Mesonic string, the valence quarks to distribute are
  // a quark and an anti-quark. Particle d will carry the quark and is assigned
  // a scaling factor of 0.7 * 1/3. On the other side of the string is a meson
  // (Particle e). This contains an anti-quark and will therefore get a scaling
  // factor of 0.7 * 1/2.
  baryon_string = 0;
  StringProcess::assign_all_scaling_factors(baryon_string, outgoing, evec_coll,
                                            coherence_factor);
  COMPARE(outgoing[0].initial_xsec_scaling_factor(), 0.5 * coherence_factor);
  COMPARE(outgoing[1].initial_xsec_scaling_factor(), 0);
  COMPARE(outgoing[2].initial_xsec_scaling_factor(), 0);
  COMPARE(outgoing[3].initial_xsec_scaling_factor(), coherence_factor / 3.);
  VERIFY(outgoing[3] == d);
  // While partile d was now the last particle in the list, if we exchange the
  // momenta of d and c, particle c will be assigned the scaling factor.
  // Even though particle c is an anti-baryon, this is correct, since the meson
  // on the other end of the string can also carry the quark instead.
  c.set_4momentum(0.938, {0., 0., -1.0});
  d.set_4momentum(0.938, {0., 0., -0.5});
  outgoing = {c, d, e, f};
  StringProcess::assign_all_scaling_factors(baryon_string, outgoing, evec_coll,
                                            coherence_factor);
  COMPARE(outgoing[0].initial_xsec_scaling_factor(), 0.5 * coherence_factor);
  COMPARE(outgoing[1].initial_xsec_scaling_factor(), 0.);
  COMPARE(outgoing[2].initial_xsec_scaling_factor(), 0.);
  COMPARE(outgoing[3].initial_xsec_scaling_factor(), coherence_factor / 3.);
  VERIFY(outgoing[3] == c);
}
