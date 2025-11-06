/*
 *
 *    Copyright (c) 2020,2022,2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include <iostream>

#include "Pythia8/Pythia.h"

#include "histogram.h"
#include "setup.h"
#include "smash/angles.h"
#include "smash/random.h"
#include "smash/scatteraction.h"

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

static std::unique_ptr<StringProcess> dummy_string_process() {
  auto sp = std::make_unique<StringProcess>(1., 1., .0, .001, .0, .0, 1., 1.,
                                            .0, .0, .5, .0, .21, .0, .0, true,
                                            1. / 3., true, 0., false);

  return sp;
}

TEST(common_setup) {
  // StringProcess to use member functions
  std::unique_ptr<StringProcess> sp = dummy_string_process();
  // Pythia object to work with
  Pythia8::Pythia pythia_interface{PYTHIA_XML_DIR, false};
  sp->common_setup_pythia(&pythia_interface, 1.0, 1.0, .3, .5, .7, .9);

  // Verify that all parameters were set accordingly
  VERIFY(pythia_interface.settings.mode("ParticleData:modeBreitWigner") == 4);
  FUZZY_COMPARE(pythia_interface.settings.parm("MultipartonInteractions:pTmin"),
                1.5);
  VERIFY(pythia_interface.settings.mode("MultipartonInteractions:nSample") ==
         100000);
  FUZZY_COMPARE(pythia_interface.settings.parm("StringPT:sigma"), .9);
  FUZZY_COMPARE(pythia_interface.settings.parm("StringFlav:probQQtoQ"), 1.0);
  FUZZY_COMPARE(pythia_interface.settings.parm("StringFlav:probStoUD"), 1.0);
  FUZZY_COMPARE(pythia_interface.settings.parm("StringFlav:popcornRate"), .3);
  FUZZY_COMPARE(pythia_interface.settings.parm("StringZ:aLund"), .5);
  FUZZY_COMPARE(pythia_interface.settings.parm("StringZ:bLund"), .7);
  VERIFY(pythia_interface.settings.word("PDF:pSet") == "13");
  VERIFY(pythia_interface.settings.word("PDF:pSetB") == "13");
  VERIFY(pythia_interface.settings.word("PDF:piSet") == "1");
  VERIFY(pythia_interface.settings.word("PDF:piSetB") == "1");
  VERIFY(pythia_interface.settings.mode("Beams:idA") == 2212);
  VERIFY(pythia_interface.settings.mode("Beams:idB") == 2212);
  FUZZY_COMPARE(pythia_interface.settings.parm("Beams:eCM"), 10.);
  VERIFY(pythia_interface.settings.flag("Random:setSeed") == 1);
  VERIFY(pythia_interface.settings.flag("Print:quiet") == 1);
  VERIFY(pythia_interface.settings.flag("HadronLevel:Decay") == 0);
  VERIFY(pythia_interface.settings.parm("Check:epTolErr") == 1e-6);
  VERIFY(pythia_interface.settings.parm("Check:epTolWarn") == 1e-8);
}

TEST(append_final) {
  // Create StringProcess to work with
  std::unique_ptr<StringProcess> sp = dummy_string_process();

  // ParticleData object to calculate final state for
  ParticleData a{ParticleType::find(0x211)};
  ParticleData b{ParticleType::find(0x111)};
  // Easy momentum for easier test values
  a.set_4momentum(0.138, 0., 0., 0.138);
  b.set_4momentum(0.138, 0., 0., -0.138);
  ParticleList intermediate = {a, b};

  // Vectors for use in tested function
  // Values make for easily calucated test values
  FourVector uString = {1., .0, .0, .0};
  ThreeVector evecLong = {.0, .0, .0};

  // Call tested function
  sp->append_final_state(intermediate, uString, evecLong);

  // Formation time is 0 due to the soft_t_form_ in the StringProcess
  // vx and vy remain 0 even with boosting
  // As vz starts at 1 it simply gets boosted to inverse_gamma
  COMPARE_ABSOLUTE_ERROR(.0, sp->get_final_state()[0].formation_time(), 1e-7);
  COMPARE_ABSOLUTE_ERROR(.0, sp->get_final_state()[1].formation_time(), 1e-7);
  COMPARE_ABSOLUTE_ERROR(.0, sp->get_final_state()[0].velocity().x1(), 1e-7);
  COMPARE_ABSOLUTE_ERROR(.0, sp->get_final_state()[1].velocity().x1(), 1e-7);
  COMPARE_ABSOLUTE_ERROR(.0, sp->get_final_state()[0].velocity().x2(), 1e-7);
  COMPARE_ABSOLUTE_ERROR(.0, sp->get_final_state()[1].velocity().x2(), 1e-7);
  COMPARE_ABSOLUTE_ERROR(.7071067812, sp->get_final_state()[0].velocity().x3(),
                         1e-7);
  COMPARE_ABSOLUTE_ERROR(-0.7071067812,
                         sp->get_final_state()[1].velocity().x3(), 1e-7);
}

TEST(initialization) {
  std::unique_ptr<StringProcess> sp = dummy_string_process();
  ParticleData a{ParticleType::find(0x2212)};
  a.set_4momentum(1., 0., 0., 1.);

  ParticleData b{ParticleType::find(0x2212)};
  b.set_4momentum(0.9380, 0., 0., -1.);

  sp->init({a, b}, 0.);

  // Check if all values are as expected
  VERIFY(sp->get_PDGs()[0] == pdg::p);
  VERIFY(sp->get_PDGs()[1] == pdg::p);
  COMPARE_ABSOLUTE_ERROR(sp->get_massA(), 1.00, 1e-3);
  COMPARE_ABSOLUTE_ERROR(sp->get_massB(), 0.938, 1e-3);
  COMPARE_ABSOLUTE_ERROR(sp->get_sqrts(), 2.78529, 1e-5);
  COMPARE_ABSOLUTE_ERROR(sp->get_tcoll(), 0.00000, 1e-5);

  // Test vectors to compare with
  FourVector plab0_test{1.41421, 0.00000, 0.00000, 1.00000};
  FourVector plab1_test{1.37107, 0.00000, 0.00000, -1.00000};
  FourVector pcom0_test{1.41421, 0.00000, 0.00000, 1.00000};
  FourVector pcom1_test{1.37107, 0.00000, 0.00000, -1.00000};
  FourVector ucom_test{1.00000, 0.00000, 0.00000, 0.00000};
  ThreeVector vcom_test{0.00000, 0.00000, 0.00000};

  // Vector comparison
  for (int i = 0; i < 4; ++i) {
    COMPARE_ABSOLUTE_ERROR(sp->get_plab()[0][i], plab0_test[i], 1e-5);
    COMPARE_ABSOLUTE_ERROR(sp->get_plab()[1][i], plab1_test[i], 1e-5);
    COMPARE_ABSOLUTE_ERROR(sp->get_pcom()[0][i], pcom0_test[i], 1e-5);
    COMPARE_ABSOLUTE_ERROR(sp->get_pcom()[1][i], pcom1_test[i], 1e-5);
    COMPARE_ABSOLUTE_ERROR(sp->get_ucom()[i], ucom_test[i], 1e-5);
    // ThreeVector comparison excludes i = 3
    if (i < 3)
      COMPARE_ABSOLUTE_ERROR(sp->get_vcom()[i], vcom_test[i], 1e-5);
  }
}

TEST(rearrange_ex) {
  // StringProcess to use member functions
  std::unique_ptr<StringProcess> sp = dummy_string_process();
  // Array for total quark numbers
  std::array<int, 5> tot_quark = {0, 0, 0, 0, 0};
  // Arrays with the excess constituents
  std::array<std::array<int, 5>, 2> exc_quark = {
      {{1, 0, 1, 0, 0}, {0, 1, 1, 0, 0}}};
  std::array<std::array<int, 5>, 2> exc_antiq = {
      {{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}}};
  // Test arrays, nothing has changed as all nquark_final are positive
  std::array<std::array<int, 5>, 2> quark_test = {
      {{1, 0, 1, 0, 0}, {0, 1, 1, 0, 0}}};
  std::array<std::array<int, 5>, 2> antiq_test = {
      {{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}}};

  // Rearrange the excess constituents
  sp->rearrange_excess(tot_quark, exc_quark, exc_antiq);
  // Verify that nothing has changed
  VERIFY(exc_quark == quark_test);
  VERIFY(exc_antiq == antiq_test);

  tot_quark = {-1, 0, 0, 0, 0};
  exc_quark = {{{-1, 0, 0, 0, 0}, {1, 0, 0, 0, 0}}};
  exc_antiq = {{{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}}};
  quark_test = {{{0, 0, 0, 0, 0}, {1, 0, 0, 0, 0}}};
  antiq_test = {{{1, 0, 0, 0, 0}, {0, 0, 0, 0, 0}}};
  sp->rearrange_excess(tot_quark, exc_quark, exc_antiq);
  VERIFY(exc_quark == quark_test);
  VERIFY(exc_antiq == antiq_test);
}

TEST(find_excess) {
  // StringProcess to use member functions
  std::unique_ptr<StringProcess> sp = dummy_string_process();
  // PDG codes for proton and neutron
  PdgCode actual = pdg::p;
  PdgCode mapped = pdg::n;
  // Arrays for excess numbers
  std::array<int, 5> q_excess = {0, 0, 0, 0, 0};
  std::array<int, 5> antiq_excess = {0, 0, 0, 0, 0};
  // Test arrays with expected outcome
  std::array<int, 5> q_test = {-1, 1, 0, 0, 0};
  std::array<int, 5> antiq_test = {0, 0, 0, 0, 0};
  // Find excess quarks
  sp->find_excess_constituent(actual, mapped, q_excess, antiq_excess);
  // Verify predicted values
  VERIFY(q_excess == q_test);
  VERIFY(antiq_excess == antiq_test);

  // Same as the above using different particles
  actual = pdg::pi_m;
  mapped = pdg::n;
  q_excess = {0, 0, 0, 0, 0};
  antiq_excess = {0, 0, 0, 0, 0};
  q_test = {-1, -1, 0, 0, 0};
  antiq_test = {0, 1, 0, 0, 0};
  sp->find_excess_constituent(actual, mapped, q_excess, antiq_excess);
  VERIFY(q_excess == q_test);
  VERIFY(antiq_excess == antiq_test);
}

TEST(restore_constituents) {
  // Create StringProcess to work with member functions
  std::unique_ptr<StringProcess> sp = dummy_string_process();
  // create Pythia class object to simulate p-p collision
  Pythia8::Pythia pythia(PYTHIA_XML_DIR, false);
  pythia.readString("Print:quiet = on");
  pythia.readString("Top:gg2ttbar = on");
  pythia.readString("Beams:eCM = 8000.");
  pythia.init();

  // Arrays which contain the excess quark flavours
  // First array ensures conversion of a top quark to a down quark
  std::array<std::array<int, 5>, 2> exc_quark = {
      {{1, 0, 0, 0, 0}, {-1, 0, 0, 0, 0}}};
  std::array<std::array<int, 5>, 2> exc_antiq = {
      {{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}}};

  // Generate first event
  pythia.next();

  // Apply restore_constituent function
  sp->restore_constituent(pythia.process, exc_quark, exc_antiq);

  // Add total system momentum to substract from
  double px = pythia.process[0].px();
  double py = pythia.process[0].py();

  bool top_exists = false;
  bool down_exists = false;

  // Loop over particles and substract momenta
  // Also check that there is no top quark and a new down quark
  for (int i = 1; i < pythia.process.size(); ++i) {
    px -= pythia.process[i].px();
    py -= pythia.process[i].py();
    if (pythia.process[i].id() == 6) {
      top_exists = true;
    } else if (pythia.process[i].id() == 1) {
      down_exists = true;
    }
  }

  // Check that the sum of momenta in x and y direction remains 0
  COMPARE_ABSOLUTE_ERROR(px, .0, 1e-7);
  COMPARE_ABSOLUTE_ERROR(py, .0, 1e-7);

  // Verify that the energy is conserved at 8000
  COMPARE_ABSOLUTE_ERROR(pythia.process[0].e(), 8000., 1e-3);

  // Verify that there is no top quark and a new down quark
  VERIFY(top_exists == false);
  VERIFY(down_exists == true);
}

TEST(replace_const) {
  // Create StringProcess to work with member functions
  std::unique_ptr<StringProcess> sp = dummy_string_process();
  // Create particle entry for an electron
  Pythia8::ParticleDataEntry entry1(11, "e", "e+");
  // Create Pythia Particle, an electron in this case
  Pythia8::Particle part1(11);
  // Set Data Entry Pointer to the corresponding entry
  part1.setPDEPtr(std::make_shared<Pythia8::ParticleDataEntry>(entry1));

  // Initialize array with excess constituents
  // Initialize array to check correct function behavior
  // Expected output is the same as the input as an electron
  // is neither quark nor diquark
  std::array<int, 5> exc_const = {0, 1, 0, 1, 0};
  std::array<int, 5> test_const = {0, 1, 0, 1, 0};

  // Call function to be tested
  sp->replace_constituent(part1, exc_const);

  // Check content of excess consituent array
  VERIFY(exc_const == test_const);
  VERIFY(part1.id() == 11);

  // Case 2: Particle is a quark but excess array is zero
  Pythia8::ParticleDataEntry entry2(1, "d", "dbar");
  Pythia8::Particle part2(1);
  exc_const = {0, 0, 0, 0, 0};
  test_const = {0, 0, 0, 0, 0};
  sp->replace_constituent(part2, exc_const);
  VERIFY(exc_const == test_const);
  VERIFY(part2.id() == 1);

  // Case 3: Particle is a quark, excess reaches 0
  Pythia8::ParticleDataEntry entry3(1, "d", "dbar");
  Pythia8::Particle part3(1);
  part3.setPDEPtr(std::make_shared<Pythia8::ParticleDataEntry>(entry3));
  exc_const = {-1, 0, 0, 1, 0};
  // Outcome is expected to be 0 as the charm converts to down
  test_const = {0, 0, 0, 0, 0};
  sp->replace_constituent(part3, exc_const);
  VERIFY(exc_const == test_const);
  VERIFY(part3.id() == 4);

  // Case 4: Particle is a diquark, excess reaches 0
  Pythia8::ParticleDataEntry entry4(2101, "ud", "udbar");
  Pythia8::Particle part4(2101);
  part4.setPDEPtr(std::make_shared<Pythia8::ParticleDataEntry>(entry4));
  exc_const = {-1, -1, 0, 1, 1};
  // Outcome is expected to be 0 again
  test_const = {0, 0, 0, 0, 0};
  sp->replace_constituent(part4, exc_const);
  VERIFY(exc_const == test_const);
  VERIFY(part4.id() == 5401);
}

TEST(find_total_number_constituent) {
  // Create Pythia event to fill with an event later
  Pythia8::Event intermediate;

  // Create Pythia object to work with
  std::unique_ptr<Pythia8::Pythia> pythia_hadron =
      std::make_unique<Pythia8::Pythia>(PYTHIA_XML_DIR, false);
  pythia_hadron->readString("Print:quiet = on");
  pythia_hadron->readString("ProcessLevel:all = off");
  pythia_hadron->readString("Top:gg2ttbar = on");
  pythia_hadron->readString("Beams:eCM = 8000.");
  pythia_hadron->init();

  // Initialize the intermediate event trivially
  intermediate.init("intermediate partons", &pythia_hadron->particleData);

  // String process to be able to call the member functions
  std::unique_ptr<StringProcess> sp = dummy_string_process();

  // Arrays for quark and antiquark content
  std::array<int, 5> nquark;
  std::array<int, 5> nantiq;

  // Call tested function
  sp->find_total_number_constituent(intermediate, nquark, nantiq);

  // All should contain 0 due to the default initialization
  for (int i = 0; i < 5; ++i) {
    VERIFY((nquark[i] == 0) && (nantiq[i] == 0));
  }
}

/**
 * Compare sampled values of the LUND function to analytical ones via:
 * \f$ f(z) = \frac{1}{z} (1 - z)^a \exp{ \left(- \frac{b m_T^2}{z} \right) }
 * Using the same framework as the tests in random.cc do.
 */
TEST(string_zlund) {
  test_distribution(
      1e7, 0.0001, []() { return StringProcess::sample_zLund(1, 1, 1); },
      [](double x) { return 1 / x * (1. - x) * exp(-1. / x); });
}

TEST(string_incoming_lightcone_momenta) {
  std::unique_ptr<StringProcess> sp = dummy_string_process();

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
  std::unique_ptr<StringProcess> sp = dummy_string_process();

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
