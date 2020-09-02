/*
 *    Copyright (c) 2016-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include <vir/test.h>  // This include has to be first

#include "setup.h"

#include "../include/smash/bremsstrahlungaction.h"
#include "../include/smash/scatteractionphoton.h"
#include "../include/smash/crosssectionsphoton.h"

using namespace smash;
using smash::Test::Momentum;

TEST(init_particle_types) {
  // enable debugging output
  create_all_loggers(Configuration(""));
  Test::create_actual_particletypes();
}

TEST(init_decay_modes) { Test::create_actual_decaymodes(); }

////
// Test photon production in binary scatterings
////

TEST(pi_rho0_pi_gamma) {
  // set up a π+ and a ρ0 in center-of-momentum-frame
  const ParticleType &type_pi = ParticleType::find(0x211);
  ParticleData pi{type_pi};
  pi.set_4momentum(type_pi.mass(),  // pole mass
                   ThreeVector(0., 0., 2.));
  const ParticleType &type_rho0 = ParticleType::find(0x113);
  ParticleData rho0{type_rho0};
  rho0.set_4momentum(type_rho0.mass(),  // pole mass
                     ThreeVector(0., 0., -2.));
  const int number_of_photons = 10000;
  ParticleList in{pi, rho0};
  const auto act =
      make_unique<ScatterActionPhoton>(in, 0.05, number_of_photons, 5.0);
  act->add_single_process();
  double tot_weight = 0.0;
  for (int i = 0; i < number_of_photons; i++) {
    act->generate_final_state();
    tot_weight += act->get_total_weight();
  }
  COMPARE_RELATIVE_ERROR(tot_weight, 0.000722419008, 0.08);
}

TEST(photon_and_hadron_reaction_type_function) {
  /*
   *creates possible photon reactions and also some that
   *produce no photons. Checks if the computed photon reaction type
   *is correct. Continues to check if the respective hadron types
   *are correct.
   */

  // setup for the particles
  const ParticleData pip{ParticleType::find(0x211)};
  const ParticleData pim{ParticleType::find(-0x211)};
  const ParticleData piz{ParticleType::find(0x111)};
  const ParticleData rhop{ParticleType::find(0x213)};
  const ParticleData rhom{ParticleType::find(-0x213)};
  const ParticleData rhoz{ParticleType::find(0x113)};
  const ParticleData eta{ParticleType::find(0x221)};
  const ParticleData p{ParticleType::find(0x2112)};

  // gets all possible outgoing hadrons
  static const ParticleTypePtr rho_z_particle_ptr =
      &ParticleType::find(pdg::rho_z);
  static const ParticleTypePtr rho_p_particle_ptr =
      &ParticleType::find(pdg::rho_p);
  static const ParticleTypePtr rho_m_particle_ptr =
      &ParticleType::find(pdg::rho_m);
  static const ParticleTypePtr pi_z_particle_ptr =
      &ParticleType::find(pdg::pi_z);
  static const ParticleTypePtr pi_m_particle_ptr =
      &ParticleType::find(pdg::pi_m);

  // puts the particles in lists
  const ParticleList l1{pip, pim}, l2{rhop, pim}, l3{p, pim}, l4{pip, eta},
      l5{pip, piz}, l6{piz, pip}, l7{pim, piz}, l8{piz, pim}, l9{pim, rhoz},
      l10{piz, rhom}, l11{piz, rhoz};

  // test if the photon reaction type is correct
  VERIFY(ScatterActionPhoton::photon_reaction_type(l1) !=
         ScatterActionPhoton::ReactionType::no_reaction);
  VERIFY(ScatterActionPhoton::photon_reaction_type(l2) !=
         ScatterActionPhoton::ReactionType::no_reaction);
  VERIFY(ScatterActionPhoton::photon_reaction_type(l3) ==
         ScatterActionPhoton::ReactionType::no_reaction);
  VERIFY(ScatterActionPhoton::photon_reaction_type(l4) ==
         ScatterActionPhoton::ReactionType::no_reaction);

  // test if the respective reactions from above result in
  // the right hadron type.

  // these also test if the produced hadron is invariant
  // if you swap the partners
  VERIFY(ScatterActionPhoton::outgoing_hadron_type(l5) == rho_p_particle_ptr);
  VERIFY(ScatterActionPhoton::outgoing_hadron_type(l6) == rho_p_particle_ptr);
  VERIFY(ScatterActionPhoton::outgoing_hadron_type(l7) == rho_m_particle_ptr);
  VERIFY(ScatterActionPhoton::outgoing_hadron_type(l8) == rho_m_particle_ptr);

  // these simply test more cases
  VERIFY(ScatterActionPhoton::outgoing_hadron_type(l9) == pi_m_particle_ptr);
  VERIFY(ScatterActionPhoton::outgoing_hadron_type(l10) == pi_m_particle_ptr);
  VERIFY(ScatterActionPhoton::outgoing_hadron_type(l1) == rho_z_particle_ptr);
  VERIFY(ScatterActionPhoton::outgoing_hadron_type(l11) == pi_z_particle_ptr);
}

TEST(check_kinematic_thresholds) {
  /*
   * Make sure pi + pi -> rho + photon process is only executed if sqrt(s) is
   * high enough to not only create final state rho, but also to assign momentum
   * to rho and photon.
   */

  Particles particles;
  ParticleData a{ParticleType::find(0x211)};   // pi+
  ParticleData b{ParticleType::find(-0x211)};  // pi-
  ParticleData c{ParticleType::find(-0x211)};  // pi-

  /*
   * Pick energies such that energy_a + energy_b > m_rho_min + really_small
   * and energy_a + energy_c < m_rho_min + really_small.
   * Hence a+b should be performed while a+c should be rejected.
   */

  double energy_a = sqrt(pion_mass * pion_mass + 0.001 * 0.001);
  double energy_b = sqrt(pion_mass * pion_mass + 0.5 * 0.5);
  double energy_c = sqrt(pion_mass * pion_mass + 0.002 * 0.002);

  a.set_4momentum(Momentum{energy_a, 0.001, 0., 0.});
  b.set_4momentum(Momentum{energy_b, -0.5, 0., 0.});
  c.set_4momentum(Momentum{energy_c, -0.002, 0., 0.});

  a = particles.insert(a);
  b = particles.insert(b);
  c = particles.insert(c);

  // create underlying hadronic interactions
  constexpr double time = 1.;
  ScatterAction act_highE(a, b, time);
  ScatterAction act_lowE(a, c, time);

  // create photon scatter action
  ParticleList in_highE = act_highE.incoming_particles();
  ParticleList in_lowE = act_lowE.incoming_particles();
  ScatterActionPhoton photonAct_highE(in_highE, time, 1, 30.0);
  ScatterActionPhoton photonAct_lowE(in_lowE, time, 1, 30.0);

  VERIFY(
      photonAct_highE.is_kinematically_possible(energy_a + energy_b, in_highE));
  VERIFY(
      !photonAct_lowE.is_kinematically_possible(energy_a + energy_c, in_lowE));
}

////
// Test photon production in Bremsstrahlung processes
////

TEST(gen_final_state) {
  // set up a π+ and a π- in center-of-momentum-frame
  const ParticleType &type_pip = ParticleType::find(0x211);
  ParticleData pip{type_pip};
  pip.set_4momentum(type_pip.mass(), ThreeVector(0., 0., 2.));
  const ParticleType &type_pim = ParticleType::find(-0x211);
  ParticleData pim{type_pim};
  pim.set_4momentum(type_pim.mass(), ThreeVector(0., 0., -2.));
  const ParticleType &type_photon = ParticleType::find(0x22);
  const int number_of_photons = 10;
  ParticleList in{pip, pim};

  // create bremsstrahlung action
  const auto act =
      make_unique<BremsstrahlungAction>(in, 0.05, number_of_photons, 20.0);
  act->add_single_process();

  // Sample photons, implicitly test sample_3body_phasespace() and
  // cross section functions
  double tot_weight = 0.0;
  for (int i = 0; i < number_of_photons; i++) {
    act->generate_final_state();
    tot_weight += act->get_total_weight();
    VERIFY(act->outgoing_particles().size() == 3);
    VERIFY(act->outgoing_particles()[0].type() == type_pip);
    VERIFY(act->outgoing_particles()[1].type() == type_pim);
    VERIFY(act->outgoing_particles()[2].type() == type_photon);
  }
  COMPARE_RELATIVE_ERROR(tot_weight, 1.84592, 1e-5);
}

TEST(bremsstrahlung_reaction_type_function) {
  const ParticleData pip{ParticleType::find(0x211)};
  const ParticleData pim{ParticleType::find(-0x211)};
  const ParticleData piz{ParticleType::find(0x111)};
  const ParticleData eta{ParticleType::find(0x221)};
  const ParticleData p{ParticleType::find(0x2112)};

  const ParticleList l1{pip, pim}, l2{piz, pim}, l3{pip, pip}, l4{piz, piz},
      l5{pim, pim}, l6{pip, piz}, l7{p, pim}, l8{pip, eta};

  VERIFY(BremsstrahlungAction::bremsstrahlung_reaction_type(l1) ==
         BremsstrahlungAction::ReactionType::pi_p_pi_m);
  VERIFY(BremsstrahlungAction::bremsstrahlung_reaction_type(l2) ==
         BremsstrahlungAction::ReactionType::pi_z_pi_m);
  VERIFY(BremsstrahlungAction::bremsstrahlung_reaction_type(l3) ==
         BremsstrahlungAction::ReactionType::pi_p_pi_p);
  VERIFY(BremsstrahlungAction::bremsstrahlung_reaction_type(l4) ==
         BremsstrahlungAction::ReactionType::pi_z_pi_z);
  VERIFY(BremsstrahlungAction::bremsstrahlung_reaction_type(l5) ==
         BremsstrahlungAction::ReactionType::pi_m_pi_m);
  VERIFY(BremsstrahlungAction::bremsstrahlung_reaction_type(l6) ==
         BremsstrahlungAction::ReactionType::pi_z_pi_p);
  VERIFY(BremsstrahlungAction::bremsstrahlung_reaction_type(l7) ==
         BremsstrahlungAction::ReactionType::no_reaction);
  VERIFY(BremsstrahlungAction::bremsstrahlung_reaction_type(l8) ==
         BremsstrahlungAction::ReactionType::no_reaction);
}

TEST(photon_cross_sections) {
  // calculate crosssections and compare to analytic values

  // (6a)
  double cross1 =
      CrosssectionsPhoton<ComputationMethod::Analytic>::xs_pi_pi_rho0(
          0.996, 0.776);
  std::cout << "Cross1:    " << cross1 << std::endl;
  std::cout << "Analytic1: " << 0.661515*0.3894 << std::endl;
  // COMPARE_ABSOLUTE_ERROR(cross1, .0, 1e-5);

  // (6b)
  double cross2 =
      CrosssectionsPhoton<ComputationMethod::Analytic>::xs_pi_pi0_rho(
          1.12, 0.9);
  std::cout << "Cross2:    " << cross2 << std::endl;
  std::cout << "Analytic2: " << 1.89737*0.3894 << std::endl;
  // COMPARE_ABSOLUTE_ERROR(cross2, .0, 1e-5);

  // (6c)
  double cross3 =
      CrosssectionsPhoton<ComputationMethod::Analytic>::xs_pi_rho0_pi(
          1.224, 0.776);
  std::cout << "Cross3:    " << cross3 << std::endl;
  std::cout << "Analytic3: " << 0.171629*0.3894 << std::endl;
  // COMPARE_ABSOLUTE_ERROR(cross3, .0, 1e-5);

  // (6f)
  double cross4 =
      CrosssectionsPhoton<ComputationMethod::Analytic>::xs_pi0_rho0_pi0(
          4.979, 0.776);
  std::cout << "Cross4:    " << cross4 << std::endl;
  std::cout << "Analytic4: " << 1.96799*0.3894 << std::endl;
  // COMPARE_ABSOLUTE_ERROR(cross4, .0, 1e-5);

  // (6d) + (6h)
  double cross5 =
      CrosssectionsPhoton<ComputationMethod::Analytic>::xs_pi0_rho_pi(
          1.103+5.058, 0.8+0.9);
  std::cout << "Cross5:    " << cross5 << std::endl;
  std::cout << "Analytic5: " << (0.0750229+1.63798)*0.3894 << std::endl;

  // COMPARE_ABSOLUTE_ERROR(cross5, .0, 1e-5);

  // (6e) + (6g)
  double cross6 =
      CrosssectionsPhoton<ComputationMethod::Analytic>::xs_pi_rho_pi0(
          2.948+0.973, 0.9+0.8);
  std::cout << "Cross6:    " << cross6 << std::endl;
  std::cout << "Analytic6: " << (0.195696+0.0552201)*0.3894 << std::endl;
  // COMPARE_ABSOLUTE_ERROR(cross6, .0, 1e-5);
}

TEST(diff_cross_section) {
  // calculate crosssections and compare with analytic values

  // (6a)
  double cross1 =
      CrosssectionsPhoton<ComputationMethod::Analytic>::xs_diff_pi_pi_rho0(
          1., -0.367612, 0.776);
  std::cout << "Cross1:    " << cross1 << std::endl;
  std::cout << "Analytic1: " << 8.79811*0.3894 << std::endl;
  // COMPARE_ABSOLUTE_ERROR(cross4, 0.033093, 1e-5);

  // (6b)
  double cross2 =
      CrosssectionsPhoton<ComputationMethod::Analytic>::xs_diff_pi_pi0_rho(
          1., -0.07066, 0.9);
  std::cout << "Cross2:    " << cross2 << std::endl;
  std::cout << "Analytic2: " << 4.85461*0.3894 << std::endl;
  // COMPARE_ABSOLUTE_ERROR(cross2, .0, 1e-5);

  // (6c)
  double cross3 =
      CrosssectionsPhoton<ComputationMethod::Analytic>::xs_diff_pi_rho0_pi(
          1., -0.316209, 0.776);
  std::cout << "Cross3:    " << cross3 << std::endl;
  std::cout << "Analytic3: " << 1.04421*0.3894 << std::endl;
  // COMPARE_ABSOLUTE_ERROR(cross3, .0, 1e-5);

  // (6f)
  double cross4 =
      CrosssectionsPhoton<ComputationMethod::Analytic>::xs_diff_pi0_rho0_pi0(
          1., -0.248786, 0.776);
  std::cout << "Cross4:    " << cross4 << std::endl;
  std::cout << "Analytic4: " << 0.13775*0.3894 << std::endl;
  // COMPARE_ABSOLUTE_ERROR(cross4, .0, 1e-5);


  // These two functions are defined in crosssectionsphoton.h
  // but were never implemented !!!

  // // (6d) + (6h)
  // double cross5 =
  //     CrosssectionsPhoton<ComputationMethod::Analytic>::xs_diff_pi0_rho_pi(
  //         1.14301, -0.423672, 0.739093);
  // std::cout << "Cross5:    " << cross5 << std::endl;
  // std::cout << "Analytic5: " << 0.070465*0.3894 << std::endl;
  // // COMPARE_ABSOLUTE_ERROR(cross5, .0, 1e-5);
  //
  // // (6e) + (6g)
  // double cross6 =
  //     CrosssectionsPhoton<ComputationMethod::Analytic>::xs_diff_pi0_rho_pi(
  //         1., -0.599925, 0.592348);
  // std::cout << "Cross6:    " << cross6 << std::endl;
  // std::cout << "Analytic6: " << 0.216896*0.3894 << std::endl;
  // // COMPARE_ABSOLUTE_ERROR(cross6, .0, 1e-5);

}
