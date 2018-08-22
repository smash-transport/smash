/*
 *
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "setup.h"

#include "../include/smash/angles.h"
#include "../include/smash/scatteraction.h"
#include "Pythia8/Pythia.h"

using namespace smash;
using smash::Test::Momentum;
using smash::Test::Position;

TEST(init_particle_types) {
  Test::create_actual_particletypes();
  Test::create_actual_decaymodes();
}

constexpr double r_x = 0.1;
const FourVector pos_a = Position{0., -r_x, 0., 0.};
const FourVector pos_b = Position{0., r_x, 0., 0.};
const FourVector middle = (pos_a + pos_b) / 2.;

TEST(sorting) {
  ParticleData a{ParticleType::find(0x111)};  // pi0
  a.set_4position(pos_a);
  a.set_4momentum(Momentum{1.1, 1.0, 0., 0.});

  ParticleData b{ParticleType::find(0x111)};  // pi0
  b.set_4position(pos_b);
  a.set_4momentum(Momentum{1.1, 1.0, 0., 0.});

  constexpr double time1 = 1.;
  ScatterAction act1(a, b, time1);
  COMPARE(act1.get_interaction_point(), middle);

  constexpr double time2 = 1.1;
  ScatterAction act2(a, b, time2);
  VERIFY(act1 < act2);
}

TEST(elastic_collision) {
  // put particles in list
  Particles particles;
  ParticleData a{ParticleType::find(0x211)};  // pi+
  a.set_4position(pos_a);
  a.set_4momentum(Momentum{1.1, 1.0, 0., 0.});
  a.set_history(3, 1, ProcessType::None, 1.2, ParticleList{});

  ParticleData b{ParticleType::find(0x211)};  // pi+
  b.set_4position(pos_b);
  b.set_4momentum(Momentum{1.1, 1.0, 0., 0.});
  b.set_history(3, 1, ProcessType::None, 1.2, ParticleList{});

  a = particles.insert(a);
  b = particles.insert(b);

  // create action
  constexpr double time = 1.;
  ScatterAction act(a, b, time);
  ScatterAction act_copy(a, b, time);
  VERIFY(act.is_valid(particles));
  VERIFY(act_copy.is_valid(particles));

  // add elastic channel
  constexpr double sigma = 10.0;
  constexpr bool strings_switch = false;
  constexpr NNbarTreatment nnbar_treatment = NNbarTreatment::NoAnnihilation;
  act.add_all_scatterings(sigma, true, Test::all_reactions_included(), 0.,
                          strings_switch, false, false, nnbar_treatment);

  // check cross section
  COMPARE(act.cross_section(), sigma);

  // generate final state
  act.generate_final_state();

  // verify that the action is indeed elastic
  COMPARE(act.get_type(), ProcessType::Elastic);

  // verify that particles didn't change in the collision
  ParticleList in = act.incoming_particles();
  const ParticleList& out = act.outgoing_particles();
  VERIFY((in[0] == out[0] && in[1] == out[1]) ||
         (in[0] == out[1] && in[1] == out[0]));

  // verify that the particles keep their positions after elastic scattering
  COMPARE(out[0].position(), pos_a);
  COMPARE(out[1].position(), pos_b);

  // perform the action
  COMPARE(particles.front().id_process(), 1u);
  uint32_t id_process = 2;
  act.perform(&particles, id_process);
  id_process++;
  // check id_process
  COMPARE(particles.front().id_process(), 2u);
  COMPARE(particles.back().id_process(), 2u);
  COMPARE(id_process, 3u);

  // action should not be valid anymore
  VERIFY(!act.is_valid(particles));
  VERIFY(!act_copy.is_valid(particles));

  // verify that the particles don't change in the particle list
  VERIFY((in[0] == particles.front() && in[1] == particles.back()) ||
         (in[0] == particles.back() && in[1] == particles.front()));
}

TEST(outgoing_valid) {
  // create a proton and a pion
  ParticleData p1{ParticleType::find(0x2212)};
  ParticleData p2{ParticleType::find(0x111)};
  // set position
  p1.set_4position(pos_a);
  p2.set_4position(pos_b);
  // set momenta
  constexpr double p_x = 0.1;
  p1.set_4momentum(p1.pole_mass(), p_x, 0., 0.);
  p2.set_4momentum(p2.pole_mass(), -p_x, 0., 0.);

  // put in particles object
  Particles particles;
  particles.insert(p1);
  particles.insert(p2);

  // get valid copies back
  ParticleList plist = particles.copy_to_vector();
  auto p1_copy = plist[0];
  auto p2_copy = plist[1];
  VERIFY(particles.is_valid(p1_copy) && particles.is_valid(p2_copy));

  // construct action
  ScatterActionPtr act;
  act = make_unique<ScatterAction>(p1_copy, p2_copy, 0.2);
  VERIFY(act != nullptr);
  COMPARE(p2_copy.type(), ParticleType::find(0x111));

  // add processes
  constexpr double elastic_parameter = 0.;  // don't include elastic scattering
  constexpr bool strings_switch = false;
  constexpr NNbarTreatment nnbar_treatment = NNbarTreatment::NoAnnihilation;
  act->add_all_scatterings(elastic_parameter, true,
                           Test::all_reactions_included(), 0., strings_switch,
                           false, false, nnbar_treatment);

  VERIFY(act->cross_section() > 0.);

  // perform actions
  VERIFY(act->is_valid(particles));
  act->generate_final_state();
  VERIFY(act->get_type() != ProcessType::Elastic);
  const uint32_t id_process = 1;
  act->perform(&particles, id_process);
  COMPARE(id_process, 1u);

  // check the outgoing particles
  const ParticleList& outgoing_particles = act->outgoing_particles();
  VERIFY(outgoing_particles.size() > 0u);  // should be at least one
  VERIFY(particles.is_valid(outgoing_particles[0]));
  VERIFY(outgoing_particles[0].id() > p1_copy.id());
  VERIFY(outgoing_particles[0].id() > p2_copy.id());
  // verify that particle is placed in the middle between the incoming ones
  COMPARE(outgoing_particles[0].position(), middle);
}

TEST(pythia_running) {
  // create two protons
  ParticleData p1{ParticleType::find(0x2212)};
  ParticleData p2{ParticleType::find(0x2212)};
  // set position
  p1.set_4position(pos_a);
  p2.set_4position(pos_b);
  // set momenta
  constexpr double p_x = 3.0;
  p1.set_4momentum(p1.pole_mass(), p_x, 0., 0.);
  p2.set_4momentum(p2.pole_mass(), -p_x, 0., 0.);

  // put in particles object
  Particles particles;
  particles.insert(p1);
  particles.insert(p2);

  // get valid copies back
  ParticleList plist = particles.copy_to_vector();
  auto p1_copy = plist[0];
  auto p2_copy = plist[1];
  VERIFY(particles.is_valid(p1_copy) && particles.is_valid(p2_copy));

  // construct action
  ScatterActionPtr act;
  act = make_unique<ScatterAction>(p1_copy, p2_copy, 0.2, false, 1.0);
  std::unique_ptr<StringProcess> string_process_interface =
      make_unique<StringProcess>(1.0, 1.0, 0.5, 0.001, 1.0, 2.5, 0.217, 0.081,
                                 0.7, 0.68, 0.98, 0.25, 1.0, true);
  act->set_string_interface(string_process_interface.get());
  VERIFY(act != nullptr);
  COMPARE(p2_copy.type(), ParticleType::find(0x2212));

  // add processes
  constexpr double elastic_parameter = 0.;  // don't include elastic scattering
  constexpr bool strings_switch = true;
  constexpr NNbarTreatment nnbar_treatment = NNbarTreatment::NoAnnihilation;
  act->add_all_scatterings(elastic_parameter, false, ReactionsBitSet(), 0.,
                           strings_switch, false, false, nnbar_treatment);

  VERIFY(act->cross_section() > 0.);

  // perform actions
  VERIFY(act->is_valid(particles));
  act->generate_final_state();
  VERIFY(act->get_type() != ProcessType::Elastic);
  COMPARE(is_string_soft_process(act->get_type()), true) << act->get_type();
  const uint32_t id_process = 1;
  act->perform(&particles, id_process);
  COMPARE(id_process, 1u);

  // check the outgoing particles
  const ParticleList& outgoing_particles = act->outgoing_particles();
  VERIFY(outgoing_particles.size() > 0u);  // should be at least one
  VERIFY(particles.is_valid(outgoing_particles[0]));
  VERIFY(outgoing_particles[0].id() > p1_copy.id());
  VERIFY(outgoing_particles[0].id() > p2_copy.id());
}

TEST(no_strings) {
  // create two protons
  // TODO(steinberg): test more pairs after Jan's restructuring is merged
  ParticleData p1{ParticleType::find(0x2212)};
  ParticleData p2{ParticleType::find(0x2212)};
  // set position
  p1.set_4position(pos_a);
  p2.set_4position(pos_b);
  // set momenta
  constexpr double p_x = 1.0;
  p1.set_4momentum(p1.pole_mass(), p_x, 0., 0.);
  p2.set_4momentum(p2.pole_mass(), -p_x, 0., 0.);

  // put in particles object
  Particles particles;
  particles.insert(p1);
  particles.insert(p2);

  // get valid copies back
  ParticleList plist = particles.copy_to_vector();
  auto p1_copy = plist[0];
  auto p2_copy = plist[1];
  VERIFY(particles.is_valid(p1_copy) && particles.is_valid(p2_copy));

  // construct action
  ScatterActionPtr act;
  act = make_unique<ScatterAction>(p1_copy, p2_copy, 0.2, false, 1.0);
  VERIFY(act != nullptr);
  COMPARE(p2_copy.type(), ParticleType::find(0x2212));

  // add processes
  constexpr double elastic_parameter = 0.;  // don't include elastic scattering
  constexpr bool strings_switch = false;
  constexpr NNbarTreatment nnbar_treatment = NNbarTreatment::NoAnnihilation;
  act->add_all_scatterings(elastic_parameter, false,
                           Test::all_reactions_included(), 0., strings_switch,
                           false, false, nnbar_treatment);

  VERIFY(act->cross_section() > 0.);

  // perform actions
  VERIFY(act->is_valid(particles));
  act->generate_final_state();
  const uint32_t id_process = 1;
  act->perform(&particles, id_process);
  COMPARE(id_process, 1u);

  // check the outgoing particles
  const ParticleList& outgoing_particles = act->outgoing_particles();
  VERIFY(outgoing_particles.size() > 0u);  // should be at least one
  VERIFY(particles.is_valid(outgoing_particles[0]));
}

TEST(update_incoming) {
  // put particles in list
  Particles particles;
  ParticleData a{ParticleType::find(0x211)};  // pi+
  a.set_4position(pos_a);
  a.set_4momentum(Momentum{1.1, 1.0, 0., 0.});

  ParticleData b{ParticleType::find(0x211)};  // pi+
  b.set_4position(pos_b);
  b.set_4momentum(Momentum{1.1, -1.0, 0., 0.});

  a = particles.insert(a);
  b = particles.insert(b);

  // create action
  constexpr double time = 0.2;
  ScatterAction act(a, b, time);
  VERIFY(act.is_valid(particles));

  // add elastic channel
  constexpr double sigma = 10.0;
  bool string_switch = true;
  NNbarTreatment nnbar_treatment = NNbarTreatment::NoAnnihilation;
  act.add_all_scatterings(sigma, true, Test::all_reactions_included(), 0.,
                          string_switch, false, false, nnbar_treatment);

  // change the position of one of the particles
  const FourVector new_position(0.1, 0., 0., 0.);
  particles.front().set_4position(new_position);

  // update the action
  act.update_incoming(particles);
  COMPARE(act.incoming_particles()[0].position(), new_position);
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
  // decompose pion+ into u, dbar
  PdgCode pdg_piplus = PdgCode(0x211);
  StringProcess::make_string_ends(pdg_piplus, id1, id2);
  VERIFY(id1 == 2 && id2 == -1);
  // decompose pion- into d, ubar
  PdgCode pdg_piminus = PdgCode(-0x211);
  StringProcess::make_string_ends(pdg_piminus, id1, id2);
  VERIFY(id1 == 1 && id2 == -2);
  // decompose proton into u, ud-diquark or d, uu-diquark
  PdgCode pdg_proton = PdgCode(0x2212);
  StringProcess::make_string_ends(pdg_proton, id1, id2);
  VERIFY((id1 == 1 && id2 == 2203) || (id1 == 2 && (id2 == 2101 || 2103)));
  // decompose anti-proton ubar, ud-antidiquark or dbar, uu-antidiquark
  PdgCode pdg_antip = PdgCode(-0x2212);
  StringProcess::make_string_ends(pdg_antip, id1, id2);
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
  c.set_4momentum(0.938, {0., 0., -1.});
  d.set_4momentum(0.938, {0., 0., -0.5});
  e.set_4momentum(0.138, {0., 0., 0.5});
  f.set_4momentum(0.138, {0., 0., 1.});
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
  // Since the string is baryonic, the proton has to carry the diquark,
  // which leads to a scaling factor of 0.7*2/3 and the faster pion (f)
  // gets the other quark and a scaling factor of 0.7*1/2
  COMPARE(outgoing[0].initial_xsec_scaling_factor(), 0.);
  COMPARE(outgoing[1].initial_xsec_scaling_factor(),
          coherence_factor * 2. / 3.);
  COMPARE(outgoing[2].initial_xsec_scaling_factor(), 0.);
  COMPARE(outgoing[3].initial_xsec_scaling_factor(), coherence_factor / 2.0);

  incoming = {e, f};  // Mesonic string
  e.set_4momentum(0.138, {0., 0., -1.0});
  f.set_4momentum(0.138, {0., 0., -0.5});
  c.set_4momentum(0.938, {0., 0., 0.5});
  d.set_4momentum(0.938, {0., 0., 1.0});
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
  c.set_4momentum(0.938, {0., 0., 1.0});
  d.set_4momentum(0.938, {0., 0., 0.5});
  outgoing = {c, d, e, f};
  StringProcess::assign_all_scaling_factors(baryon_string, outgoing, evec_coll,
                                            coherence_factor);
  COMPARE(outgoing[0].initial_xsec_scaling_factor(), 0.5 * coherence_factor);
  COMPARE(outgoing[1].initial_xsec_scaling_factor(), 0.);
  COMPARE(outgoing[2].initial_xsec_scaling_factor(), 0.);
  COMPARE(outgoing[3].initial_xsec_scaling_factor(), coherence_factor / 3.);
  VERIFY(outgoing[3] == c);
}

TEST(pdg_map_for_pythia) {
  int pdgid_mapped = 0;

  // pi+ is mapped onto pi+
  PdgCode pdg_piplus = PdgCode(0x211);
  pdgid_mapped = StringProcess::pdg_map_for_pythia(pdg_piplus);
  VERIFY(pdgid_mapped == 211);

  // pi0 is mapped onto pi0
  PdgCode pdg_pi0 = PdgCode(0x111);
  pdgid_mapped = StringProcess::pdg_map_for_pythia(pdg_pi0);
  VERIFY(pdgid_mapped == 111);

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
