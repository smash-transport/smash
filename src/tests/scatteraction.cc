/*
 *
 *    Copyright (c) 2015-2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/scatteraction.h"

#include <algorithm>

#include "Pythia8/Pythia.h"

#include "setup.h"
#include "smash/angles.h"
#include "smash/random.h"
#include "smash/scatteractionmulti.h"

using namespace smash;
using smash::Test::Momentum;
using smash::Test::Position;

TEST(init_particle_types) {
  Test::create_actual_particletypes();
  Test::create_actual_decaymodes();
  ParticleType::check_consistency();
  sha256::Hash hash;
  hash.fill(0);
  IsoParticleType::tabulate_integrals(hash, "");
}

constexpr double r_x = 0.1;
const FourVector pos_a = Position{0., -r_x, 0., 0.};
const FourVector pos_b = Position{0., r_x, 0., 0.};
const ThreeVector middle = ((pos_a + pos_b) / 2.).threevec();

TEST(sorting) {
  ParticleData a{ParticleType::find(0x111)};  // pi0
  a.set_4position(pos_a);
  a.set_4momentum(Momentum{1.1, 1.0, 0., 0.});

  ParticleData b{ParticleType::find(0x111)};  // pi0
  b.set_4position(pos_b);
  b.set_4momentum(Momentum{1.1, -1.0, 0., 0.});

  constexpr double time1 = 1.;
  ScatterAction act1(a, b, time1);
  COMPARE(act1.get_interaction_point(), FourVector(time1, middle));

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
  b.set_4momentum(Momentum{1.1, -1.0, 0., 0.});
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
  act.add_all_scatterings(Test::default_finder_parameters(
      sigma, NNbarTreatment::NoAnnihilation, Test::all_reactions_included(),
      false, false, false));

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
  act = std::make_unique<ScatterAction>(p1_copy, p2_copy, 0.2);
  VERIFY(act != nullptr);
  COMPARE(p2_copy.type(), ParticleType::find(0x111));

  // add processes

  double elastic_parameter = 0;  // don't include elastic scattering
  bool strings_switch = false;
  act->add_all_scatterings(Test::default_finder_parameters(
      elastic_parameter, NNbarTreatment::NoAnnihilation,
      Test::all_reactions_included(), strings_switch, false, false));

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
  // at the time of the collision
  ThreeVector interaction_point =
      middle + 0.2 * (p1.velocity() + p2.velocity()) / 2.0;
  COMPARE(outgoing_particles[0].position(), FourVector(0.2, interaction_point));
}

// TEST(spin_flip_in_elastic_scatterings) {
//   // put particles in list
//   Particles particles;
//   ParticleData a{ParticleType::find(smash::pdg::p)};  // Proton
//   a.set_4position(pos_a);
//   a.set_4momentum(Momentum{1.1, 1.0, 0., 0.});
//   a.set_spin_projection(1);
//   a.set_history(3, 1, ProcessType::None, 1.2, ParticleList{});

//   ParticleData b{ParticleType::find(0x1114)};  // Delta-
//   b.set_4position(pos_b);
//   b.set_4momentum(Momentum{1.1, -1.0, 0., 0.});
//   b.set_spin_projection(-3);
//   b.set_history(3, 1, ProcessType::None, 1.2, ParticleList{});

//   a = particles.insert(a);
//   b = particles.insert(b);

//   // create action and switch spin interactions on
//   constexpr double time = 1.;
//   ScatterAction act(a, b, time, false, 1.0, -1.0, false, true);

//   // add elastic channel
//   constexpr double sigma = 10.0;
//   act.add_all_scatterings(Test::default_finder_parameters(
//       sigma, NNbarTreatment::NoAnnihilation, Test::all_reactions_included(),
//       false, false, false, TotalCrossSectionStrategy::BottomUp, true));

//   // generate final state
//   act.generate_final_state();

//   // verify that the action is indeed elastic
//   COMPARE(act.get_type(), ProcessType::Elastic);

//   // verify that particle's spin projections were flipped
//   ParticleList in = act.incoming_particles();
//   ParticleList out = act.outgoing_particles();
//   VERIFY((in[0].spin_projection() == -out[0].spin_projection() &&
//           in[1].spin_projection() == -out[1].spin_projection()));
// }

TEST(cross_sections_symmetric) {
  // create a list of all particles
  const auto& all_types = ParticleType::list_all();
  int ntypes = all_types.size();
  int64_t seed = random::generate_63bit_seed();
  random::set_seed(seed);
  for (int i = 0; i < 42; i++) {
    // create a random pair of particles
    ParticleData p1{ParticleType::find(
        all_types[random::uniform_int(0, ntypes - 1)].pdgcode())};
    ParticleData p2{ParticleType::find(
        all_types[random::uniform_int(0, ntypes - 1)].pdgcode())};
    // set position
    constexpr double p_x = 3.0;
    p1.set_4position(pos_a);
    p2.set_4position(pos_b);
    // set momenta
    p1.set_4momentum(p1.pole_mass(), p_x, 0., 0.);
    p2.set_4momentum(p2.pole_mass(), -p_x, 0., 0.);

    // put in particles object
    Particles particles;
    particles.insert(p1);
    particles.insert(p2);

    // construct actions
    ScatterActionPtr act12, act21;
    act12 = std::make_unique<ScatterAction>(p1, p2, 0.2, false, 1.0);
    act21 = std::make_unique<ScatterAction>(p2, p1, 0.2, false, 1.0);
    auto string_process_interface = Test::default_string_process_interface();
    act12->set_string_interface(string_process_interface.get());
    act21->set_string_interface(string_process_interface.get());
    VERIFY(act12 != nullptr);
    VERIFY(act21 != nullptr);

    // add processes

    double elastic_parameter = -10.;  // no added elastic cross sections
    bool use_AQM = true;
    bool strings_with_probability = true;
    NNbarTreatment nnbar_treatment = NNbarTreatment::Strings;
    ScatterActionsFinderParameters finder_params =
        Test::default_finder_parameters(elastic_parameter, nnbar_treatment,
                                        Test::all_reactions_included(), true,
                                        use_AQM, strings_with_probability);

    act12->add_all_scatterings(finder_params);
    act21->add_all_scatterings(finder_params);

    VERIFY(act12->cross_section() >= 0.);
    VERIFY(act21->cross_section() >= 0.);

    // fuzzyness needs to be slightly larger than 1 for some compilers
    vir::test::setFuzzyness<double>(5);

    // check symmetry of the cross-section, i.e. xsec_AB == xsec_BA
    FUZZY_COMPARE(act12->cross_section(), act21->cross_section())
        << "Colliding " << p1.pdgcode() << " with " << p2.pdgcode()
        << " does not yield the same cross-section as " << p2.pdgcode()
        << " with " << p1.pdgcode() << "\nRandom seed used for test: " << seed;
  }
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
  act = std::make_unique<ScatterAction>(p1_copy, p2_copy, 0.2, false, 1.0);
  auto string_process_interface = Test::default_string_process_interface();
  act->set_string_interface(string_process_interface.get());
  VERIFY(act != nullptr);
  COMPARE(p2_copy.type(), ParticleType::find(0x2212));

  // add processes

  double elastic_parameter = 0;  // don't include elastic scattering
  ReactionsBitSet included_2to2 = ReactionsBitSet();

  act->add_all_scatterings(Test::default_finder_parameters(
      elastic_parameter, NNbarTreatment::NoAnnihilation, included_2to2, true,
      false, false));

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
  const auto& proton = ParticleType::find(pdg::p);
  const auto& pi_z = ParticleType::find(pdg::pi_z);
  const auto& pi_m = ParticleType::find(pdg::pi_m);
  const auto& pi_p = ParticleType::find(pdg::pi_p);
  const auto& K_p = ParticleType::find(pdg::K_p);
  const auto& K_m = ParticleType::find(pdg::K_m);
  const std::vector<std::pair<const ParticleType&, const ParticleType&>> pairs =
      {{proton, proton}, {proton, pi_z}, {proton, pi_m}, {proton, pi_p},
       {pi_p, pi_m},     {proton, K_m},  {proton, K_p}};
  for (const auto& p : pairs) {
    ParticleData p1{p.first};
    ParticleData p2{p.second};

    // set up particles
    p1.set_4position(pos_a);
    p2.set_4position(pos_b);
    const double m1 = p1.pole_mass();
    const double m2 = p2.pole_mass();
    // 0.2 GeV above the threshold, the cross section should be non-zero without
    // strings. At much higher energies, it is expected to be zero.
    const double sqrts = m1 + m2 + 0.2;
    const double p_x = plab_from_s(sqrts * sqrts, m1, m2);
    p1.set_4momentum(p1.pole_mass(), p_x, 0., 0.);
    p2.set_4momentum(p2.pole_mass(), 0., 0., 0.);
    std::cout << "Testing following pair:\n" << p1 << "\n" << p2 << std::endl;

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

    act = std::make_unique<ScatterAction>(p1_copy, p2_copy, 0.2, false, 1.0);
    VERIFY(act != nullptr);

    // add processes
    double elastic_parameter = 0;  // don't include elastic scattering
    bool strings_switch = false;

    act->add_all_scatterings(Test::default_finder_parameters(
        elastic_parameter, NNbarTreatment::NoAnnihilation,
        Test::all_reactions_included(), strings_switch, false, false));

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
  act.add_all_scatterings(Test::default_finder_parameters());

  // change the position of one of the particles
  const FourVector new_position(0.1, 0., 0., 0.);
  particles.front().set_4position(new_position);

  // update the action
  act.update_incoming(particles);
  COMPARE(act.incoming_particles()[0].position(), new_position);
}

static bool collisionbranches_equal(const CollisionBranchPtr& b1,
                                    const CollisionBranchPtr& b2) {
  bool same_particle_number = b1->particle_number() == b2->particle_number();
  bool same_weight = (std::abs(b1->weight() - b2->weight()) < really_small);
  if (b1->get_type() == ProcessType::StringSoftSingleDiffractiveAX) {
    return same_weight && same_particle_number &&
           b2->get_type() == ProcessType::StringSoftSingleDiffractiveXB;
  } else if (b1->get_type() == ProcessType::StringSoftSingleDiffractiveXB) {
    return same_weight && same_particle_number &&
           b2->get_type() == ProcessType::StringSoftSingleDiffractiveAX;
  } else {
    return same_weight && same_particle_number &&
           b1->get_type() == b2->get_type();
  }
}

TEST(particle_ordering) {
  // test if particle order in a scatteraction matters.
  //
  // create a list of all particles
  const auto& all_types = ParticleType::list_all();
  int ntypes = all_types.size();
  int64_t seed = random::generate_63bit_seed();
  random::set_seed(seed);
  for (int i = 0; i < 42 + 1; i++) {
    // create a random pair of particles
    ParticleData p1{ParticleType::find(
        all_types[random::uniform_int(0, ntypes - 1)].pdgcode())};
    ParticleData p2{ParticleType::find(
        all_types[random::uniform_int(0, ntypes - 1)].pdgcode())};

    constexpr double p_x = 3.0;
    // set momenta
    p1.set_4momentum(p1.pole_mass(), p_x, 0., 0.);
    p2.set_4momentum(p2.pole_mass(), -p_x, 0., 0.);

    // put in particles object
    Particles particles;
    particles.insert(p1);
    particles.insert(p2);

    // construct actions
    ScatterActionPtr act12, act21;
    act12 = std::make_unique<ScatterAction>(p1, p2, 0.2, false, 1.0);
    act21 = std::make_unique<ScatterAction>(p2, p1, 0.2, false, 1.0);
    auto string_process_interface = Test::default_string_process_interface();
    act12->set_string_interface(string_process_interface.get());
    act21->set_string_interface(string_process_interface.get());
    VERIFY(act12 != nullptr);
    VERIFY(act21 != nullptr);

    // add processes
    double elastic_parameter = -10.;  // no added elastic cross sections
    bool use_AQM = true;
    bool strings_with_probability = true;
    NNbarTreatment nnbar_treatment = NNbarTreatment::Strings;
    ScatterActionsFinderParameters finder_params =
        Test::default_finder_parameters(elastic_parameter, nnbar_treatment,
                                        Test::all_reactions_included(), true,
                                        use_AQM, strings_with_probability);

    act12->add_all_scatterings(finder_params);
    act21->add_all_scatterings(finder_params);

    VERIFY(act12->cross_section() >= 0.);
    VERIFY(act21->cross_section() >= 0.);

    const auto& branch12 = act12->collision_channels();
    const auto& branch21 = act21->collision_channels();
    VERIFY(branch12.size() == branch21.size());
    for (const CollisionBranchPtr& branchptr1 : branch12) {
      VERIFY(std::find_if(branch21.begin(), branch21.end(),
                          [&branchptr1](const CollisionBranchPtr& branchptr2) {
                            return collisionbranches_equal(branchptr1,
                                                           branchptr2);
                          }) != branch21.end());
    }
  }
}

TEST_CATCH(set_parametrized_total_bottomup, std::logic_error) {
  ParticleData particle{ParticleType::find(0x111)};  // pi0
  ScatterActionPtr act_bottomup = std::make_unique<ScatterAction>(
      particle, particle, 0.1, false, 1.0, -1.0, false);
  act_bottomup->set_parametrized_total_cross_section(
      Test::default_finder_parameters());
}

TEST_CATCH(add_branches_only_once, std::logic_error) {
  ParticleData particle{ParticleType::find(0x111)};  // pi0
  ScatterActionPtr act_bottomup = std::make_unique<ScatterAction>(
      particle, particle, 0.1, false, 1.0, -1.0, false);
  auto finder_parameters_bottomup = Test::default_finder_parameters();

  act_bottomup->add_all_scatterings(finder_parameters_bottomup);
  act_bottomup->add_all_scatterings(finder_parameters_bottomup);
}

TEST(top_down_sum_matches_parametrization) {
  const auto& all_types = ParticleType::list_all();
  int ntypes = all_types.size();
  int64_t seed = random::generate_63bit_seed();
  random::set_seed(seed);
  for (int i = 0; i < 42; i++) {
    // create a random pair of particles
    ParticleData p1{ParticleType::find(
        all_types[random::uniform_int(0, ntypes - 1)].pdgcode())};
    ParticleData p2{ParticleType::find(
        all_types[random::uniform_int(0, ntypes - 1)].pdgcode())};
    // set position
    p1.set_4position(pos_a);
    p2.set_4position(pos_b);
    // set momenta
    double p_x = random::uniform(0., 0.5);
    p1.set_4momentum(p1.pole_mass(), p_x, 0., 0.);
    p2.set_4momentum(p2.pole_mass(), -p_x, 0., 0.);
    std::cout << "Testing following pair:\n" << p1 << "\n" << p2 << std::endl;

    ScatterActionPtr act_topdown =
        std::make_unique<ScatterAction>(p1, p2, 0.1, false, 1.0, -1.0, true);
    VERIFY(std::isnan(act_topdown->cross_section()));

    auto string_process_interface = Test::default_string_process_interface();
    act_topdown->set_string_interface(string_process_interface.get());

    auto finder_parameters_topdown = Test::default_finder_parameters(
        0, NNbarTreatment::Strings, Test::all_reactions_included(), true, true,
        true, TotalCrossSectionStrategy::TopDown);
    act_topdown->set_parametrized_total_cross_section(
        finder_parameters_topdown);

    // fuzzyness needs to be slightly larger than 1 ulp for some compilers
    vir::test::setFuzzyness<double>(5);
    if (p1.is_hadron() && p2.is_hadron()) {
      VERIFY(act_topdown->cross_section() > really_small);

      act_topdown->add_all_scatterings(finder_parameters_topdown);

      double sum_partials = 0;
      for (const auto& proc : act_topdown->collision_channels()) {
        sum_partials += proc->weight();
      }
      FUZZY_COMPARE(sum_partials, act_topdown->cross_section());
    } else {
      FUZZY_COMPARE(act_topdown->cross_section(), 0.);
    }
  }
}
