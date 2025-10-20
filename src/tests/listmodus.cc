/*
 *
 *    Copyright (c) 2017-2020,2022-2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/listmodus.h"

#include <filesystem>
#include <string>

#include "setup.h"
#include "smash/oscaroutput.h"
#include "smash/particles.h"

using namespace smash;
static const double accuracy = 5.e-4;
static const std::filesystem::path testoutputpath =
    std::filesystem::absolute(SMASH_TEST_OUTPUT_PATH);
static const auto parameters = Test::default_parameters();

static std::filesystem::path create_particlefile(
    const OutputParameters out_par, const int file_number,
    std::vector<ParticleList> &init_particle_vec,
    const int particles_per_event = 10, const int n_events = 1) {
  std::unique_ptr<OutputInterface> osc2013final =
      create_oscar_output("Oscar2013", "Particles", testoutputpath, out_par);
  VERIFY(bool(osc2013final));

  const std::filesystem::path outputfilename = "particle_lists.oscar";
  const std::filesystem::path outputfilepath = testoutputpath / outputfilename;

  // Create random particles
  for (int event = 0; event < n_events; event++) {
    Particles particles;
    for (int i = 0; i < particles_per_event; i++) {
      particles.insert(Test::smashon_random());
    }

    /*
    std::cout << "Initial particles:" << std::endl;
    for (const auto &p : particles) {
      std::cout << p << std::endl;
    }
    */

    init_particle_vec.push_back(particles.copy_to_vector());

    // Print them to file in OSCAR 2013 format
    const double impact_parameter = 2.34;  // just a dummy value here
    const bool empty_event = false;        // just a dummy value as well
    EventInfo default_event_info =
        Test::default_event_info(impact_parameter, empty_event);
    osc2013final->at_eventend(particles, {event, 0}, default_event_info);
  }

  // release and let destructor rename the file
  osc2013final.reset();

  VERIFY(std::filesystem::exists(outputfilepath));
  // Rename the oscar file to match listmodus format
  const std::filesystem::path inputfilepath =
      testoutputpath / ("event" + std::to_string(file_number));
  std::filesystem::rename(outputfilepath, inputfilepath);

  VERIFY(std::filesystem::exists(inputfilepath));

  return inputfilepath;
}

static void create_non_oscar_particlefile(
    const int file_number, std::vector<ParticleList> &init_particle_vec) {
  // Write oscar output, but write only one event and remove all comment lines.
  // This mimics the output of some hydro codes
  OutputParameters out_par = OutputParameters();
  out_par.part_only_final = OutputOnlyFinal::Yes;
  out_par.part_extended = false;
  auto input_path =
      create_particlefile(out_par, file_number, init_particle_vec);

  std::fstream file{input_path};
  const std::filesystem::path tmp_file_str = "tmp_file";
  const std::filesystem::path tmp_path = testoutputpath / tmp_file_str;

  std::ofstream tmp_file;
  tmp_file.open(tmp_path);
  VERIFY(tmp_file.good());
  std::string line;
  while (getline(file, line)) {
    if (line[0] != '#') {
      tmp_file << line << '\n';
    }
  }
  std::filesystem::rename(tmp_path, input_path);
}

static void create_particle_file_with_multiple_particles_at_same_position(
    const OutputParameters out_par) {
  std::unique_ptr<OutputInterface> oscar_output =
      create_oscar_output("Oscar2013", "Particles", testoutputpath, out_par);

  const std::filesystem::path outputfilename = "particle_lists.oscar";
  const std::filesystem::path outputfilepath = testoutputpath / outputfilename;

  // Create some random particles with many at same 4-position
  for (int event = 0; event < 2; event++) {
    Particles particles;
    for (int i = 0; i < 3; i++) {
      particles.insert(Test::smashon_random());
    }
    auto particle = Test::smashon_random();
    for (int i = 0; i < event * 3; i++) {
      particles.insert(particle);
      particle.boost_momentum({0.1, 0.2, 0.3});
    }
    particle = Test::smashon_random();
    for (int i = 0; i < event + 3; i++) {
      particles.insert(particle);
      particle.boost_momentum({0.01, 0.02, 0.03});
    }

    // Print them to file in OSCAR 2013 format
    const double impact_parameter = 2.34;  // just a dummy value here
    const bool empty_event = false;        // just a dummy value as well
    EventInfo default_event_info =
        Test::default_event_info(impact_parameter, empty_event);
    oscar_output->at_eventend(particles, {event, 0}, default_event_info);
  }

  // release and let destructor rename the file
  oscar_output.reset();

  VERIFY(std::filesystem::exists(outputfilepath));
  // Rename the oscar file to match listmodus format
  const std::filesystem::path inputfilepath = testoutputpath / "event0";
  std::filesystem::rename(outputfilepath, inputfilepath);
}

static ListModus create_list_modus_for_test() {
  Configuration config{R"(
    Modi:
      List:
        File_Directory: ToBeSet
        File_Prefix: event
    )"};
  config.set_value(InputKeys::modi_list_fileDirectory, testoutputpath.string());
  return ListModus(std::move(config), parameters);
}

static ListBoxModus create_list_box_modus_for_test() {
  Configuration config{R"(
    Modi:
      ListBox:
        File_Directory: ToBeSet
        File_Prefix: event
        Length: 3
    )"};
  config.set_value(InputKeys::modi_listBox_fileDirectory,
                   testoutputpath.string());
  return ListBoxModus(std::move(config), parameters);
}

static ListModus create_list_modus_with_single_file_for_test() {
  Configuration config{R"(
    Modi:
      List:
        File_Directory: ToBeSet
        Filename: event0
    )"};
  config.set_value(InputKeys::modi_list_fileDirectory, testoutputpath.string());
  return ListModus(std::move(config), parameters);
}

TEST(directory_is_created) {
  std::filesystem::create_directories(testoutputpath);
  VERIFY(std::filesystem::exists(testoutputpath));
}

TEST(create_particle_types) { Test::create_stable_smashon_particletypes(); }

static void compare_fourvector(const FourVector &a, const FourVector &b) {
  COMPARE_ABSOLUTE_ERROR(a.x0(), b.x0(), accuracy);
  COMPARE_ABSOLUTE_ERROR(a.x1(), b.x1(), accuracy);
  COMPARE_ABSOLUTE_ERROR(a.x2(), b.x2(), accuracy);
  COMPARE_ABSOLUTE_ERROR(a.x3(), b.x3(), accuracy);
}

TEST(list_from_non_oscar_output) {
  std::vector<ParticleList> init_particles;
  create_non_oscar_particlefile(0, init_particles);
  ListModus list_modus = create_list_modus_for_test();

  // Read the file with list modus
  Particles particles_read;
  list_modus.initial_conditions(&particles_read, parameters);

  /*
  std::cout << "Particles from list modus:" << std::endl;
  for (const auto &p : particles_read) {
    std::cout << p << std::endl;
  }
  */

  // Scroll particles back to the earliest time, as list modus should do
  double earliest_t = 1.e8;
  for (const auto &particle : init_particles[0]) {
    if (particle.position().x0() < earliest_t) {
      earliest_t = particle.position().x0();
    }
  }
  for (auto &particle : init_particles[0]) {
    const double t = particle.position().x0();
    const FourVector u(1.0, particle.velocity());
    particle.set_formation_time(t);
    particle.set_4position(particle.position() + u * (earliest_t - t));
  }

  COMPARE(particles_read.size(), init_particles[0].size());
  ParticleList p_init = init_particles[0];
  ParticleList p_fin = particles_read.copy_to_vector();
  for (size_t i = 0; i < p_fin.size(); i++) {
    ParticleData a = p_init.back();
    ParticleData b = p_fin.back();
    p_init.pop_back();
    p_fin.pop_back();
    compare_fourvector(a.momentum(), b.momentum());
    compare_fourvector(a.position(), b.position());
    COMPARE(a.id(), b.id());
    COMPARE_ABSOLUTE_ERROR(a.formation_time(), b.formation_time(), accuracy);
    COMPARE(a.pdgcode(), b.pdgcode());
  }
}

TEST(listbox_creation_from_non_oscar_output) {
  std::vector<ParticleList> init_particles;
  create_non_oscar_particlefile(0, init_particles);
  ListBoxModus list_box_modus = create_list_box_modus_for_test();

  // Read the file with list modus
  Particles particles_read;
  list_box_modus.initial_conditions(&particles_read, parameters);
}

TEST(multiple_file_non_oscar_output) {
  ListModus list_modus = create_list_modus_for_test();

  std::vector<ParticleList> init_particles;

  constexpr size_t max_events = 10;
  for (size_t i = 0; i < max_events; i++) {
    create_non_oscar_particlefile(i, init_particles);
  }
  COMPARE(init_particles.size(), max_events);

  // Read particles with list modus
  for (size_t current_event = 0; current_event < max_events; current_event++) {
    Particles particles_read;
    list_modus.initial_conditions(&particles_read, parameters);

    /*
    std::cout << "Particles from list modus:" << std::endl;
    for (const auto &p : particles_read) {
      std::cout << p << std::endl;
    }
    */

    // Scroll particles back to the earliest time, as list modus should do
    double earliest_t = 1.e8;
    for (const auto &particle : init_particles[current_event]) {
      if (particle.position().x0() < earliest_t) {
        earliest_t = particle.position().x0();
      }
    }
    for (auto &particle : init_particles[current_event]) {
      const double t = particle.position().x0();
      const FourVector u(1.0, particle.velocity());
      particle.set_formation_time(t);
      particle.set_4position(particle.position() + u * (earliest_t - t));
    }

    COMPARE(particles_read.size(), init_particles[current_event].size());
    ParticleList p_init = init_particles[current_event];
    ParticleList p_fin = particles_read.copy_to_vector();
    for (size_t i = 0; i < p_fin.size(); i++) {
      ParticleData a = p_init.back();
      ParticleData b = p_fin.back();
      p_init.pop_back();
      p_fin.pop_back();
      compare_fourvector(a.momentum(), b.momentum());
      compare_fourvector(a.position(), b.position());
      COMPARE(a.id(), b.id());
      COMPARE_ABSOLUTE_ERROR(a.formation_time(), b.formation_time(), accuracy);
      COMPARE(a.pdgcode(), b.pdgcode());
    }
  }
}

TEST(list_from_oscar2013_output) {
  // Create OSCAR 2013 output
  OutputParameters out_par = OutputParameters();
  out_par.part_only_final = OutputOnlyFinal::Yes;
  out_par.part_extended = false;
  std::vector<ParticleList> init_particles;
  create_particlefile(out_par, 0, init_particles, 10, 1);
  ListModus list_modus = create_list_modus_for_test();

  // Read the file with list modus
  Particles particles_read;
  list_modus.initial_conditions(&particles_read, parameters);

  /*
  std::cout << "Particles from list modus:" << std::endl;
  for (const auto &p : particles_read) {
    std::cout << p << std::endl;
  }
  */

  // Scroll particles back to the earliest time, as list modus should do
  double earliest_t = 1.e8;
  for (const auto &particle : init_particles[0]) {
    if (particle.position().x0() < earliest_t) {
      earliest_t = particle.position().x0();
    }
  }
  for (auto &particle : init_particles[0]) {
    const double t = particle.position().x0();
    const FourVector u(1.0, particle.velocity());
    particle.set_formation_time(t);
    particle.set_4position(particle.position() + u * (earliest_t - t));
  }

  COMPARE(particles_read.size(), init_particles[0].size());
  ParticleList p_init = init_particles[0];
  ParticleList p_fin = particles_read.copy_to_vector();
  for (size_t i = 0; i < p_fin.size(); i++) {
    ParticleData a = p_init.back();
    ParticleData b = p_fin.back();
    p_init.pop_back();
    p_fin.pop_back();
    compare_fourvector(a.momentum(), b.momentum());
    compare_fourvector(a.position(), b.position());
    COMPARE(a.id(), b.id());
    COMPARE_ABSOLUTE_ERROR(a.formation_time(), b.formation_time(), accuracy);
    COMPARE(a.pdgcode(), b.pdgcode());
  }
}

TEST(multiple_files_one_event) {
  OutputParameters out_par = OutputParameters();
  out_par.part_only_final = OutputOnlyFinal::Yes;
  out_par.part_extended = false;
  std::vector<ParticleList> init_particles;
  constexpr int events_per_file = 1;
  constexpr int particles_per_event = 2;
  constexpr int n_files = 2;
  for (int i = 0; i < n_files; i++) {
    create_particlefile(out_par, i, init_particles, particles_per_event,
                        events_per_file);
  }
  ListModus list_modus = create_list_modus_for_test();

  for (int i = 0; i < events_per_file * n_files; i++) {
    Particles particles_read;
    list_modus.initial_conditions(&particles_read, parameters);

    // Scroll particles back to the earliest time, as list modus should do
    double earliest_t = 1.e8;
    for (const auto &particle : init_particles[i]) {
      if (particle.position().x0() < earliest_t) {
        earliest_t = particle.position().x0();
      }
    }
    for (auto &particle : init_particles[i]) {
      const double t = particle.position().x0();
      const FourVector u(1.0, particle.velocity());
      particle.set_formation_time(t);
      particle.set_4position(particle.position() + u * (earliest_t - t));
    }

    COMPARE(particles_read.size(), init_particles[i].size());
    ParticleList p_init = init_particles[i];
    ParticleList p_fin = particles_read.copy_to_vector();
    for (size_t j = 0; j < p_fin.size(); j++) {
      ParticleData a = p_init.back();
      ParticleData b = p_fin.back();
      p_init.pop_back();
      p_fin.pop_back();
      compare_fourvector(a.momentum(), b.momentum());
      compare_fourvector(a.position(), b.position());
      COMPARE(a.id(), b.id());
      COMPARE_ABSOLUTE_ERROR(a.formation_time(), b.formation_time(), accuracy);
      COMPARE(a.pdgcode(), b.pdgcode());
    }
  }
}

TEST(multiple_files_multiple_events) {
  OutputParameters out_par = OutputParameters();
  out_par.part_only_final = OutputOnlyFinal::Yes;
  out_par.part_extended = false;

  std::vector<ParticleList> init_particles;
  constexpr int events_per_file = 5;
  constexpr int particles_per_event = 10;
  constexpr int n_files = 5;

  for (int i = 0; i < n_files; i++) {
    create_particlefile(out_par, i, init_particles, particles_per_event,
                        events_per_file);
  }
  ListModus list_modus = create_list_modus_for_test();

  for (int i = 0; i < events_per_file * n_files; i++) {
    Particles particles_read;
    list_modus.initial_conditions(&particles_read, parameters);

    // Scroll particles back to the earliest time, as ListModus should do
    double earliest_t = 1.e8;
    for (const auto &particle : init_particles[i]) {
      if (particle.position().x0() < earliest_t) {
        earliest_t = particle.position().x0();
      }
    }
    for (auto &particle : init_particles[i]) {
      const double t = particle.position().x0();
      const FourVector u(1.0, particle.velocity());
      particle.set_formation_time(t);
      particle.set_4position(particle.position() + u * (earliest_t - t));
    }

    COMPARE(particles_read.size(), init_particles[i].size());
    ParticleList p_init = init_particles[i];
    ParticleList p_fin = particles_read.copy_to_vector();
    for (size_t j = 0; j < p_fin.size(); j++) {
      ParticleData a = p_init.back();
      ParticleData b = p_fin.back();
      p_init.pop_back();
      p_fin.pop_back();
      compare_fourvector(a.momentum(), b.momentum());
      compare_fourvector(a.position(), b.position());
      COMPARE(a.id(), b.id());
      COMPARE_ABSOLUTE_ERROR(a.formation_time(), b.formation_time(), accuracy);
      COMPARE(a.pdgcode(), b.pdgcode());
    }
  }
}

TEST(multiple_events_in_file) {
  OutputParameters out_par = OutputParameters();
  out_par.part_only_final = OutputOnlyFinal::Yes;
  out_par.part_extended = false;
  constexpr int max_events = 2;
  constexpr int particles_per_event = 10;
  std::vector<ParticleList> init_particles;
  create_particlefile(out_par, 0, init_particles, particles_per_event,
                      max_events);
  ListModus list_modus = create_list_modus_with_single_file_for_test();

  for (int cur_event = 0; cur_event < max_events; cur_event++) {
    // Read the file with list modus
    Particles particles_read;
    list_modus.initial_conditions(&particles_read, parameters);

    /*
    std::cout << "Particles from list modus:" << std::endl;
    for (const auto &p : particles_read) {
      std::cout << p << std::endl;
    }
    */

    // Scroll particles back to the earliest time, as list modus should do
    double earliest_t = 1.e8;
    for (const auto &particle : init_particles[cur_event]) {
      if (particle.position().x0() < earliest_t) {
        earliest_t = particle.position().x0();
      }
    }
    for (auto &particle : init_particles[cur_event]) {
      const double t = particle.position().x0();
      const FourVector u(1.0, particle.velocity());
      particle.set_formation_time(t);
      particle.set_4position(particle.position() + u * (earliest_t - t));
    }

    COMPARE(particles_read.size(), init_particles[cur_event].size());
    ParticleList p_init = init_particles[cur_event];
    ParticleList p_fin = particles_read.copy_to_vector();
    for (size_t i = 0; i < p_fin.size(); i++) {
      ParticleData a = p_init.back();
      ParticleData b = p_fin.back();
      p_init.pop_back();
      p_fin.pop_back();
      compare_fourvector(a.momentum(), b.momentum());
      compare_fourvector(a.position(), b.position());
      COMPARE(a.id(), b.id());
      COMPARE_ABSOLUTE_ERROR(a.formation_time(), b.formation_time(), accuracy);
      COMPARE(a.pdgcode(), b.pdgcode());
    }
  }
}

TEST(try_create_particle_func) {
  ListModus list_modus = create_list_modus_for_test();
  Particles particles;
  ParticleList plist_init, plist_fin;
  const int npart = 10;
  const double m0 = Test::smashon_mass;

  for (int i = 0; i < npart; i++) {
    ParticleData smashon = Test::smashon_random();
    plist_init.push_back(smashon);
    FourVector r = smashon.position(), p = smashon.momentum();
    PdgCode pdg = smashon.pdgcode();
    list_modus.try_create_particle(particles, pdg, r.x0(), r.x1(), r.x2(),
                                   r.x3(), m0, p.x0(), p.x1(), p.x2(), p.x3());
  }
  plist_fin = particles.copy_to_vector();
  for (int i = 0; i < npart; i++) {
    ParticleData a = plist_init.back();
    ParticleData b = plist_fin.back();
    plist_init.pop_back();
    plist_fin.pop_back();
    compare_fourvector(a.momentum(), b.momentum());
    compare_fourvector(a.position(), b.position());
    COMPARE_ABSOLUTE_ERROR(b.formation_time(), a.position().x0(), accuracy);
    COMPARE(a.pdgcode(), b.pdgcode());
  }

  // Create stable smashons some with a mass discrepancy and some off-shell
  // because of their energy. Test if they are returned on-shell with pole mass.
  particles.reset();
  plist_init.clear();
  plist_fin.clear();
  for (int i = 0; i < 2 * npart; i++) {
    ParticleData smashon = Test::smashon_random();
    plist_init.push_back(smashon);
    FourVector r = smashon.position(), p = smashon.momentum();
    PdgCode pdg = smashon.pdgcode();
    const auto creation_mass = m0 + (i < npart);
    const auto creation_energy = p.x0() + (i >= npart);
    list_modus.try_create_particle(particles, pdg, r.x0(), r.x1(), r.x2(),
                                   r.x3(), creation_mass, creation_energy,
                                   p.x1(), p.x2(), p.x3());
  }
  plist_fin = particles.copy_to_vector();
  for (int i = 0; i < npart; i++) {
    ParticleData a = plist_init.back();
    ParticleData b = plist_fin.back();
    plist_init.pop_back();
    plist_fin.pop_back();
    // Test smashon should be on mass shell with pole mass
    COMPARE_ABSOLUTE_ERROR(a.momentum().abs(), m0, accuracy);
    // The smashon read from the list should have the pole mass,
    // but still obey E^2 - p^2 = m^2.
    COMPARE_ABSOLUTE_ERROR(b.momentum().abs(), m0, accuracy);
    compare_fourvector(a.position(), b.position());
    COMPARE_ABSOLUTE_ERROR(b.formation_time(), a.position().x0(), accuracy);
    COMPARE(a.pdgcode(), b.pdgcode());
  }
}

// Test that try_create_particle correctly stores the provided spin vector
// when spin interactions are enabled.
TEST(try_create_particle_with_spin_func) {
  // Create a ListModus with spin interactions enabled
  ExperimentParameters parameters_with_spin{
      std::make_unique<UniformClock>(0., 0.1, 300.0),
      std::make_unique<UniformClock>(0., 1., 300.0),
      1,  // ensembles
      1,  // testparticles
      DerivativesMode::CovariantGaussian,
      RestFrameDensityDerivativesMode::Off,
      FieldDerivativesMode::ChainRule,
      SmearingMode::CovariantGaussian,
      1.0,
      4.0,
      0.333333,
      2.0,
      CollisionCriterion::Geometric,
      true,
      Test::all_reactions_included(),
      Test::no_multiparticle_reactions(),
      false,  // strings
      1.0,
      NNbarTreatment::NoAnnihilation,
      0.,
      false,
      -1.0,
      200.0,
      2.5,
      1.0,
      false,
      false,
      true,
      SpinInteractionType::On,
      std::nullopt};

  Configuration config{R"(
    Modi:
      List:
        File_Directory: ToBeSet
        File_Prefix: event
    )"};
  config.set_value(InputKeys::modi_list_fileDirectory, testoutputpath.string());
  config.set_value(
      InputKeys::modi_list_optionalQuantities,
      std::vector<std::string>{"spin0", "spinx", "spiny", "spinz"});
  ListModus list_modus = ListModus(std::move(config), parameters_with_spin);
  Particles particles;
  ParticleData smashon = Test::smashon_random();
  FourVector r = smashon.position(), p = smashon.momentum();
  PdgCode pdg = smashon.pdgcode();
  // Spin vector components as 4-vector and optional values
  FourVector spin_vec(0.1, 0.2, 0.3, 0.4);
  std::vector<std::string> opt_vals = {"0.1", "0.2", "0.3", "0.4"};
  list_modus.try_create_particle(particles, pdg, r.x0(), r.x1(), r.x2(), r.x3(),
                                 Test::smashon_mass, p.x0(), p.x1(), p.x2(),
                                 p.x3(), opt_vals);

  ParticleList plist = particles.copy_to_vector();
  ParticleData created = plist.back();
  compare_fourvector(created.spin_vector(), spin_vec);
}

TEST_CATCH(create_particle_with_nan, std::invalid_argument) {
  ListModus list_modus = create_list_modus_for_test();
  Particles particles;
  const double m0 = Test::smashon_mass;
  ParticleData smashon = Test::smashon_random();
  FourVector r = smashon.position(), p = smashon.momentum();
  PdgCode pdg = smashon.pdgcode();

  // Create a particle with either a NAN value in the position
  // to trigger an invalid_argument error.
  list_modus.try_create_particle(particles, pdg, NAN, r.x1(), r.x2(), r.x3(),
                                 m0, p.x0(), p.x1(), p.x2(), p.x3());
}

TEST_CATCH(create_particles_at_same_position, ListModus::InvalidEvents) {
  const OutputParameters out_par = OutputParameters();
  create_particle_file_with_multiple_particles_at_same_position(out_par);
  ListModus list_modus = create_list_modus_with_single_file_for_test();
}
