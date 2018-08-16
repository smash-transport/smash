/*
 *
 *    Copyright (c) 2017-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "setup.h"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <string>

#include "../include/smash/listmodus.h"
#include "../include/smash/oscaroutput.h"
#include "../include/smash/particles.h"

using namespace smash;
static const double accuracy = 5.e-4;
static const bf::path testoutputpath = bf::absolute(SMASH_TEST_OUTPUT_PATH);

static bf::path create_particlefile(
    const OutputParameters out_par, const int file_number,
    std::vector<ParticleList> &init_particle_vec,
    const int particles_per_event = 10, const int n_events = 1) {
  std::unique_ptr<OutputInterface> osc2013final =
      create_oscar_output("Oscar2013", "Particles", testoutputpath, out_par);
  VERIFY(bool(osc2013final));

  const bf::path outputfilename = "particle_lists.oscar";
  const bf::path outputfilepath = testoutputpath / outputfilename;

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
    osc2013final->at_eventend(particles, event, impact_parameter);
  }

  // release and let destructor rename the file
  osc2013final.reset();

  VERIFY(bf::exists(outputfilepath));
  // Rename the oscar file to match listmodus format
  const std::string pathstring = "event" + std::to_string(file_number);
  const bf::path listinputfile = pathstring;
  const bf::path inputfilepath = testoutputpath / listinputfile;
  std::rename(outputfilepath.native().c_str(), inputfilepath.native().c_str());

  VERIFY(bf::exists(inputfilepath));

  return inputfilepath;
}

static void create_non_oscar_particlefile(
    const int file_number, std::vector<ParticleList> &init_particle_vec) {
  // Write oscar output, but write only one event and remove all comment lines.
  // This mimics the output of some hydro codes
  OutputParameters out_par = OutputParameters();
  out_par.part_only_final = true;
  out_par.part_extended = false;
  auto input_path =
      create_particlefile(out_par, file_number, init_particle_vec);

  bf::fstream file{input_path};
  const bf::path tmp_file_str = "tmp_file";
  const bf::path tmp_path = testoutputpath / tmp_file_str;

  bf::ofstream tmp_file;
  tmp_file.open(tmp_path);
  VERIFY(tmp_file.good());
  std::string line;
  while (getline(file, line)) {
    if (line[0] != '#') {
      tmp_file << line << '\n';
    }
  }
  std::rename(tmp_path.native().c_str(), input_path.native().c_str());
}

TEST(directory_is_created) {
  bf::create_directories(testoutputpath);
  VERIFY(bf::exists(testoutputpath));
}

TEST(create_particle_types) { Test::create_smashon_particletypes(); }

static void compare_fourvector(const FourVector &a, const FourVector &b) {
  COMPARE_ABSOLUTE_ERROR(a.x0(), b.x0(), accuracy);
  COMPARE_ABSOLUTE_ERROR(a.x1(), b.x1(), accuracy);
  COMPARE_ABSOLUTE_ERROR(a.x2(), b.x2(), accuracy);
  COMPARE_ABSOLUTE_ERROR(a.x3(), b.x3(), accuracy);
}

TEST(list_from_non_oscar_output) {
  std::vector<ParticleList> init_particles;
  create_non_oscar_particlefile(0, init_particles);

  // Create list modus
  std::string list_conf_str = "List:\n";
  list_conf_str += "    File_Directory: \"";
  list_conf_str += testoutputpath.native() + "\"\n";
  list_conf_str += "    File_Prefix: \"event\"\n";
  list_conf_str += "    Shift_Id: 0\n";
  auto config = Configuration(list_conf_str.c_str());
  auto par = Test::default_parameters();
  ListModus list_modus(config, par);

  // Read the file with list modus
  Particles particles_read;
  list_modus.initial_conditions(&particles_read, par);

  /*
  std::cout << "Particles from list modus:" << std::endl;
  for (const auto &p : particles_read) {
    std::cout << p << std::endl;
  }
  */

  // Scroll particles back to the earliest time, as list modus is supposed to do
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

TEST(multiple_file_non_oscar_output) {
  // Create list modus
  std::string list_conf_str = "List:\n";
  list_conf_str += "    File_Directory: \"";
  list_conf_str += testoutputpath.native() + "\"\n";
  list_conf_str += "    File_Prefix: \"event\"\n";
  list_conf_str += "    Shift_Id: 0\n";
  auto config = Configuration(list_conf_str.c_str());
  auto par = Test::default_parameters();
  ListModus list_modus(config, par);

  std::vector<ParticleList> init_particles;

  constexpr size_t max_events = 10;
  for (size_t i = 0; i < max_events; i++) {
    create_non_oscar_particlefile(i, init_particles);
  }
  COMPARE(init_particles.size(), max_events);

  // Read particles with list modus
  for (size_t current_event = 0; current_event < max_events; current_event++) {
    Particles particles_read;
    list_modus.initial_conditions(&particles_read, par);

    /*
    std::cout << "Particles from list modus:" << std::endl;
    for (const auto &p : particles_read) {
      std::cout << p << std::endl;
    }
    */

    // Scroll particles back to the earliest time, as list modus is supposed to
    // do
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
  out_par.part_only_final = true;
  out_par.part_extended = false;

  std::vector<ParticleList> init_particles;
  create_particlefile(out_par, 0, init_particles, 10, 1);

  // Create list modus
  std::string list_conf_str = "List:\n";
  list_conf_str += "    File_Directory: \"";
  list_conf_str += testoutputpath.native() + "\"\n";
  list_conf_str += "    File_Prefix: \"event\"\n";
  list_conf_str += "    Shift_Id: 0\n";
  auto config = Configuration(list_conf_str.c_str());
  auto par = Test::default_parameters();
  ListModus list_modus(config, par);

  // Read the file with list modus
  Particles particles_read;
  list_modus.initial_conditions(&particles_read, par);

  /*
  std::cout << "Particles from list modus:" << std::endl;
  for (const auto &p : particles_read) {
    std::cout << p << std::endl;
  }
  */

  // Scroll particles back to the earliest time, as list modus is supposed to do
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
  out_par.part_only_final = true;
  out_par.part_extended = false;

  std::vector<ParticleList> init_particles;
  constexpr int events_per_file = 1;
  constexpr int particles_per_event = 2;
  constexpr int n_files = 2;

  for (int i = 0; i < n_files; i++) {
    create_particlefile(out_par, i, init_particles, particles_per_event,
                        events_per_file);
  }
  // Create list modus
  std::string list_conf_str = "List:\n";
  list_conf_str += "    File_Directory: \"";
  list_conf_str += testoutputpath.native() + "\"\n";
  list_conf_str += "    File_Prefix: \"event\"\n";
  list_conf_str += "    Shift_Id: 0\n";
  auto config = Configuration(list_conf_str.c_str());
  auto par = Test::default_parameters();
  ListModus list_modus(config, par);

  for (int i = 0; i < events_per_file * n_files; i++) {
    Particles particles_read;
    list_modus.initial_conditions(&particles_read, par);

    // Scroll particles back to the earliest time, as list modus is supposed to
    // do
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
  out_par.part_only_final = true;
  out_par.part_extended = false;

  std::vector<ParticleList> init_particles;
  constexpr int events_per_file = 5;
  constexpr int particles_per_event = 10;
  constexpr int n_files = 5;

  for (int i = 0; i < n_files; i++) {
    create_particlefile(out_par, i, init_particles, particles_per_event,
                        events_per_file);
  }
  // Create list modus
  std::string list_conf_str = "List:\n";
  list_conf_str += "    File_Directory: \"";
  list_conf_str += testoutputpath.native() + "\"\n";
  list_conf_str += "    File_Prefix: \"event\"\n";
  list_conf_str += "    Shift_Id: 0\n";
  auto config = Configuration(list_conf_str.c_str());
  auto par = Test::default_parameters();
  ListModus list_modus(config, par);

  for (int i = 0; i < events_per_file * n_files; i++) {
    Particles particles_read;
    list_modus.initial_conditions(&particles_read, par);

    // Scroll particles back to the earliest time, as list modus is supposed to
    // do
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
  out_par.part_only_final = true;
  out_par.part_extended = false;

  constexpr int max_events = 2;
  constexpr int particles_per_event = 10;

  std::vector<ParticleList> init_particles;

  create_particlefile(out_par, 0, init_particles, particles_per_event,
                      max_events);

  // Create list modus
  std::string list_conf_str = "List:\n";
  list_conf_str += "    File_Directory: \"";
  list_conf_str += testoutputpath.native() + "\"\n";
  list_conf_str += "    File_Prefix: \"event\"\n";
  list_conf_str += "    Shift_Id: 0\n";
  auto config = Configuration(list_conf_str.c_str());
  auto par = Test::default_parameters();
  ListModus list_modus(config, par);

  for (int cur_event = 0; cur_event < max_events; cur_event++) {
    // Read the file with list modus
    Particles particles_read;
    list_modus.initial_conditions(&particles_read, par);

    /*
    std::cout << "Particles from list modus:" << std::endl;
    for (const auto &p : particles_read) {
      std::cout << p << std::endl;
    }
    */

    // Scroll particles back to the earliest time, as list modus is supposed
    // to do
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
  // Create list modus
  std::string list_conf_str = "List:\n";
  list_conf_str += "    File_Directory: \"";
  list_conf_str += testoutputpath.native() + "\"\n";
  list_conf_str += "    File_Prefix: \"event\"\n";
  list_conf_str += "    Shift_Id: 0\n";
  auto config = Configuration(list_conf_str.c_str());
  auto par = Test::default_parameters();
  ListModus list_modus(config, par);

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

  // Create smashons off-shell and check if they are returned on-shell
  particles.reset();
  plist_init.clear();
  plist_fin.clear();
  for (int i = 0; i < npart; i++) {
    ParticleData smashon = Test::smashon_random();
    plist_init.push_back(smashon);
    FourVector r = smashon.position(), p = smashon.momentum();
    PdgCode pdg = smashon.pdgcode();
    list_modus.try_create_particle(particles, pdg, r.x0(), r.x1(), r.x2(),
                                   r.x3(), m0 + 0.1, p.x0(), p.x1(), p.x2(),
                                   p.x3());
  }
  plist_fin = particles.copy_to_vector();
  for (int i = 0; i < npart; i++) {
    ParticleData a = plist_init.back();
    ParticleData b = plist_fin.back();
    plist_init.pop_back();
    plist_fin.pop_back();
    // Test smashon should be on mass shell with pole mass
    COMPARE_ABSOLUTE_ERROR(a.momentum().abs(), m0, accuracy);
    // The smashon read from the list should have the mass m created by user,
    // but still obey E^2 - p^2 = m^2.
    COMPARE_ABSOLUTE_ERROR(b.momentum().abs(), m0 + 0.1, accuracy);
    compare_fourvector(a.position(), b.position());
    COMPARE_ABSOLUTE_ERROR(b.formation_time(), a.position().x0(), accuracy);
    COMPARE(a.pdgcode(), b.pdgcode());
  }
}
