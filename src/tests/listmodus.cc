/*
 *
 *    Copyright (c) 2017
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

#include "../include/listmodus.h"
#include "../include/oscaroutput.h"
#include "../include/particles.h"

using namespace smash;
static const double accuracy = 5.e-5;
static const bf::path testoutputpath = bf::absolute(SMASH_TEST_OUTPUT_PATH);

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

TEST(list_from_oscar2013_output) {
  // Create OSCAR 2013 output
  OutputParameters out_par = OutputParameters();
  out_par.part_only_final = true;
  out_par.part_extended = false;

  std::unique_ptr<OutputInterface> osc2013final =
      create_oscar_output("Oscar2013", "Particles", testoutputpath, out_par);
  VERIFY(bool(osc2013final));

  const bf::path outputfilename = "particle_lists.oscar";
  const bf::path outputfilepath = testoutputpath / outputfilename;
  VERIFY(bf::exists(outputfilepath));

  // Create random particles
  Particles particles;
  constexpr size_t N = 10;
  for (size_t i = 0; i < N; i++) {
    particles.insert(Test::smashon_random());
  }

  std::cout << "Initial particles:" << std::endl;
  for (const auto &p : particles) {
    std::cout << p << std::endl;
  }

  // Print them to file in OSCAR 2013 format
  const int event_id = 0;
  const double impact_parameter = 2.34;  // just a dummy value here
  osc2013final->at_eventend(particles, event_id, impact_parameter);

  // Rename the oscar file to match listmodus format
  const bf::path listinputfile = "event0";
  const bf::path inputfilepath = testoutputpath / listinputfile;
  std::rename(outputfilepath.native().c_str(), inputfilepath.native().c_str());

  // Create list modus
  std::string list_conf_str = "List:\n";
  list_conf_str += "    File_Directory: \"";
  list_conf_str += testoutputpath.native() + "\"\n";
  list_conf_str += "    File_Prefix: \"event\"\n";
  list_conf_str += "    Shift_Id: 0\n";
  list_conf_str += "    Start_Time: 0.0\n";
  auto config = Configuration(list_conf_str.c_str());
  auto par = Test::default_parameters();
  ListModus list_modus(config, par);

  // Read the file with list modus
  Particles particles_read;
  list_modus.initial_conditions(&particles_read, par);

  std::cout << "Particles from list modus:" << std::endl;
  for (const auto &p : particles_read) {
    std::cout << p << std::endl;
  }

  // Scroll particles back to the earliest time, as list modus is supposed to do
  double earliest_t = 1.e8;
  for (const auto &particle : particles) {
    if (particle.position().x0() < earliest_t) {
      earliest_t = particle.position().x0();
    }
  }
  for (auto &particle : particles) {
    const double t = particle.position().x0();
    const FourVector u(1.0, particle.velocity());
    particle.set_formation_time(t);
    particle.set_4position(particle.position() + u * (earliest_t - t));
  }

  COMPARE(particles_read.size(), particles.size());
  ParticleList p_init = particles.copy_to_vector();
  ParticleList p_fin = particles_read.copy_to_vector();
  for (size_t i = 0; i < N; i++) {
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
