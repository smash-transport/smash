/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "setup.h"

#include <include/config.h>
#include <array>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <map>
#include <string>
#include <vector>

#include "../include/oscaroutput.h"
#include "../include/outputinterface.h"
#include "../include/particles.h"
#include "../include/processbranch.h"
#include "../include/random.h"
#include "../include/scatteraction.h"

using namespace Smash;

static const double accuracy = 1.0e-4;
static const int data_elements = 11;
static const bf::path testoutputpath = bf::absolute(SMASH_TEST_OUTPUT_PATH);
static auto random_value = Random::make_uniform_distribution(-15.0, +15.0);

TEST(directory_is_created) {
  bf::create_directories(testoutputpath);
  VERIFY(bf::exists(testoutputpath));
}

static void compare_fourvector(const std::array<std::string, 4> &stringarray,
                               const FourVector &fourvector) {
  COMPARE_ABSOLUTE_ERROR(std::atof(stringarray.at(0).c_str()), fourvector.x0(),
                         accuracy);
  COMPARE_ABSOLUTE_ERROR(std::atof(stringarray.at(1).c_str()), fourvector.x1(),
                         accuracy);
  COMPARE_ABSOLUTE_ERROR(std::atof(stringarray.at(2).c_str()), fourvector.x2(),
                         accuracy);
  COMPARE_ABSOLUTE_ERROR(std::atof(stringarray.at(3).c_str()), fourvector.x3(),
                         accuracy);
}

static void compare_particledata(
    const std::array<std::string, data_elements> &datastring,
    const ParticleData &particle, const int id) {
  std::array<std::string, 4> position_string;
  for (int i = 0; i < 4; i++) {
    position_string.at(i) = datastring.at(i);
  }
  compare_fourvector(position_string, particle.position());
  COMPARE(double(std::atof(datastring.at(4).c_str())), Test::smashon_mass)
      << datastring.at(4);
  std::array<std::string, 4> momentum_string;
  for (int i = 0; i < 4; i++) {
    momentum_string.at(i) = datastring.at(i + 5);
  }
  compare_fourvector(momentum_string, particle.momentum());
  COMPARE(datastring.at(9), Test::smashon_pdg_string);
  COMPARE(std::atoi(datastring.at(10).c_str()), id);
}

TEST(full2013_format) {
  // Set options
  const bf::path configfilename = "oscar_2013.yaml";
  const bf::path configfilepath = testoutputpath / configfilename;
  bf::ofstream(configfilepath) << "Oscar_Collisions:\n"
                                  "    Enable:          True\n"
                                  "    Print_Start_End: True\n"
                                  "    2013_Format:     True\n";
  VERIFY(bf::exists(configfilepath));

  std::unique_ptr<OutputInterface> osc2013full = create_oscar_output(
      testoutputpath, Configuration{testoutputpath, configfilename});
  VERIFY(bool(osc2013full));

  const bf::path outputfilename = "full_event_history.oscar";
  const bf::path outputfilepath = testoutputpath / outputfilename;
  VERIFY(bf::exists(outputfilepath));

  Test::create_smashon_particletypes();

  Particles particles;
  const ParticleData p1 = particles.insert(Test::smashon_random());
  const ParticleData p2 = particles.insert(Test::smashon_random());

  int event_id = 0;
  /* Initial state output */
  osc2013full->at_eventstart(particles, event_id);

  /* Create elastic interaction (smashon + smashon). */
  ScatterActionPtr action = make_unique<ScatterAction>(p1, p2, 0.);
  action->add_all_processes(10., true, true, 0., true,
                            NNbarTreatment::NoAnnihilation);
  action->generate_final_state();
  ParticleList final_particles = action->outgoing_particles();
  osc2013full->at_interaction(*action, 0.);

  /* Final state output */
  action->perform(&particles, 1);
  osc2013full->at_eventend(particles, event_id);

  bf::fstream outputfile;
  outputfile.open(outputfilepath, std::ios_base::in);
  VERIFY(outputfile.good());
  if (outputfile.good()) {
    std::string line;
    /* Check header */
    std::getline(outputfile, line);
    COMPARE(line,
            "#!OSCAR2013 full_event_history t x y z mass p0 px py pz pdg ID");
    std::getline(outputfile, line);
    COMPARE(line, "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none");
    std::getline(outputfile, line);
    COMPARE(line, "# " VERSION_MAJOR);
    /* Check initial particle list description line */
    std::string initial_line =
        "# event " + std::to_string(event_id + 1) + " in " + std::to_string(2);
    std::getline(outputfile, line);
    COMPARE(line, initial_line);
    /* Check initial particle data lines item by item */
    for (const ParticleData &data : action->incoming_particles()) {
      std::array<std::string, data_elements> datastring;
      for (int j = 0; j < data_elements; j++) {
        outputfile >> datastring.at(j);
      }
      compare_particledata(datastring, data, data.id());
    }
    /* Get the dangling newline character */
    outputfile.get();
    /* Check interaction block */
    std::getline(outputfile, line);
    std::string interaction_line = "# interaction in " + std::to_string(2) +
                                   " out " +
                                   std::to_string(final_particles.size());
    // Allow additional fields
    COMPARE(line.substr(0, interaction_line.size()), interaction_line);
    for (const ParticleData &data : action->incoming_particles()) {
      std::array<std::string, data_elements> datastring;
      for (int j = 0; j < data_elements; j++) {
        outputfile >> datastring.at(j);
      }
      compare_particledata(datastring, data, data.id());
    }
    for (ParticleData &data : final_particles) {
      std::array<std::string, data_elements> datastring;
      for (int j = 0; j < data_elements; j++) {
        outputfile >> datastring.at(j);
      }
      compare_particledata(datastring, data, data.id());
    }
    /* Get the dangling newline character */
    outputfile.get();
    /* Check final particle list */
    std::getline(outputfile, line);
    std::string final_line = "# event " + std::to_string(event_id + 1) +
                             " out " + std::to_string(particles.size());
    COMPARE(line, final_line);
    for (ParticleData &data : particles) {
      std::array<std::string, data_elements> datastring;
      for (int j = 0; j < data_elements; j++) {
        outputfile >> datastring.at(j);
      }
      compare_particledata(datastring, data, data.id());
    }
    /* Get the dangling newline character */
    outputfile.get();
    /* Check for event end line */
    std::getline(outputfile, line);
    std::string end_line = "# event " + std::to_string(event_id + 1) + " end 0";
    COMPARE(line, end_line);
  }
  outputfile.close();
  VERIFY(bf::remove(outputfilepath));
  VERIFY(bf::remove(configfilepath));
}

TEST(final2013_format) {
  // Set options
  const bf::path configfilename = "oscar_2013.yaml";
  const bf::path configfilepath = testoutputpath / configfilename;
  bf::ofstream(configfilepath) << "Oscar_Particlelist:\n"
                                  "    Enable:          True\n"
                                  "    Only_Final:      True\n"
                                  "    2013_Format:     True\n";
  VERIFY(bf::exists(configfilepath));

  std::unique_ptr<OutputInterface> osc2013final = create_oscar_output(
      testoutputpath, Configuration{testoutputpath, configfilename});
  VERIFY(bool(osc2013final));

  const bf::path outputfilename = "particle_lists.oscar";
  const bf::path outputfilepath = testoutputpath / outputfilename;
  VERIFY(bf::exists(outputfilepath));

  Particles particles;

  /* Create 2 particles */
  const ParticleData p1 = particles.insert(Test::smashon_random());
  const ParticleData p2 = particles.insert(Test::smashon_random());

  int event_id = 0;

  /* Initial state output (note that this should not do anything!) */
  osc2013final->at_eventstart(particles, event_id);

  /* Create interaction ("elastic scattering") */
  ScatterActionPtr action = make_unique<ScatterAction>(p1, p2, 0.);
  action->add_all_processes(10., true, true, 0., true,
                            NNbarTreatment::NoAnnihilation);
  action->generate_final_state();

  /* As with initial state output, this should not do anything */
  osc2013final->at_interaction(*action, 0.);

  /* Final state output; this is the only thing we expect to find in file */
  action->perform(&particles, 1);
  osc2013final->at_eventend(particles, event_id);
  COMPARE(action->outgoing_particles(), particles.copy_to_vector());

  bf::fstream outputfile;
  outputfile.open(outputfilepath, std::ios_base::in);
  VERIFY(outputfile.good());
  if (outputfile.good()) {
    std::string line;
    /* Check header */
    std::getline(outputfile, line);
    COMPARE(line, "#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID");
    std::getline(outputfile, line);
    COMPARE(line, "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none");
    std::getline(outputfile, line);
    COMPARE(line, "# " VERSION_MAJOR);
    /* Check final particle list */
    std::getline(outputfile, line);
    std::string final_line = "# event " + std::to_string(event_id + 1) +
                             " out " + std::to_string(particles.size());
    COMPARE(line, final_line);
    for (ParticleData &data : particles) {
      std::array<std::string, data_elements> datastring;
      for (int j = 0; j < data_elements; j++) {
        outputfile >> datastring.at(j);
      }
      compare_particledata(datastring, data, data.id());
    }
    /* Get the dangling newline character */
    outputfile.get();
    /* Check for event end line */
    std::getline(outputfile, line);
    std::string end_line = "# event " + std::to_string(event_id + 1) + " end 0";
    COMPARE(line, end_line);
  }
  outputfile.close();
  VERIFY(bf::remove(outputfilepath));
  VERIFY(bf::remove(configfilepath));
}
