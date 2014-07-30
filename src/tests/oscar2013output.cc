/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include <cstdio>
#include <cstring>
#include <array>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <boost/filesystem.hpp>
#include "../include/outputinterface.h"
#include "../include/oscarfullhistoryoutput.h"
#include "../include/oscarparticlelistoutput.h"
#include "../include/particles.h"
#include "../include/random.h"

using namespace Smash;

static const float accuracy = 1.0e-4;
static const int data_elements = 11;
static const bf::path testoutputpath = bf::absolute(SMASH_TEST_OUTPUT_PATH);
static auto random_value = Random::make_uniform_distribution(-15.0, +15.0);

TEST(directory_is_created) {
  bf::create_directory(testoutputpath);
  VERIFY(bf::exists(testoutputpath));
}

static const float mass_smashon = 0.123;
static std::string mass_str = std::to_string(mass_smashon);
static std::string width_str = "1.200";
static std::string pdg_str = "661";
static std::string smashon_str = "smashon " + mass_str + " "
    + width_str + " " + pdg_str + "\n";
static const int zero = 0;

static ParticleData create_smashon_particle() {
  ParticleData particle = ParticleData{ParticleType::find(0x661)};
  particle.set_momentum(mass_smashon, random_value(), random_value(),
                        random_value());
  particle.set_position(FourVector(random_value(), random_value(),
                                   random_value(), random_value()));
  return particle;
}

static void compare_fourvector(const std::array<std::string,4> &stringarray,
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
  std::array<std::string,4> position_string;
  for (int i = 0; i < 4 ; i++) {
    position_string.at(i) = datastring.at(i);
  }
  compare_fourvector(position_string, particle.position());
  COMPARE(float(std::atof(datastring.at(4).c_str())), mass_smashon);
  std::array<std::string,4> momentum_string;
  for (int i = 0; i < 4 ; i++) {
    momentum_string.at(i) = datastring.at(i + 5);
  }
  compare_fourvector(momentum_string, particle.momentum());
  COMPARE(datastring.at(9), pdg_str);
  COMPARE(std::atoi(datastring.at(10).c_str()), id);
}

TEST(full2013_format) {
  // Set options
  std::string configfilename = "oscar_2013.yaml";
  std::ofstream configfile;
  configfile.open((testoutputpath / configfilename).native().c_str());
  configfile << "Print_start_end:" << "\t" << "True" << std::endl;
  configfile << "2013_format:" << "\t" << "True" << std::endl;
  configfile.close();
  VERIFY(bf::exists(testoutputpath / configfilename));
  Configuration&& op{testoutputpath, configfilename};
  OscarFullHistoryOutput *osc2013full
                   = new OscarFullHistoryOutput(testoutputpath, std::move(op));
  std::string outputfilename = "full_event_history.oscar";
  VERIFY(bf::exists(testoutputpath / outputfilename));

  std::cout << "Testing full format" << std::endl;

  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n" + smashon_str);

  Particles particles;

  ParticleData particle = create_smashon_particle();
  particles.add_data(particle);

  ParticleData second_particle = create_smashon_particle();
  particles.add_data(second_particle);

  int event_id = 0;
  /* Initial state output */
  osc2013full->at_eventstart(particles, event_id);
  /* Create interaction ("resonance formation") */
  std::vector<ParticleData> initial_particles, final_particles;
  initial_particles.push_back(particles.data(0));
  initial_particles.push_back(particles.data(1));
  particles.remove(0);
  particles.remove(1);
  ParticleData final_particle = create_smashon_particle();
  particles.add_data(final_particle);
  final_particles.push_back(particles.data(particles.id_max()));
  osc2013full->at_interaction(initial_particles, final_particles);
  /* Final state output */
  osc2013full->at_eventend(particles, event_id);

  std::fstream outputfile;
  outputfile.open((testoutputpath / outputfilename)
                  .native().c_str(), std::ios_base::in);
  VERIFY(outputfile.good());
  if (outputfile.good()) {
    std::string line, item;
    /* Check header */
    std::string output_header = "";
    std::string header = "#!OSCAR2013 "
                         "full_event_history "
                         "t x y z mass p0 px py pz pdg ID\n"
                         "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none\n";
    do {
      std::getline(outputfile, line);
      output_header += line + '\n';
    } while (line != "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none");
    COMPARE(output_header, header);
    /* Check initial particle list description line */
    std::string initial_line = "# event " + std::to_string(event_id + 1)
      + " in " + std::to_string(initial_particles.size());
    std::getline(outputfile, line);
    COMPARE(line, initial_line);
    /* Check initial particle data lines item by item */
    for (ParticleData &data : initial_particles) {
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
    std::string interaction_line = "# interaction in "
      + std::to_string(initial_particles.size())
      + " out " + std::to_string(final_particles.size());
    COMPARE(line, interaction_line);
    for (ParticleData &data : initial_particles) {
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
    std::string final_line = "# event " + std::to_string(event_id + 1)
      + " out " + std::to_string(particles.size());
    COMPARE(line, final_line);
    for (ParticleData &data : particles.data()) {
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
    std::string end_line = "# event " + std::to_string(event_id + 1)
      + " end";
    COMPARE(line, end_line);
  }
}


TEST(final2013_format) {
  // Set options
  std::string configfilename = "oscar_2013.yaml";
  std::ofstream configfile;
  configfile.open((testoutputpath / configfilename).native().c_str(),
                  std::ios::in);
  configfile << "Only_final:" << "\t" << "True" << std::endl;
  configfile << "2013_format:" << "\t" << "True" << std::endl;
  configfile.close();
  VERIFY(bf::exists(testoutputpath / configfilename));

  Configuration&& op{testoutputpath, configfilename};

  OscarParticleListOutput *osc2013final
    = new OscarParticleListOutput(testoutputpath, std::move(op));
  std::string outputfilename = "particle_lists.oscar";
  VERIFY(bf::exists(testoutputpath / outputfilename));

  std::cout << "Testing final format" << std::endl;

  Particles particles;

  /* Create 5 particles */
  for (int i = 0; i < 5; i++) {
    ParticleData particle = create_smashon_particle();
    particles.add_data(particle);
  }
  int event_id = 0;

  /* Initial state output (note that this should not do anything!) */
  osc2013final->at_eventstart(particles, event_id);
  /* Create interaction ("elastic scattering") */
  std::vector<ParticleData> initial_particles, final_particles;
  initial_particles.push_back(particles.data(0));
  initial_particles.push_back(particles.data(1));
  /* Change the momenta */
  particles.data(0).set_momentum(mass_smashon, random_value(), random_value(),
                        random_value());
  particles.data(1).set_momentum(mass_smashon, random_value(), random_value(),
                                 random_value());
  final_particles.push_back(particles.data(0));
  final_particles.push_back(particles.data(1));
  /* As with initial state output, this should not do anything */
  osc2013final->at_interaction(initial_particles, final_particles);
  /* Final state output; this is the only thing we expect to find in file */
  osc2013final->at_eventend(particles, event_id);

  std::fstream outputfile;
  outputfile.open((testoutputpath / outputfilename)
                  .native().c_str(), std::ios_base::in);
  VERIFY(outputfile.good());
  if (outputfile.good()) {
    std::string line, item;
    /* Check header */
    std::string output_header = "";
    std::string header = "#!OSCAR2013 "
                         "particle_lists "
                         "t x y z mass p0 px py pz pdg ID\n"
                         "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none\n";
    do {
      std::getline(outputfile, line);
      output_header += line + '\n';
    } while (line != "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none");
    COMPARE(output_header, header);
    /* Check final particle list */
    std::getline(outputfile, line);
    std::string final_line = "# event " + std::to_string(event_id + 1)
      + " out " + std::to_string(particles.size());
    COMPARE(line, final_line);
    for (ParticleData &data : particles.data()) {
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
    std::string end_line = "# event " + std::to_string(event_id + 1)
      + " end";
    COMPARE(line, end_line);
  }
}
