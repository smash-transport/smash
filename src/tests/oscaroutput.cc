/*
 *
 *    Copyright (c) 2014-2015
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
#include <string>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <include/config.h>
#include "../include/outputinterface.h"
#include "../include/oscaroutput.h"
#include "../include/particles.h"
#include "../include/processbranch.h"
#include "../include/random.h"
#include "../include/configuration.h"

using namespace Smash;

static const float accuracy = 1.0e-4;
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
  particle.set_4momentum(mass_smashon, random_value(), random_value(),
                         random_value());
  particle.set_4position(FourVector(random_value(), random_value(),
                                    random_value(), random_value()));
  return particle;
}

static void compare_fourvector(const std::array<std::string,4> &stringarray,
                               const FourVector &fourvector) {
  COMPARE_ABSOLUTE_ERROR(std::atof(stringarray.at(0).c_str()), fourvector.x1(),
                         accuracy);
  COMPARE_ABSOLUTE_ERROR(std::atof(stringarray.at(1).c_str()), fourvector.x2(),
                         accuracy);
  COMPARE_ABSOLUTE_ERROR(std::atof(stringarray.at(2).c_str()), fourvector.x3(),
                         accuracy);
  COMPARE_ABSOLUTE_ERROR(std::atof(stringarray.at(3).c_str()), fourvector.x0(),
                         accuracy);
}

static void compare_particledata(const std::array<std::string,12> &datastring,
                                 const ParticleData &particle, const int id) {
  COMPARE(std::atoi(datastring.at(0).c_str()), id);
  COMPARE(datastring.at(1), pdg_str);
  COMPARE(std::atoi(datastring.at(2).c_str()), zero);
  std::array<std::string,4> momentum_string;
  for (int i = 0; i < 4 ; i++) {
    momentum_string.at(i) = datastring.at(i + 3);
  }
  compare_fourvector(momentum_string, particle.momentum());
  COMPARE(float(std::atof(datastring.at(7).c_str())), mass_smashon);
  std::array<std::string,4> position_string;
  for (int i = 0; i < 4 ; i++) {
    position_string.at(i) = datastring.at(i + 8);
  }
  compare_fourvector(position_string, particle.position());
}

TEST(fullhistory_format) {
  // Set options
  std::string configfilename = "oscar_1999.yaml";
  bf::ofstream(testoutputpath / configfilename)
      << "Oscar_Collisions:\n"
         "    Enable:          True\n"
         "    Print_Start_End: True\n"
         "    2013_Format:     False\n";
  VERIFY(bf::exists(testoutputpath / configfilename));

  std::unique_ptr<OutputInterface> oscfull = create_oscar_output(
      testoutputpath, Configuration{testoutputpath, configfilename});
  VERIFY(bool(oscfull));

  std::string outputfilename = "full_event_history.oscar";
  VERIFY(bf::exists(testoutputpath / outputfilename));

  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n" + smashon_str);

  Particles particles;

  ParticleData particle = create_smashon_particle();
  particles.add_data(particle);

  ParticleData second_particle = create_smashon_particle();
  particles.add_data(second_particle);

  int event_id = 0;
  /* Initial state output */
  oscfull->at_eventstart(particles, event_id);
  /* Create interaction ("resonance formation") */
  std::vector<ParticleData> initial_particles, final_particles;
  initial_particles.push_back(particles.data(0));
  initial_particles.push_back(particles.data(1));
  particles.remove(0);
  particles.remove(1);
  ParticleData final_particle = create_smashon_particle();
  particles.add_data(final_particle);
  final_particles.push_back(particles.data(particles.id_max()));
  oscfull->at_interaction(initial_particles, final_particles, 0.0, 0.0,
      ProcessType::None);
  /* Final state output */
  oscfull->at_eventend(particles, event_id);

  std::fstream outputfile;
  outputfile.open((testoutputpath / outputfilename)
                  .native().c_str(), std::ios_base::in);
  if (outputfile.good()) {
    std::string line, item;
    /* Check header */
    std::string output_header = "";
    std::string header = "# OSC1999A\n"
                         "# full_event_history\n"
                         "# " VERSION_MAJOR "\n"
                         "# Block format:\n"
                         "# nin nout event_number\n"
                         "# id pdg 0 px py pz p0 mass x y z t\n"
                         "# End of event: 0 0 event_number\n"
                         "#\n";
    do {
      std::getline(outputfile, line);
      output_header += line + '\n';
    } while (line != "#");
    COMPARE(output_header, header);
    /* Check initial particle list description line item by item */
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), 0);
    outputfile >> item;
    COMPARE(std::stoul(item), initial_particles.size());
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), event_id + 1);
    /* Check initial particle data lines item by item */
    for (ParticleData &data : initial_particles) {
      std::array<std::string,12> datastring;
      for (int j = 0; j < 12; j++) {
        outputfile >> datastring.at(j);
      }
      compare_particledata(datastring, data, data.id());
    }
    /* Check interaction block */
    outputfile >> item;
    COMPARE(std::stoul(item), initial_particles.size());
    outputfile >> item;
    // Additional fields are allowed: take rest of the line
    std::getline(outputfile, line);
    COMPARE(std::stoul(item), final_particles.size());
    for (ParticleData &data : initial_particles) {
      std::array<std::string,12> datastring;
      for (int j = 0; j < 12; j++) {
        outputfile >> datastring.at(j);
      }
      compare_particledata(datastring, data, data.id());
    }
    for (ParticleData &data : final_particles) {
      std::array<std::string,12> datastring;
      for (int j = 0; j < 12; j++) {
        outputfile >> datastring.at(j);
      }
      compare_particledata(datastring, data, data.id());
    }
    /* Check final particle list */
    outputfile >> item;
    COMPARE(std::stoul(item), final_particles.size());
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), 0);
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), event_id + 1);
    for (ParticleData &data : particles.data()) {
      std::array<std::string,12> datastring;
      for (int j = 0; j < 12; j++) {
        outputfile >> datastring.at(j);
      }
      compare_particledata(datastring, data, data.id());
    }
    /* Check for null interaction */
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), 0);
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), 0);
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), event_id + 1);
  }
}


TEST(particlelist_format) {
  // Set options
  std::string configfilename = "oscar_1999.yaml";
  bf::ofstream(testoutputpath / configfilename)
      << "Oscar_Particlelist:\n"
         "    Enable:          True\n"
         "    Only_Final:      True\n"
         "    2013_Format:     False\n";
  VERIFY(bf::exists(testoutputpath / configfilename));

  std::unique_ptr<OutputInterface> oscfinal = create_oscar_output(
      testoutputpath, Configuration{testoutputpath, configfilename});
  VERIFY(bool(oscfinal));

  std::string outputfilename = "particle_lists.oscar";
  VERIFY(bf::exists(testoutputpath / outputfilename));

  Particles particles;
  /* Create 5 particles */
  for (int i = 0; i < 5; i++) {
    ParticleData particle = create_smashon_particle();
    particles.add_data(particle);
  }
  int event_id = 0;

  /* Initial state output (note that this should not do anything!) */
  oscfinal->at_eventstart(particles, event_id);
  /* Create interaction ("elastic scattering") */
  std::vector<ParticleData> initial_particles, final_particles;
  initial_particles.push_back(particles.data(0));
  initial_particles.push_back(particles.data(1));
  /* Change the momenta */
  particles.data(0).set_4momentum(mass_smashon, random_value(),
                                  random_value(), random_value());
  particles.data(1).set_4momentum(mass_smashon, random_value(),
                                  random_value(), random_value());
  final_particles.push_back(particles.data(0));
  final_particles.push_back(particles.data(1));
  /* As with initial state output, this should not do anything */
  oscfinal->at_interaction(initial_particles, final_particles, 0.0, 0.0,
      ProcessType::None);
  /* Final state output; this is the only thing we expect to find in file */
  oscfinal->at_eventend(particles, event_id);

  std::fstream outputfile;
  outputfile.open((testoutputpath / outputfilename)
                  .native().c_str(), std::ios_base::in);
  if (outputfile.good()) {
    std::string line, item;
    /* Check header */
    std::string output_header = "";
    std::string header = "# OSC1999A\n"
                         "# final_id_p_x\n"
                         "# " VERSION_MAJOR "\n"
                         "# Block format:\n"
                         "# nin nout event_number\n"
                         "# id pdg 0 px py pz p0 mass x y z t\n"
                         "# End of event: 0 0 event_number\n"
                         "#\n";
    do {
      std::getline(outputfile, line);
      output_header += line + '\n';
    } while (line != "#");
    COMPARE(output_header, header);
    /* Check final particle list */
    outputfile >> item;
    COMPARE(std::stoul(item), particles.size());
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), 0);
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), event_id + 1);
    for (ParticleData &data : particles.data()) {
      std::array<std::string,12> datastring;
      for (int j = 0; j < 12; j++) {
        outputfile >> datastring.at(j);
      }
      compare_particledata(datastring, data, data.id());
    }
    /* Check for null interaction */
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), 0);
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), 0);
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), event_id + 1);
  }
}
