/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "setup.h"
#include <array>
#include <string>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <include/config.h>

#include "../include/configuration.h"
#include "../include/outputinterface.h"
#include "../include/oscaroutput.h"
#include "../include/particles.h"
#include "../include/processbranch.h"
#include "../include/random.h"
#include "../include/scatteraction.h"

using namespace Smash;

static const double accuracy = 1.0e-4;
static const bf::path testoutputpath = bf::absolute(SMASH_TEST_OUTPUT_PATH);
static auto random_value = Random::make_uniform_distribution(-15.0, +15.0);

TEST(directory_is_created) {
  bf::create_directories(testoutputpath);
  VERIFY(bf::exists(testoutputpath));
}

static void compare_fourvector(const std::array<std::string, 4> &stringarray,
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

static void compare_particledata(const std::array<std::string, 12> &datastring,
                                 const ParticleData &particle, const int id) {
  COMPARE(std::atoi(datastring.at(0).c_str()), id);
  COMPARE(datastring.at(1), Test::smashon_pdg_string);
  COMPARE(std::atoi(datastring.at(2).c_str()), 0);
  std::array<std::string, 4> momentum_string;
  for (int i = 0; i < 4; i++) {
    momentum_string.at(i) = datastring.at(i + 3);
  }
  compare_fourvector(momentum_string, particle.momentum());
  COMPARE(double(std::atof(datastring.at(7).c_str())), Test::smashon_mass);
  std::array<std::string, 4> position_string;
  for (int i = 0; i < 4; i++) {
    position_string.at(i) = datastring.at(i + 8);
  }
  compare_fourvector(position_string, particle.position());
}

TEST(fullhistory_format) {
  // Set options
  const bf::path configfilename = "oscar_1999.yaml";
  const bf::path configfilepath = testoutputpath / configfilename;
  bf::ofstream(configfilepath) << "Oscar_Collisions:\n"
                                  "    Enable:          True\n"
                                  "    Print_Start_End: True\n"
                                  "    2013_Format:     False\n";
  VERIFY(bf::exists(configfilepath));

  std::unique_ptr<OutputInterface> oscfull = create_oscar_output(
      testoutputpath, Configuration{testoutputpath, configfilename});
  VERIFY(bool(oscfull));

  const bf::path outputfilepath = testoutputpath / "full_event_history.oscar";
  VERIFY(bf::exists(outputfilepath));

  Test::create_smashon_particletypes();

  Particles particles;
  const ParticleData p1 = particles.insert(Test::smashon_random());
  const ParticleData p2 = particles.insert(Test::smashon_random());

  int event_id = 0;
  /* Initial state output */
  oscfull->at_eventstart(particles, event_id);

  /* Create elastic interaction (smashon + smashon). */
  ScatterActionPtr action = make_unique<ScatterAction>(p1, p2, 0.);
  action->add_all_processes(10., true, 0., true,
                            NNbarTreatment::NoAnnihilation);
  action->generate_final_state();
  ParticleList final_particles = action->outgoing_particles();
  oscfull->at_interaction(*action, 0.);

  /* Final state output */
  oscfull->at_eventend(particles, event_id);

  // const std::string outputfilename
  //    = (testoutputpath / outputfilename).native();
  bf::fstream outputfile;
  outputfile.open(outputfilepath, std::ios_base::in);
  if (outputfile.good()) {
    std::string line, item;
    /* Check header */
    std::string output_header = "";
    std::string header =
        "# OSC1999A\n"
        "# full_event_history\n"
        "# " VERSION_MAJOR
        "\n"
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
    COMPARE(std::stoul(item), 2u);
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), event_id + 1);
    /* Check initial particle data lines item by item */
    for (const ParticleData &data : action->incoming_particles()) {
      std::array<std::string, 12> datastring;
      for (int j = 0; j < 12; j++) {
        outputfile >> datastring.at(j);
      }
      compare_particledata(datastring, data, data.id());
    }
    /* Check interaction block */
    outputfile >> item;
    COMPARE(std::stoul(item), 2u);
    outputfile >> item;
    // Additional fields are allowed: take rest of the line
    std::getline(outputfile, line);
    COMPARE(std::stoul(item), final_particles.size());
    for (const ParticleData &data : action->incoming_particles()) {
      std::array<std::string, 12> datastring;
      for (int j = 0; j < 12; j++) {
        outputfile >> datastring.at(j);
      }
      compare_particledata(datastring, data, data.id());
    }
    for (ParticleData &data : final_particles) {
      std::array<std::string, 12> datastring;
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
    for (ParticleData &data : particles) {
      std::array<std::string, 12> datastring;
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
  outputfile.close();
  VERIFY(bf::remove(outputfilepath));
}

TEST(particlelist_format) {
  // Set options
  const bf::path configfilename = "oscar_1999.yaml";
  const bf::path configfilepath = testoutputpath / configfilename;
  bf::ofstream(configfilepath) << "Oscar_Particlelist:\n"
                                  "    Enable:          True\n"
                                  "    Only_Final:      True\n"
                                  "    2013_Format:     False\n";
  VERIFY(bf::exists(configfilepath));

  std::unique_ptr<OutputInterface> oscfinal = create_oscar_output(
      testoutputpath, Configuration{testoutputpath, configfilename});
  VERIFY(bool(oscfinal));

  const bf::path outputfilepath = testoutputpath / "particle_lists.oscar";
  VERIFY(bf::exists(outputfilepath));

  Particles particles;
  /* Create 2 particles */
  const ParticleData p1 = particles.insert(Test::smashon_random());
  const ParticleData p2 = particles.insert(Test::smashon_random());

  int event_id = 0;

  /* Initial state output (note that this should not do anything!) */
  oscfinal->at_eventstart(particles, event_id);

  /* Create interaction ("elastic scattering") */
  ScatterActionPtr action = make_unique<ScatterAction>(p1, p2, 0.);
  action->add_all_processes(10., true, 0., true,
                            NNbarTreatment::NoAnnihilation);
  action->generate_final_state();

  /* As with initial state output, this should not do anything */
  oscfinal->at_interaction(*action, 0.);

  /* Final state output; this is the only thing we expect to find in file */
  action->perform(&particles, 1);
  oscfinal->at_eventend(particles, event_id);

  bf::fstream outputfile;
  outputfile.open(outputfilepath, std::ios_base::in);
  if (outputfile.good()) {
    std::string line, item;
    /* Check header */
    std::string output_header = "";
    std::string header =
        "# OSC1999A\n"
        "# final_id_p_x\n"
        "# " VERSION_MAJOR
        "\n"
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
    for (ParticleData &data : particles) {
      std::array<std::string, 12> datastring;
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
  outputfile.close();
  VERIFY(bf::remove(outputfilepath));
  VERIFY(bf::remove(configfilepath));
}
