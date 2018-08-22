/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "setup.h"

#include <smash/config.h>
#include <array>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <string>
#include <vector>

#include "../include/smash/configuration.h"
#include "../include/smash/oscaroutput.h"
#include "../include/smash/outputinterface.h"
#include "../include/smash/particles.h"
#include "../include/smash/processbranch.h"
#include "../include/smash/random.h"
#include "../include/smash/scatteraction.h"

using namespace smash;

static const double accuracy = 1.0e-4;
static const bf::path testoutputpath = bf::absolute(SMASH_TEST_OUTPUT_PATH);
static auto random_value = random::make_uniform_distribution(-15.0, +15.0);

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
  /* Create elastic interaction (smashon + smashon). */
  Test::create_smashon_particletypes();
  Particles particles;
  const ParticleData p1 = particles.insert(Test::smashon_random());
  const ParticleData p2 = particles.insert(Test::smashon_random());
  ScatterActionPtr action = make_unique<ScatterAction>(p1, p2, 0.);
  action->add_all_scatterings(10., true, Test::all_reactions_included(), 0.,
                              true, false, false,
                              NNbarTreatment::NoAnnihilation);
  action->generate_final_state();
  const ParticleList final_particles = action->outgoing_particles();

  const int event_id = 0;
  const double impact_parameter = 3.7;

  const bf::path outputfilepath =
      testoutputpath / "full_event_history.oscar1999";
  bf::path outputfilepath_unfinished = outputfilepath;
  outputfilepath_unfinished += ".unfinished";

  {
    OutputParameters out_par = OutputParameters();
    out_par.coll_printstartend = true;
    out_par.coll_extended = false;

    std::unique_ptr<OutputInterface> oscfull =
        create_oscar_output("Oscar1999", "Collisions", testoutputpath, out_par);
    VERIFY(bool(oscfull));
    VERIFY(bf::exists(outputfilepath_unfinished));

    /* Initial state output */
    oscfull->at_eventstart(particles, event_id);

    oscfull->at_interaction(*action, 0.);

    /* Final state output */
    oscfull->at_eventend(particles, event_id, impact_parameter);
  }
  VERIFY(!bf::exists(outputfilepath_unfinished));
  VERIFY(bf::exists(outputfilepath));

  {
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
          "# End of event: 0 0 event_number impact_parameter\n"
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
      for (const ParticleData &data : final_particles) {
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
      outputfile >> item;
      COMPARE(std::stod(item.c_str()), impact_parameter);
    }
  }
  VERIFY(bf::remove(outputfilepath));
}

TEST(particlelist_format) {
  /* Create 2 particles */
  Particles particles;
  const ParticleData p1 = particles.insert(Test::smashon_random());
  const ParticleData p2 = particles.insert(Test::smashon_random());

  /* Create interaction ("elastic scattering") */
  ScatterActionPtr action = make_unique<ScatterAction>(p1, p2, 0.);
  action->add_all_scatterings(10., true, Test::all_reactions_included(), 0.,
                              true, false, false,
                              NNbarTreatment::NoAnnihilation);
  action->generate_final_state();
  const int event_id = 0;
  double impact_parameter = 2.4;

  const bf::path outputfilepath = testoutputpath / "particle_lists.oscar1999";
  bf::path outputfilepath_unfinished = outputfilepath;
  outputfilepath_unfinished += ".unfinished";
  {
    OutputParameters out_par = OutputParameters();
    out_par.part_only_final = true;
    out_par.part_extended = false;

    std::unique_ptr<OutputInterface> oscfinal =
        create_oscar_output("Oscar1999", "Particles", testoutputpath, out_par);
    VERIFY(bool(oscfinal));
    VERIFY(bf::exists(outputfilepath_unfinished));

    /* Initial state output (note that this should not do anything!) */
    oscfinal->at_eventstart(particles, event_id);

    /* As with initial state output, this should not do anything */
    oscfinal->at_interaction(*action, 0.);

    /* Final state output; this is the only thing we expect to find in file */
    action->perform(&particles, 1);
    oscfinal->at_eventend(particles, event_id, impact_parameter);
  }
  VERIFY(!bf::exists(outputfilepath_unfinished));
  VERIFY(bf::exists(outputfilepath));

  {
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
          "# End of event: 0 0 event_number impact_parameter\n"
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
      outputfile >> item;
      COMPARE(std::stod(item.c_str()), impact_parameter);
    }
  }
  VERIFY(bf::remove(outputfilepath));
}
