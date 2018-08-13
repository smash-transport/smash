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
#include <map>
#include <string>
#include <vector>

#include "../include/smash/oscaroutput.h"
#include "../include/smash/outputinterface.h"
#include "../include/smash/particles.h"
#include "../include/smash/processbranch.h"
#include "../include/smash/random.h"
#include "../include/smash/scatteraction.h"

using namespace smash;

static const double accuracy = 1.0e-4;
static const int data_elements = 12;
static const int data_elements_extended = 20;
static const bf::path testoutputpath = bf::absolute(SMASH_TEST_OUTPUT_PATH);
static auto random_value = random::make_uniform_distribution(-15.0, +15.0);

TEST(directory_is_created) {
  bf::create_directories(testoutputpath);
  VERIFY(bf::exists(testoutputpath));
}

TEST(extended_has_more_fields) {
  VERIFY(data_elements_extended > data_elements);
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
  COMPARE(std::atoi(datastring.at(11).c_str()), particle.type().charge());
}

static void compare_extended_particledata(
    const std::array<std::string, data_elements_extended> &datastring,
    const ParticleData &particle, const int id) {
  std::array<std::string, data_elements> smaller_datastring;
  for (size_t i = 0; i < data_elements; i++) {
    smaller_datastring[i] = datastring[i];
  }
  compare_particledata(smaller_datastring, particle, id);
  const auto h = particle.get_history();
  COMPARE(std::atoi(datastring.at(12).c_str()), h.collisions_per_particle);
  COMPARE(std::atoi(datastring.at(13).c_str()), particle.formation_time());
  COMPARE(std::atoi(datastring.at(14).c_str()),
          particle.cross_section_scaling_factor());
  COMPARE(std::atoi(datastring.at(15).c_str()), static_cast<int>(h.id_process));
  COMPARE(std::atoi(datastring.at(16).c_str()),
          static_cast<int>(h.process_type));
  COMPARE_ABSOLUTE_ERROR(std::atof(datastring.at(17).c_str()),
                         h.time_last_collision, accuracy);
  COMPARE(datastring.at(18), h.p1.string());
  COMPARE(datastring.at(19), h.p2.string());
}

TEST(full2013_format) {
  /* Create elastic interaction (smashon + smashon). */
  Test::create_smashon_particletypes();
  Particles particles;
  const ParticleData p1 = particles.insert(Test::smashon_random());
  const ParticleData p2 = particles.insert(Test::smashon_random());
  ScatterActionPtr action = make_unique<ScatterAction>(p1, p2, 0.);
  action->add_all_scatterings(10., true, Test::all_reactions_included(), 0.,
                          true, false, false, NNbarTreatment::NoAnnihilation);
  action->generate_final_state();
  ParticleList final_particles = action->outgoing_particles();
  const double impact_parameter = 1.783;
  const int event_id = 0;

  const bf::path outputfilename = "full_event_history.oscar";
  const bf::path outputfilepath = testoutputpath / outputfilename;
  bf::path outputfilepath_unfinished = outputfilepath;
  outputfilepath_unfinished += ".unfinished";
  {
    OutputParameters out_par = OutputParameters();
    out_par.coll_printstartend = true;
    out_par.coll_extended = false;

    std::unique_ptr<OutputInterface> osc2013full =
        create_oscar_output("Oscar2013", "Collisions", testoutputpath, out_par);
    VERIFY(bool(osc2013full));
    VERIFY(bf::exists(outputfilepath_unfinished));

    osc2013full->at_eventstart(particles, event_id);
    osc2013full->at_interaction(*action, 0.);
    action->perform(&particles, 1);
    osc2013full->at_eventend(particles, event_id, impact_parameter);
  }
  VERIFY(!bf::exists(outputfilepath_unfinished));
  VERIFY(bf::exists(outputfilepath));

  {
    bf::fstream outputfile;
    outputfile.open(outputfilepath, std::ios_base::in);
    VERIFY(outputfile.good());
    if (outputfile.good()) {
      std::string line;
      /* Check header */
      std::getline(outputfile, line);
      COMPARE(line,
              "#!OSCAR2013 full_event_history t x y z mass p0 px py pz"
              " pdg ID charge");
      std::getline(outputfile, line);
      COMPARE(line, "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none none");
      std::getline(outputfile, line);
      COMPARE(line, "# " VERSION_MAJOR);
      /* Check initial particle list description line */
      std::string initial_line = "# event " + std::to_string(event_id + 1) +
                                 " in " + std::to_string(2);
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
      std::string end_line = "# event " + std::to_string(event_id + 1) +
                             " end 0" + " impact   1.783";
      COMPARE(line, end_line);
    }
  }
  VERIFY(bf::remove(outputfilepath));
}

TEST(final2013_format) {
  // Set options
  const bf::path configfilename = "oscar_2013.yaml";
  const bf::path configfilepath = testoutputpath / configfilename;
  bf::ofstream(configfilepath) << "    Only_Final:      True\n";
  VERIFY(bf::exists(configfilepath));

  /* Create 2 particles */
  Particles particles;
  const ParticleData p1 = particles.insert(Test::smashon_random());
  const ParticleData p2 = particles.insert(Test::smashon_random());
  const int event_id = 0;
  const double impact_parameter = 2.34;

  /* Create interaction ("elastic scattering") */
  ScatterActionPtr action = make_unique<ScatterAction>(p1, p2, 0.);
  action->add_all_scatterings(10., true, Test::all_reactions_included(), 0.,
                          true, false, false, NNbarTreatment::NoAnnihilation);
  action->generate_final_state();

  const bf::path outputfilename = "particle_lists.oscar";
  const bf::path outputfilepath = testoutputpath / outputfilename;
  bf::path outputfilepath_unfinished = outputfilepath;
  outputfilepath_unfinished += ".unfinished";
  {
    OutputParameters out_par = OutputParameters();
    out_par.part_only_final = true;
    out_par.part_extended = false;

    std::unique_ptr<OutputInterface> osc2013final =
        create_oscar_output("Oscar2013", "Particles", testoutputpath, out_par);
    VERIFY(bool(osc2013final));
    VERIFY(bf::exists(outputfilepath_unfinished));
    /* Initial state output (note that this should not do anything!) */
    osc2013final->at_eventstart(particles, event_id);
    /* As with initial state output, this should not do anything */
    osc2013final->at_interaction(*action, 0.);
    /* Final state output; this is the only thing we expect to find in file */
    action->perform(&particles, 1);
    osc2013final->at_eventend(particles, event_id, impact_parameter);
  }
  VERIFY(!bf::exists(outputfilepath_unfinished));
  VERIFY(bf::exists(outputfilepath));

  COMPARE(action->outgoing_particles(), particles.copy_to_vector());

  {
    bf::fstream outputfile;
    outputfile.open(outputfilepath, std::ios_base::in);
    VERIFY(outputfile.good());
    if (outputfile.good()) {
      std::string line;
      /* Check header */
      std::getline(outputfile, line);
      COMPARE(line,
              "#!OSCAR2013 particle_lists t x y z mass p0 px py pz"
              " pdg ID charge");
      std::getline(outputfile, line);
      COMPARE(line, "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none none");
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
      std::string end_line = "# event " + std::to_string(event_id + 1) +
                             " end 0" + " impact   2.340";
      COMPARE(line, end_line);
    }
  }
  VERIFY(bf::remove(outputfilepath));
}

TEST(full_extended_oscar) {
  const bf::path outputfilename = "full_event_history.oscar";
  const bf::path outputfilepath = testoutputpath / outputfilename;
  bf::path outputfilepath_unfinished = outputfilepath;
  outputfilepath_unfinished += ".unfinished";

  /* Create elastic interaction (smashon + smashon). */
  Particles particles;
  const ParticleData p1 = particles.insert(Test::smashon_random());
  const ParticleData p2 = particles.insert(Test::smashon_random());
  ScatterActionPtr action = make_unique<ScatterAction>(p1, p2, 0.);
  action->add_all_scatterings(10., true, Test::all_reactions_included(), 0.,
                          true, false, false, NNbarTreatment::NoAnnihilation);
  action->generate_final_state();
  ParticleList final_particles = action->outgoing_particles();
  const int event_id = 0;
  const double impact_parameter = 1.783;

  {
    OutputParameters out_par = OutputParameters();
    out_par.coll_printstartend = true;
    out_par.coll_extended = true;

    std::unique_ptr<OutputInterface> osc2013full =
        create_oscar_output("Oscar2013", "Collisions", testoutputpath, out_par);
    VERIFY(bool(osc2013full));
    VERIFY(bf::exists(outputfilepath_unfinished));

    /* Initial state output */
    osc2013full->at_eventstart(particles, event_id);
    osc2013full->at_interaction(*action, 0.);
    /* Final state output */
    action->perform(&particles, 1);
    osc2013full->at_eventend(particles, event_id, impact_parameter);
  }
  VERIFY(!bf::exists(outputfilepath_unfinished));
  VERIFY(bf::exists(outputfilepath));

  {
    bf::fstream outputfile;
    outputfile.open(outputfilepath, std::ios_base::in);
    VERIFY(outputfile.good());
    std::string line;
    /* Check header */
    std::getline(outputfile, line);
    COMPARE(line,
            "#!OSCAR2013Extended full_event_history"
            " t x y z mass p0 px py pz pdg ID charge ncoll"
            " form_time xsecfac proc_id_origin proc_type_origin"
            " time_last_coll pdg_mother1 pdg_mother2");
    std::getline(outputfile, line);
    COMPARE(line,
            "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none"
            " none none fm none none none fm none none");
    std::getline(outputfile, line);
    COMPARE(line, "# " VERSION_MAJOR);
    /* Check initial particle list description line */
    std::string initial_line =
        "# event " + std::to_string(event_id + 1) + " in " + std::to_string(2);
    std::getline(outputfile, line);
    COMPARE(line, initial_line);
    /* Check initial particle data lines item by item */
    for (const ParticleData &data : action->incoming_particles()) {
      std::array<std::string, data_elements_extended> datastring;
      for (int j = 0; j < data_elements_extended; j++) {
        outputfile >> datastring.at(j);
      }
      compare_extended_particledata(datastring, data, data.id());
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
      std::array<std::string, data_elements_extended> datastring;
      for (int j = 0; j < data_elements_extended; j++) {
        outputfile >> datastring.at(j);
      }
      compare_extended_particledata(datastring, data, data.id());
    }
    for (ParticleData &data : final_particles) {
      std::array<std::string, data_elements_extended> datastring;
      for (int j = 0; j < data_elements_extended; j++) {
        outputfile >> datastring.at(j);
      }
      compare_extended_particledata(datastring, data, data.id());
    }
    /* Get the dangling newline character */
    outputfile.get();
    /* Check final particle list */
    std::getline(outputfile, line);
    std::string final_line = "# event " + std::to_string(event_id + 1) +
                             " out " + std::to_string(particles.size());
    COMPARE(line, final_line);
    for (ParticleData &data : particles) {
      std::array<std::string, data_elements_extended> datastring;
      for (int j = 0; j < data_elements_extended; j++) {
        outputfile >> datastring.at(j);
      }
      compare_extended_particledata(datastring, data, data.id());
    }
    /* Get the dangling newline character */
    outputfile.get();
    /* Check for event end line */
    std::getline(outputfile, line);
    std::string end_line = "# event " + std::to_string(event_id + 1) +
                           " end 0" + " impact   1.783";
    COMPARE(line, end_line);
  }
  VERIFY(bf::remove(outputfilepath));
}
