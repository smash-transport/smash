/*
 *
 *    Copyright (c) 2014-2020,2022-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include <array>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include "setup.h"
#include "smash/config.h"
#include "smash/fluidizationaction.h"
#include "smash/oscaroutput.h"
#include "smash/outputinterface.h"
#include "smash/particles.h"
#include "smash/processbranch.h"
#include "smash/random.h"
#include "smash/scatteraction.h"

using namespace smash;

static const double accuracy = 1.0e-4;
static const int data_elements = 12;
static const int data_elements_extended = 22;
static const std::filesystem::path testoutputpath =
    std::filesystem::absolute(SMASH_TEST_OUTPUT_PATH);
static auto random_value = random::make_uniform_distribution(-15.0, +15.0);
static const EventLabel event_id = {0, 0};

TEST(directory_is_created) {
  std::filesystem::create_directories(testoutputpath);
  VERIFY(std::filesystem::exists(testoutputpath));
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
  COMPARE_ABSOLUTE_ERROR(std::stod(datastring.at(13)),
                         particle.formation_time(), accuracy);
  COMPARE(std::stod(datastring.at(14)), particle.xsec_scaling_factor());
  COMPARE(std::atoi(datastring.at(15).c_str()), static_cast<int>(h.id_process));
  COMPARE(std::atoi(datastring.at(16).c_str()),
          static_cast<int>(h.process_type));
  COMPARE_ABSOLUTE_ERROR(std::stod(datastring.at(17)), h.time_last_collision,
                         accuracy);
  COMPARE(datastring.at(18), h.p1.string());
  COMPARE(datastring.at(19), h.p2.string());
  COMPARE(std::atoi(datastring.at(20).c_str()),
          particle.type().baryon_number());
  COMPARE(std::atoi(datastring.at(21).c_str()), particle.type().strangeness());
  // COMPARE(std::atoi(datastring.at(22).c_str()), particle.spin_projection());
}

TEST(full2013_format) {
  /* Create elastic interaction (smashon + smashon). */
  Test::create_smashon_particletypes();
  Particles particles;
  const ParticleData p1 = particles.insert(Test::smashon_random());
  const ParticleData p2 = particles.insert(Test::smashon_random());
  ScatterActionPtr action = std::make_unique<ScatterAction>(p1, p2, 0.);
  action->add_all_scatterings(Test::default_finder_parameters());
  action->generate_final_state();
  ParticleList final_particles = action->outgoing_particles();
  const double impact_parameter = 1.783;
  const bool empty_event = false;
  EventInfo event = Test::default_event_info(impact_parameter, empty_event);

  const std::filesystem::path outputfilename = "full_event_history.oscar";
  const std::filesystem::path outputfilepath = testoutputpath / outputfilename;
  std::filesystem::path outputfilepath_unfinished = outputfilepath;
  outputfilepath_unfinished += ".unfinished";
  {
    OutputParameters out_par = OutputParameters();
    out_par.coll_printstartend = true;
    out_par.coll_extended = false;

    std::unique_ptr<OutputInterface> osc2013full =
        create_oscar_output("Oscar2013", "Collisions", testoutputpath, out_par);
    VERIFY(bool(osc2013full));
    VERIFY(std::filesystem::exists(outputfilepath_unfinished));

    osc2013full->at_eventstart(particles, event_id, event);
    osc2013full->at_interaction(*action, 0.);
    action->perform(&particles, 1);
    osc2013full->at_eventend(particles, event_id, event);
  }
  VERIFY(!std::filesystem::exists(outputfilepath_unfinished));
  VERIFY(std::filesystem::exists(outputfilepath));

  {
    std::fstream outputfile;
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
      COMPARE(line, "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e");
      std::getline(outputfile, line);
      COMPARE(line, "# " SMASH_VERSION);
      /* Check initial particle list description line */
      std::string initial_line =
          "# event " + std::to_string(event_id.event_number) + " ensemble " +
          std::to_string(event_id.ensemble_number) + " in " + std::to_string(2);
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
      std::string final_line =
          "# event " + std::to_string(event_id.event_number) + " ensemble " +
          std::to_string(event_id.ensemble_number) + " out " +
          std::to_string(particles.size());
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
      std::string end_line =
          "# event " + std::to_string(event_id.event_number) + " ensemble " +
          std::to_string(event_id.ensemble_number) + " end 0" +
          " impact   1.783 scattering_projectile_target yes";
      COMPARE(line, end_line);
    }
  }
  VERIFY(std::filesystem::remove(outputfilepath));
}

TEST(final2013_format) {
  // Set options
  const std::filesystem::path configfilename = "oscar_2013.yaml";
  const std::filesystem::path configfilepath = testoutputpath / configfilename;
  std::ofstream(configfilepath) << "    Only_Final:      Yes\n";
  VERIFY(std::filesystem::exists(configfilepath));

  /* Create 2 particles */
  Particles particles;
  const ParticleData p1 = particles.insert(Test::smashon_random());
  const ParticleData p2 = particles.insert(Test::smashon_random());
  const double impact_parameter = 2.34;
  const bool empty_event = true;
  EventInfo event = Test::default_event_info(impact_parameter, empty_event);

  /* Create interaction ("elastic scattering") */
  ScatterActionPtr action = std::make_unique<ScatterAction>(p1, p2, 0.);
  action->add_all_scatterings(Test::default_finder_parameters());
  action->generate_final_state();

  const std::filesystem::path outputfilename = "particle_lists.oscar";
  const std::filesystem::path outputfilepath = testoutputpath / outputfilename;
  std::filesystem::path outputfilepath_unfinished = outputfilepath;
  outputfilepath_unfinished += ".unfinished";
  {
    OutputParameters out_par = OutputParameters();
    out_par.part_only_final = OutputOnlyFinal::Yes;
    out_par.part_extended = false;

    std::unique_ptr<OutputInterface> osc2013final =
        create_oscar_output("Oscar2013", "Particles", testoutputpath, out_par);
    VERIFY(bool(osc2013final));
    VERIFY(std::filesystem::exists(outputfilepath_unfinished));
    /* Initial state output (note that this should not do anything!) */
    osc2013final->at_eventstart(particles, event_id, event);
    /* As with initial state output, this should not do anything */
    osc2013final->at_interaction(*action, 0.);
    /* Final state output; this is the only thing we expect to find in file */
    action->perform(&particles, 1);
    osc2013final->at_eventend(particles, event_id, event);
  }
  VERIFY(!std::filesystem::exists(outputfilepath_unfinished));
  VERIFY(std::filesystem::exists(outputfilepath));

  COMPARE(action->outgoing_particles(), particles.copy_to_vector());

  {
    std::fstream outputfile;
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
      COMPARE(line, "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e");
      std::getline(outputfile, line);
      COMPARE(line, "# " SMASH_VERSION);
      /* Check final particle list */
      std::getline(outputfile, line);
      std::string final_line =
          "# event " + std::to_string(event_id.event_number) + " ensemble " +
          std::to_string(event_id.ensemble_number) + " out " +
          std::to_string(particles.size());
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
      std::string end_line =
          "# event " + std::to_string(event_id.event_number) + " ensemble " +
          std::to_string(event_id.ensemble_number) + " end 0" +
          " impact   2.340 scattering_projectile_target no";
      COMPARE(line, end_line);
    }
  }
  VERIFY(std::filesystem::remove(outputfilepath));
}

TEST(full_extended_oscar) {
  const std::filesystem::path outputfilename = "full_event_history.oscar";
  const std::filesystem::path outputfilepath = testoutputpath / outputfilename;
  std::filesystem::path outputfilepath_unfinished = outputfilepath;
  outputfilepath_unfinished += ".unfinished";

  /* Create elastic interaction (smashon + smashon). */
  Particles particles;
  const ParticleData p1 = particles.insert(Test::smashon_random());
  const ParticleData p2 = particles.insert(Test::smashon_random());
  ScatterActionPtr action = std::make_unique<ScatterAction>(p1, p2, 0.);
  action->add_all_scatterings(Test::default_finder_parameters());
  action->generate_final_state();
  ParticleList final_particles = action->outgoing_particles();
  const double impact_parameter = 1.783;
  const bool empty_event = false;
  EventInfo event = Test::default_event_info(impact_parameter, empty_event);

  {
    OutputParameters out_par = OutputParameters();
    out_par.coll_printstartend = true;
    out_par.coll_extended = true;

    std::unique_ptr<OutputInterface> osc2013full =
        create_oscar_output("Oscar2013", "Collisions", testoutputpath, out_par);
    VERIFY(bool(osc2013full));
    VERIFY(std::filesystem::exists(outputfilepath_unfinished));

    /* Initial state output */
    osc2013full->at_eventstart(particles, event_id, event);
    osc2013full->at_interaction(*action, 0.);
    /* Final state output */
    action->perform(&particles, 1);
    osc2013full->at_eventend(particles, event_id, event);
  }
  VERIFY(!std::filesystem::exists(outputfilepath_unfinished));
  VERIFY(std::filesystem::exists(outputfilepath));

  {
    std::fstream outputfile;
    outputfile.open(outputfilepath, std::ios_base::in);
    VERIFY(outputfile.good());
    std::string line;
    /* Check header */
    std::getline(outputfile, line);
    COMPARE(
        line,
        "#!OSCAR2013Extended full_event_history"
        " t x y z mass p0 px py pz pdg ID charge ncoll"
        " form_time xsecfac proc_id_origin proc_type_origin"
        " time_last_coll pdg_mother1 pdg_mother2 baryon_number strangeness");
    std::getline(outputfile, line);
    COMPARE(line,
            "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none"
            " e none fm none none none fm none none none none");
    std::getline(outputfile, line);
    COMPARE(line, "# " SMASH_VERSION);
    /* Check initial particle list description line */
    std::string initial_line =
        "# event " + std::to_string(event_id.event_number) + " ensemble " +
        std::to_string(event_id.ensemble_number) + " in " + std::to_string(2);
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
    std::string final_line =
        "# event " + std::to_string(event_id.event_number) + " ensemble " +
        std::to_string(event_id.ensemble_number) + " out " +
        std::to_string(particles.size());
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
    std::string end_line = "# event " + std::to_string(event_id.event_number) +
                           " ensemble " +
                           std::to_string(event_id.ensemble_number) + " end 0" +
                           " impact   1.783 scattering_projectile_target yes";
    COMPARE(line, end_line);
  }
  VERIFY(std::filesystem::remove(outputfilepath));
}

TEST(initial_conditions_2013_format) {
  // Create 1 particle
  Particles particles;
  ParticleData p1 = particles.insert(Test::smashon_random());
  p1.set_4position(FourVector(2.3, 1.35722, 1.42223, 1.5));  // tau = 1.74356

  // Create action ("hypersurface crossing")
  ActionPtr action = std::make_unique<FluidizationAction>(p1, p1, 0.0);
  action->generate_final_state();

  const bool empty_event = false;
  const double impact_parameter = 1.783;
  EventInfo event = Test::default_event_info(impact_parameter, empty_event);

  const std::filesystem::path outputfilename = "SMASH_IC.oscar";
  const std::filesystem::path outputfilepath = testoutputpath / outputfilename;
  std::filesystem::path outputfilepath_unfinished = outputfilepath;
  outputfilepath_unfinished += ".unfinished";
  {
    OutputParameters out_par = OutputParameters();
    out_par.ic_extended = false;

    std::unique_ptr<OutputInterface> osc2013full = create_oscar_output(
        "Oscar2013", "Initial_Conditions", testoutputpath, out_par);
    VERIFY(bool(osc2013full));
    VERIFY(std::filesystem::exists(outputfilepath_unfinished));

    osc2013full->at_eventstart(particles, event_id, event);
    action->perform(&particles, 1);
    osc2013full->at_interaction(*action, 0.);
    osc2013full->at_eventend(particles, event_id, event);
  }
  VERIFY(!std::filesystem::exists(outputfilepath_unfinished));
  VERIFY(std::filesystem::exists(outputfilepath));

  {
    std::fstream outputfile;
    outputfile.open(outputfilepath, std::ios_base::in);
    VERIFY(outputfile.good());
    if (outputfile.good()) {
      std::string line;
      /* Check header */
      std::getline(outputfile, line);
      COMPARE(line,
              "#!OSCAR2013 SMASH_IC t x y z mass p0 px py pz"
              " pdg ID charge");
      std::getline(outputfile, line);
      COMPARE(line, "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e");
      std::getline(outputfile, line);
      COMPARE(line, "# " SMASH_VERSION);
      /* Check initial particle list description line */
      std::string initial_line =
          "# event " + std::to_string(event_id.event_number) + " ensemble " +
          std::to_string(event_id.ensemble_number) + " start";
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
      /* Check for event end line */
      std::getline(outputfile, line);
      std::string end_line =
          "# event " + std::to_string(event_id.event_number) + " ensemble " +
          std::to_string(event_id.ensemble_number) + " end";
      COMPARE(line, end_line);
    }
  }
  VERIFY(std::filesystem::remove(outputfilepath));
}
