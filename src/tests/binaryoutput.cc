/*
 *
 *    Copyright (c) 2014-2020,2022-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/binaryoutput.h"

#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "setup.h"
#include "smash/clock.h"
#include "smash/config.h"
#include "smash/file.h"
#include "smash/fluidizationaction.h"
#include "smash/outputinterface.h"
#include "smash/processbranch.h"
#include "smash/scatteraction.h"
#include "smash/scatteractionsfinderparameters.h"

using namespace smash;

static const std::filesystem::path testoutputpath =
    std::filesystem::absolute(SMASH_TEST_OUTPUT_PATH);

static const EventLabel event_id = {0, 0};

TEST(directory_is_created) {
  std::filesystem::create_directories(testoutputpath);
  VERIFY(std::filesystem::exists(testoutputpath));
}

TEST(init_particletypes) { Test::create_smashon_particletypes(); }

static const uint16_t current_format_version = 10;

/* A set of convenient functions to read binary */

static void read_binary(std::string &s, const FilePtr &file) {
  size_t size = s.size();
  COMPARE(std::fread(&size, sizeof(std::int32_t), 1, file.get()), 1u);
  std::vector<char> buf(size);
  COMPARE(std::fread(&buf[0], sizeof(char), size, file.get()),
          static_cast<size_t>(size));
  s.assign(&buf[0], size);
}

static void read_binary(FourVector &v, const FilePtr &file) {
  COMPARE(std::fread(v.begin(), sizeof(*v.begin()), 4, file.get()), 4u);
}

static void read_binary(std::uint16_t &x, const FilePtr &file) {
  COMPARE(std::fread(&x, sizeof(x), 1, file.get()), 1u);
}

static void read_binary(std::int32_t &x, const FilePtr &file) {
  COMPARE(std::fread(&x, sizeof(x), 1, file.get()), 1u);
}

static void read_binary(double &x, const FilePtr &file) {
  COMPARE(std::fread(&x, sizeof(x), 1, file.get()), 1u);
}

/* Function to read and compare particle */
static bool compare_particle(const ParticleData &p, const FilePtr &file) {
  int32_t id, pdgcode, charge;
  double mass;
  FourVector pos, mom;
  read_binary(pos, file);
  read_binary(mass, file);
  read_binary(mom, file);
  read_binary(pdgcode, file);
  read_binary(id, file);
  read_binary(charge, file);

  return (p.id() == id) && (p.pdgcode().get_decimal() == pdgcode) &&
         (pos == p.position()) &&
         (mom == p.momentum() && charge == p.type().charge());
}

/* Reads and compares particle in case of extended format */
static void compare_particle_extended(const ParticleData &p,
                                      const FilePtr &file) {
  VERIFY(compare_particle(p, file));
  int32_t collisions_per_particle, id_process, process_type, p1pdg, p2pdg,
      baryon_number, strangeness;
  double formation_time, xs_scaling_factor, time_last_collision;
  const auto h = p.get_history();
  read_binary(collisions_per_particle, file);
  read_binary(formation_time, file);
  read_binary(xs_scaling_factor, file);
  read_binary(id_process, file);
  read_binary(process_type, file);
  read_binary(time_last_collision, file);
  read_binary(p1pdg, file);
  read_binary(p2pdg, file);
  read_binary(baryon_number, file);
  read_binary(strangeness, file);
  COMPARE(collisions_per_particle, h.collisions_per_particle);
  COMPARE(formation_time, p.formation_time());
  COMPARE(xs_scaling_factor, p.xsec_scaling_factor());
  COMPARE(id_process, static_cast<int>(h.id_process));
  COMPARE(process_type, static_cast<int>(h.process_type));
  COMPARE(time_last_collision, h.time_last_collision);
  COMPARE(p1pdg, h.p1.get_decimal());
  COMPARE(p2pdg, h.p2.get_decimal());
  COMPARE(baryon_number, p.type().baryon_number());
  COMPARE(strangeness, p.type().strangeness());
}

/* function to read and compare particle block header */
static bool compare_particles_block_header(const EventLabel &ev,
                                           const int32_t npart,
                                           const FilePtr &file) {
  int32_t npart_read, ev_read, ens_read;
  char c_read;
  COMPARE(std::fread(&c_read, sizeof(char), 1, file.get()), 1u);
  read_binary(ev_read, file);
  read_binary(ens_read, file);
  read_binary(npart_read, file);
  return (c_read == 'p') && (npart_read == npart) &&
         (ev_read == ev.event_number) && (ens_read == ev.ensemble_number);
}

/* function to read and compare collision block header */
static bool compare_interaction_block_header(const int32_t nin,
                                             const int32_t nout,
                                             const Action &action, double rho,
                                             const FilePtr &file) {
  int32_t nin_read, nout_read, process_type_read;
  double rho_read, weight_read, partial_weight_read;
  char c_read;
  auto process_type = static_cast<int32_t>(action.get_type());
  COMPARE(std::fread(&c_read, sizeof(char), 1, file.get()), 1u);
  read_binary(nin_read, file);
  read_binary(nout_read, file);
  COMPARE(std::fread(&rho_read, sizeof(double), 1, file.get()), 1u);
  COMPARE(std::fread(&weight_read, sizeof(double), 1, file.get()), 1u);
  COMPARE(std::fread(&partial_weight_read, sizeof(double), 1, file.get()), 1u);
  read_binary(process_type_read, file);
  return (c_read == 'i') && (nin_read == nin) && (nout_read == nout) &&
         (rho_read == rho) && (weight_read == action.get_total_weight()) &&
         (partial_weight_read == action.get_partial_weight()) &&
         (process_type_read == process_type);
}

/* function to read and compare event end line */
static bool compare_final_block_header(const EventLabel &ev,
                                       const double impact_parameter,
                                       const bool empty_event,
                                       const FilePtr &file) {
  int32_t ev_read, ens_read;
  char c_read;
  double b_read;
  char empty_event_read;
  COMPARE(std::fread(&c_read, sizeof(char), 1, file.get()), 1u);
  read_binary(ev_read, file);
  read_binary(ens_read, file);
  COMPARE(std::fread(&b_read, sizeof(double), 1, file.get()), 1u);
  COMPARE(std::fread(&empty_event_read, sizeof(char), 1, file.get()), 1u);
  return (c_read == 'f') && (ev_read == ev.event_number) &&
         (ens_read == ev.ensemble_number) && (b_read == impact_parameter) &&
         (empty_event_read == empty_event);
}

static bool compare_initial_conditions_interaction_block_header(
    const int32_t npart, const FilePtr &file) {
  int32_t npart_read;
  char c_read;
  COMPARE(std::fread(&c_read, sizeof(char), 1, file.get()), 1u);
  read_binary(npart_read, file);
  // std::cout << c_read << std::endl;
  // std::cout << npart_read << " " << npart << std::endl;
  return (c_read == 'p') && (npart_read == npart);
}

/* function to check we reached the end of the file */
static bool check_end_of_file(const FilePtr &file) {
  return std::feof(file.get()) == 0;
}

TEST(fullhistory_format) {
  /* create two smashon particles */
  Particles particles;
  const ParticleData p1 = particles.insert(Test::smashon_random());
  const ParticleData p2 = particles.insert(Test::smashon_random());

  /* Create elastic interaction (smashon + smashon). */
  const double impact_parameter = 1.473;
  const bool empty_event = false;
  EventInfo event = Test::default_event_info(impact_parameter, empty_event);
  ScatterActionPtr action = std::make_unique<ScatterAction>(p1, p2, 0.);
  action->add_all_scatterings(Test::default_finder_parameters());
  action->generate_final_state();
  ParticleList final_particles = action->outgoing_particles();
  const double rho = 0.123;

  const std::filesystem::path collisionsoutputfilepath =
      testoutputpath / "collisions_oscar2013.bin";
  std::filesystem::path collisionsoutputfilepath_unfinished =
      collisionsoutputfilepath;
  collisionsoutputfilepath_unfinished += ".unfinished";
  {
    /* Set the most verbose option */
    OutputParameters output_par = OutputParameters();
    output_par.coll_printstartend = true;
    output_par.coll_extended = false;
    output_par.quantities["Collisions"] = {};
    /* Create an instance of binary output */
    auto bin_output = create_binary_output("Oscar2013_bin", "Collisions",
                                           testoutputpath, output_par);
    VERIFY(std::filesystem::exists(collisionsoutputfilepath_unfinished));

    /* Write initial state output: the two smashons we created */
    bin_output->at_eventstart(particles, event_id, event);
    bin_output->at_interaction(*action, rho);

    /* Final state output */
    action->perform(&particles, 1);
    bin_output->at_eventend(particles, event_id, event);
  }
  VERIFY(!std::filesystem::exists(collisionsoutputfilepath_unfinished));
  VERIFY(std::filesystem::exists(collisionsoutputfilepath));

  /*
   * Now we have an artificially generated binary output.
   * Let us try if we can read and understand it.
   */

  {
    FilePtr binF = fopen(collisionsoutputfilepath.native(), "rb");
    VERIFY(binF.get());
    // Header
    std::vector<char> buf(4);
    std::string magic, smash_version;
    uint16_t format_version_number, extended_format;

    COMPARE(std::fread(&buf[0], 1, 4, binF.get()), 4u);  // magic number
    magic.assign(&buf[0], 4);
    read_binary(format_version_number, binF);  // format version number
    read_binary(extended_format, binF);        // whether extended format
    read_binary(smash_version, binF);          // smash version

    COMPARE(magic, "SMSH");
    COMPARE(format_version_number, current_format_version);
    COMPARE(extended_format, 0);
    COMPARE(smash_version, SMASH_VERSION);

    // particles at event start: expect two smashons
    VERIFY(compare_particles_block_header(event_id, 2, binF));
    VERIFY(compare_particle(p1, binF));
    VERIFY(compare_particle(p2, binF));

    // interaction: 2 smashons -> 2 smashons
    VERIFY(compare_interaction_block_header(2, 2, *action, rho, binF));
    VERIFY(compare_particle(p1, binF));
    VERIFY(compare_particle(p2, binF));
    VERIFY(compare_particle(final_particles[0], binF));
    VERIFY(compare_particle(final_particles[1], binF));

    // paricles at event end: two smashons
    VERIFY(compare_particles_block_header(event_id, 2, binF));
    VERIFY(compare_particle(final_particles[0], binF));
    VERIFY(compare_particle(final_particles[1], binF));

    // event end line
    VERIFY(compare_final_block_header(event_id, impact_parameter, empty_event,
                                      binF));
    VERIFY(check_end_of_file(binF));
  }

  VERIFY(std::filesystem::remove(collisionsoutputfilepath));
}

TEST(particles_format) {
  /* create two smashon particles */
  const auto particles =
      Test::create_particles(2, [] { return Test::smashon_random(); });
  const double impact_parameter = 4.382;
  const bool empty_event = false;
  EventInfo event = Test::default_event_info(impact_parameter, empty_event);
  const ParticleList initial_particles = particles->copy_to_vector();

  const std::filesystem::path particleoutputpath =
      testoutputpath / "particles_oscar2013.bin";
  std::filesystem::path particleoutputpath_unfinished = particleoutputpath;
  particleoutputpath_unfinished += ".unfinished";
  {
    /* Set the most verbose option */
    OutputParameters output_par = OutputParameters();
    output_par.part_extended = false;
    output_par.part_only_final = OutputOnlyFinal::No;

    output_par.quantities["Particles"] = {};
    /* Create an instance of binary output */
    auto bin_output = create_binary_output("Oscar2013_bin", "Particles",
                                           testoutputpath, output_par);
    VERIFY(bool(bin_output));
    VERIFY(std::filesystem::exists(particleoutputpath_unfinished));

    /* Write initial state output: the two smashons we created */
    bin_output->at_eventstart(*particles, event_id, event);
    /* Interaction smashon + smashon -> smashon */
    ParticleList final_state = {Test::smashon_random()};
    particles->replace(initial_particles, final_state);

    DensityParameters dens_par(Test::default_parameters());
    bin_output->at_intermediate_time(*particles, nullptr, dens_par, event_id,
                                     event);

    /* Final state output */
    bin_output->at_eventend(*particles, event_id, event);
  }
  const ParticleList final_particles = particles->copy_to_vector();
  VERIFY(!std::filesystem::exists(particleoutputpath_unfinished));
  VERIFY(std::filesystem::exists(particleoutputpath));

  /*
   * Now we have an artificially generated binary output.
   * Let us try if we can read and understand it.
   */

  {
    FilePtr binF = fopen(particleoutputpath.native(), "rb");
    VERIFY(binF.get());
    // Header
    std::vector<char> buf(4);
    std::string magic, smash_version;
    uint16_t format_version_number, extended_format;

    COMPARE(std::fread(&buf[0], 1, 4, binF.get()), 4u);  // magic number
    magic.assign(&buf[0], 4);
    read_binary(format_version_number, binF);  // format version number
    read_binary(extended_format, binF);        // whether extended format
    read_binary(smash_version, binF);          // smash version

    COMPARE(magic, "SMSH");
    COMPARE(format_version_number, current_format_version);
    COMPARE(extended_format, 0);
    COMPARE(smash_version, SMASH_VERSION);

    int32_t npart;
    // particles at event start: expect two smashons
    npart = 2;
    VERIFY(compare_particles_block_header(event_id, npart, binF));
    VERIFY(compare_particle(initial_particles[0], binF));
    VERIFY(compare_particle(initial_particles[1], binF));

    // Periodic output: already after interaction. One smashon expected.
    npart = 1;
    VERIFY(compare_particles_block_header(event_id, npart, binF));
    VERIFY(compare_particle(final_particles[0], binF));

    // particles at event end
    npart = 1;
    VERIFY(compare_particles_block_header(event_id, npart, binF));
    VERIFY(compare_particle(final_particles[0], binF));

    // after end of event
    VERIFY(compare_final_block_header(event_id, impact_parameter, empty_event,
                                      binF));
    VERIFY(check_end_of_file(binF));
  }

  VERIFY(std::filesystem::remove(particleoutputpath));
}

TEST(extended) {
  /* create two smashon particles */
  Particles particles;
  const ParticleData p1 = particles.insert(Test::smashon_random());
  const ParticleData p2 = particles.insert(Test::smashon_random());

  /* Create elastic interaction (smashon + smashon). */
  ScatterActionPtr action = std::make_unique<ScatterAction>(p1, p2, 0.);
  action->add_all_scatterings(Test::default_finder_parameters());
  action->generate_final_state();
  ParticleList final_particles = action->outgoing_particles();
  const double rho = 0.123;

  const double impact_parameter = 1.473;
  const bool empty_event = true;
  EventInfo event = Test::default_event_info(impact_parameter, empty_event);

  const std::filesystem::path collisionsoutputfilepath =
      testoutputpath / "collisions_oscar2013_extended.bin";
  std::filesystem::path collisionsoutputfilepath_unfinished =
      collisionsoutputfilepath;
  collisionsoutputfilepath_unfinished += ".unfinished";
  {
    OutputParameters output_par = OutputParameters();
    output_par.coll_printstartend = true;
    output_par.coll_extended = true;
    output_par.quantities["Collisions"] = {};
    /* Create an instance of binary output */
    auto bin_output = create_binary_output("Oscar2013_bin", "Collisions",
                                           testoutputpath, output_par);
    VERIFY(std::filesystem::exists(collisionsoutputfilepath_unfinished));

    /* Write initial state output: the two smashons we created */
    bin_output->at_eventstart(particles, event_id, event);
    bin_output->at_interaction(*action, rho);

    /* Final state output */
    action->perform(&particles, 1);
    bin_output->at_eventend(particles, event_id, event);
  }
  VERIFY(!std::filesystem::exists(collisionsoutputfilepath_unfinished));
  VERIFY(std::filesystem::exists(collisionsoutputfilepath));

  /*
   * Now we have an artificially generated binary output.
   * Let us try if we can read and understand it.
   */

  {
    FilePtr binF = fopen(collisionsoutputfilepath.native(), "rb");
    VERIFY(binF.get());
    // Header
    std::vector<char> buf(4);
    std::string magic, smash_version;
    uint16_t format_version_number, extended_version;

    COMPARE(std::fread(&buf[0], 1, 4, binF.get()), 4u);  // magic number
    magic.assign(&buf[0], 4);
    read_binary(format_version_number, binF);  // format version number
    read_binary(extended_version, binF);       // whether extended format
    read_binary(smash_version, binF);          // smash version

    COMPARE(magic, "SMSH");
    COMPARE(static_cast<int>(format_version_number), current_format_version);
    COMPARE(extended_version, 1);
    COMPARE(smash_version, SMASH_VERSION);

    // particles at event atart: expect two smashons
    VERIFY(compare_particles_block_header(event_id, 2, binF));
    compare_particle_extended(p1, binF);
    compare_particle_extended(p2, binF);

    // interaction: 2 smashons -> 2 smashons
    VERIFY(compare_interaction_block_header(2, 2, *action, rho, binF));
    compare_particle_extended(p1, binF);
    compare_particle_extended(p2, binF);
    compare_particle_extended(final_particles[0], binF);
    compare_particle_extended(final_particles[1], binF);

    // paricles at event end: two smashons
    VERIFY(compare_particles_block_header(event_id, 2, binF));
    for (const auto &particle : particles) {
      compare_particle_extended(particle, binF);
    }

    // event end line
    VERIFY(compare_final_block_header(event_id, impact_parameter, empty_event,
                                      binF));
    VERIFY(check_end_of_file(binF));
  }

  VERIFY(std::filesystem::remove(collisionsoutputfilepath));
}

TEST(initial_conditions_format) {
  // Create 1 particle
  Particles particles;
  ParticleData p1 = particles.insert(Test::smashon_random());
  p1.set_4position(FourVector(2.3, 1.35722, 1.42223, 1.5));  // tau = 1.74356

  // Create and perform action ("hypersurface crossing")
  ActionPtr action = std::make_unique<FluidizationAction>(p1, p1, 0.0);
  action->generate_final_state();
  action->perform(&particles, 1);

  const bool empty_event = false;
  const double impact_parameter = 0.0;
  EventInfo event = Test::default_event_info(impact_parameter, empty_event);

  const std::filesystem::path particleoutputpath =
      testoutputpath / "SMASH_IC_oscar2013.bin";
  std::filesystem::path particleoutputpath_unfinished = particleoutputpath;
  particleoutputpath_unfinished += ".unfinished";

  {
    OutputParameters output_par = OutputParameters();
    output_par.part_extended = false;
    double density = 0.0;
    /* Create an instance of binary output */
    auto bin_output = create_binary_output(
        "Oscar2013_bin", "Initial_Conditions", testoutputpath, output_par);
    VERIFY(bool(bin_output));
    VERIFY(std::filesystem::exists(particleoutputpath_unfinished));

    /* Write event start information: This should do nothing for IC output */
    bin_output->at_interaction(*action, density);

    /* Write particle line for hypersurface crossing */
    bin_output->at_interaction(*action, density);

    /* Event end output */
    bin_output->at_eventend(particles, event_id, event);
  }
  VERIFY(!std::filesystem::exists(particleoutputpath_unfinished));
  VERIFY(std::filesystem::exists(particleoutputpath));

  /* Read the afore created output */
  {
    FilePtr binF = fopen(particleoutputpath.native(), "rb");
    VERIFY(binF.get());
    // Header
    std::vector<char> buf(4);
    std::string magic, smash_version;
    uint16_t format_version_number, extended_format;

    COMPARE(std::fread(&buf[0], 1, 4, binF.get()), 4u);  // magic number
    magic.assign(&buf[0], 4);
    read_binary(format_version_number, binF);  // format version number
    read_binary(extended_format, binF);        // whether extended format
    read_binary(smash_version, binF);          // smash version

    COMPARE(magic, "SMSH");
    COMPARE(static_cast<int>(format_version_number), current_format_version);
    COMPARE(extended_format, 0);
    COMPARE(smash_version, SMASH_VERSION);

    int32_t npart = 1;  // expect one particle in output

    VERIFY(compare_initial_conditions_interaction_block_header(npart, binF));
    VERIFY(compare_particle(p1, binF));

    VERIFY(check_end_of_file(binF));
  }
  VERIFY(std::filesystem::remove(particleoutputpath));
}

/* Function to read and compare particle with custom quantities */
static bool compare_particle_custom(const ParticleData &p,
                                    const FilePtr &file) {
  int charge, strangeness;
  FourVector pos;
  read_binary(pos, file);
  read_binary(charge, file);
  read_binary(strangeness, file);

  return (pos == p.position()) && (charge == p.type().charge()) &&
         (strangeness == p.type().strangeness());
}

TEST(custom) {
  /* create two smashon particles */
  const auto particles =
      Test::create_particles(2, [] { return Test::smashon_random(); });
  const double impact_parameter = 4.382;
  const bool empty_event = false;
  EventInfo event = Test::default_event_info(impact_parameter, empty_event);
  const ParticleList initial_particles = particles->copy_to_vector();

  const std::filesystem::path particleoutputpath =
      testoutputpath / "particles_custom.bin";
  std::filesystem::path particleoutputpath_unfinished = particleoutputpath;
  particleoutputpath_unfinished += ".unfinished";
  {
    /* Set the most verbose option */
    OutputParameters output_par = OutputParameters();
    output_par.part_extended = false;
    output_par.part_only_final = OutputOnlyFinal::No;
    output_par.quantities["Particles"] = {"t", "x",      "y",
                                          "z", "charge", "strangeness"};

    /* Create an instance of binary output */
    auto bin_output =
        create_binary_output("Binary", "Particles", testoutputpath, output_par);
    VERIFY(bool(bin_output));
    VERIFY(std::filesystem::exists(particleoutputpath_unfinished));

    /* Write initial state output: the two smashons we created */
    bin_output->at_eventstart(*particles, event_id, event);
    /* Interaction smashon + smashon -> smashon */
    ParticleList final_state = {Test::smashon_random()};
    particles->replace(initial_particles, final_state);

    DensityParameters dens_par(Test::default_parameters());
    bin_output->at_intermediate_time(*particles, nullptr, dens_par, event_id,
                                     event);

    /* Final state output */
    bin_output->at_eventend(*particles, event_id, event);
  }
  const ParticleList final_particles = particles->copy_to_vector();
  VERIFY(!std::filesystem::exists(particleoutputpath_unfinished));
  VERIFY(std::filesystem::exists(particleoutputpath));

  /*
   * Now we have an artificially generated binary output.
   * Let us try if we can read and understand it.
   */

  {
    FilePtr binF = fopen(particleoutputpath.native(), "rb");
    VERIFY(binF.get());
    // Header
    std::vector<char> buf(4);
    std::string magic, smash_version;
    uint16_t format_version_number, extended_version;

    COMPARE(std::fread(&buf[0], 1, 4, binF.get()), 4u);  // magic number
    magic.assign(&buf[0], 4);
    VERIFY(std::fread(&format_version_number, sizeof(format_version_number), 1,
                      binF.get()) == 1);
    VERIFY(std::fread(&extended_version, sizeof(extended_version), 1,
                      binF.get()) == 1);
    read_binary(smash_version, binF);  // smash version

    COMPARE(magic, "SMSH");
    COMPARE(static_cast<int>(format_version_number), current_format_version);
    COMPARE(extended_version, 2);
    COMPARE(smash_version, SMASH_VERSION);

    int npart;
    // particles at event start: expect two smashons
    npart = 2;
    VERIFY(compare_particles_block_header(event_id, npart, binF));
    VERIFY(compare_particle_custom(initial_particles[0], binF));
    VERIFY(compare_particle_custom(initial_particles[1], binF));

    // Periodic output: already after interaction. One smashon expected.
    npart = 1;
    VERIFY(compare_particles_block_header(event_id, npart, binF));
    VERIFY(compare_particle_custom(final_particles[0], binF));

    // particles at event end
    npart = 1;
    VERIFY(compare_particles_block_header(event_id, npart, binF));
    VERIFY(compare_particle_custom(final_particles[0], binF));

    // after end of event
    VERIFY(compare_final_block_header(event_id, impact_parameter, empty_event,
                                      binF));
    VERIFY(check_end_of_file(binF));
  }

  VERIFY(std::filesystem::remove(particleoutputpath));
}
