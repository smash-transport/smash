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
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../include/binaryoutputcollisions.h"
#include "../include/binaryoutputparticles.h"
#include "../include/clock.h"
#include "../include/outputinterface.h"
#include "../include/processbranch.h"
#include "../include/scatteraction.h"

using namespace Smash;

static const bf::path testoutputpath = bf::absolute(SMASH_TEST_OUTPUT_PATH);

TEST(directory_is_created) {
  bf::create_directories(testoutputpath);
  VERIFY(bf::exists(testoutputpath));
}

TEST(init_particletypes) { Test::create_smashon_particletypes(); }

static const int current_format_version = 5;

/* A set of convenient functions to read binary */

static void read_binary(std::string &s, FILE *file) {
  std::int32_t size = s.size();
  COMPARE(std::fread(&size, sizeof(std::int32_t), 1, file), 1u);
  std::vector<char> buf(size);
  COMPARE(std::fread(&buf[0], 1, size, file), static_cast<size_t>(size));
  s.assign(&buf[0], size);
}

static void read_binary(FourVector &v, FILE *file) {
  COMPARE(std::fread(v.begin(), sizeof(*v.begin()), 4, file), 4u);
}

static void read_binary(std::int32_t &x, FILE *file) {
  COMPARE(std::fread(&x, sizeof(x), 1, file), 1u);
}

static void read_binary(double &x, FILE *file) {
  COMPARE(std::fread(&x, sizeof(x), 1, file), 1u);
}

/* Function to read and compare particle */
static bool compare_particle(const ParticleData &p, FILE *file) {
  int id, pdgcode, charge;
  double mass;
  FourVector pos, mom;
  read_binary(pos, file);
  read_binary(mass, file);
  read_binary(mom, file);
  read_binary(pdgcode, file);
  read_binary(id, file);
  read_binary(charge, file);
  // std::cout << p.id() << " " << id << std::endl;
  // std::cout << p.pdgcode().get_decimal() << " " << pdgcode << std::endl;
  return (p.id() == id) && (p.pdgcode().get_decimal() == pdgcode) &&
         (pos == p.position()) && (mom == p.momentum());
}

/* function to read and compare particle block header */
static bool compare_particles_block_header(const int &npart, FILE *file) {
  int npart_read;
  char c_read;
  COMPARE(std::fread(&c_read, sizeof(char), 1, file), 1u);
  read_binary(npart_read, file);
  // std::cout << c_read << std::endl;
  // std::cout << npart_read << " " << npart << std::endl;
  return (c_read == 'p') && (npart_read == npart);
}

/* function to read and compare collision block header */
static bool compare_interaction_block_header(const int &nin, const int &nout,
                                             const Action &action, double rho,
                                             FILE *file) {
  int nin_read, nout_read, process_type_read;
  double rho_read, weight_read, partial_weight_read;
  char c_read;
  int process_type = static_cast<int>(action.get_type());
  COMPARE(std::fread(&c_read, sizeof(char), 1, file), 1u);
  read_binary(nin_read, file);
  read_binary(nout_read, file);
  COMPARE(std::fread(&rho_read, sizeof(double), 1, file), 1u);
  COMPARE(std::fread(&weight_read, sizeof(double), 1, file), 1u);
  COMPARE(std::fread(&partial_weight_read, sizeof(double), 1, file), 1u);
  read_binary(process_type_read, file);
  // std::cout << c_read << std::endl;
  // std::cout << nin_read << " " << nin << std::endl;
  // std::cout << nout_read << " " << nout << std::endl;
  // std::cout << rho << std::endl;
  return (c_read == 'i') && (nin_read == nin) && (nout_read == nout) &&
         (rho_read == rho) && (weight_read == action.raw_weight_value()) &&
         (partial_weight_read == action.partial_weight()) &&
         (process_type_read == process_type);
}

/* function to read and compare event end line */
static bool compare_final_block_header(const int &ev,
                                       const double &impact_parameter,
                                       FILE *file) {
  int ev_read;
  char c_read;
  double b_read;
  COMPARE(std::fread(&c_read, sizeof(char), 1, file), 1u);
  read_binary(ev_read, file);
  COMPARE(std::fread(&b_read, sizeof(double), 1, file), 1u);
  return (c_read == 'f') && (ev_read == ev) && (b_read == impact_parameter);
}

TEST(fullhistory_format) {
  /* Set the most verbose option */
  OutputParameters output_par = OutputParameters();
  output_par.coll_printstartend = true;
  output_par.coll_extended = false;

  /* Create an instance of binary output */
  std::unique_ptr<BinaryOutputCollisions> bin_output =
      make_unique<BinaryOutputCollisions>(testoutputpath,
        "Collisions", output_par);
  const bf::path collisionsoutputfilepath =
      testoutputpath / "collisions_binary.bin";
  VERIFY(bf::exists(collisionsoutputfilepath));

  /* create two smashon particles */
  Particles particles;
  const ParticleData p1 = particles.insert(Test::smashon_random());
  const ParticleData p2 = particles.insert(Test::smashon_random());

  int event_id = 0;
  /* Write initial state output: the two smashons we created */
  bin_output->at_eventstart(particles, event_id);

  /* Create elastic interaction (smashon + smashon). */
  ScatterActionPtr action = make_unique<ScatterAction>(p1, p2, 0.);
  action->add_all_processes(10., true, true, 0., true,
                            NNbarTreatment::NoAnnihilation);
  action->generate_final_state();
  ParticleList final_particles = action->outgoing_particles();
  const double rho = 0.123;
  bin_output->at_interaction(*action, rho);

  /* Final state output */
  action->perform(&particles, 1);
  const double impact_parameter = 1.473;
  bin_output->at_eventend(particles, event_id, impact_parameter);

  /*
   * Now we have an artificially generated binary output.
   * Let us try if we can read and understand it.
   */

  // Open file as a binary
  FILE *binF;
  const auto filename = collisionsoutputfilepath.native();
  binF = fopen(filename.c_str(), "rb");
  VERIFY(binF);
  // Header
  std::vector<char> buf(4);
  std::string magic, smash_version;
  int format_version_number;

  COMPARE(std::fread(&buf[0], 1, 4, binF), 4u);  // magic number
  magic.assign(&buf[0], 4);
  read_binary(format_version_number, binF);  // format version number
  read_binary(smash_version, binF);          // smash version

  COMPARE(magic, "SMSH");
  COMPARE(format_version_number, current_format_version);
  COMPARE(smash_version, VERSION_MAJOR);

  // particles at event atart: expect two smashons
  VERIFY(compare_particles_block_header(2, binF));
  VERIFY(compare_particle(p1, binF));
  VERIFY(compare_particle(p2, binF));

  // interaction: 2 smashons -> 2 smashons
  VERIFY(compare_interaction_block_header(2, 2, *action, rho, binF));
  VERIFY(compare_particle(p1, binF));
  VERIFY(compare_particle(p2, binF));
  VERIFY(compare_particle(final_particles[0], binF));
  VERIFY(compare_particle(final_particles[1], binF));

  // paricles at event end: two smashons
  VERIFY(compare_particles_block_header(2, binF));
  VERIFY(compare_particle(final_particles[0], binF));
  VERIFY(compare_particle(final_particles[1], binF));

  // event end line
  VERIFY(compare_final_block_header(event_id, impact_parameter, binF));

  // remove file
  VERIFY(!std::fclose(binF));
  VERIFY(!std::remove(filename.c_str()));
}

TEST(particles_format) {
  /* Set the most verbose option */
  OutputParameters output_par = OutputParameters();
  output_par.part_extended = false;
  output_par.part_only_final = false;

  /* Create an instance of binary output */
  std::unique_ptr<BinaryOutputParticles> bin_output =
      make_unique<BinaryOutputParticles>(testoutputpath,
          "Particles", output_par);
  VERIFY(bf::exists(testoutputpath / "particles_binary.bin"));

  /* create two smashon particles */
  const auto particles =
      Test::create_particles(2, [] { return Test::smashon_random(); });

  int event_id = 0;
  /* Write initial state output: the two smashons we created */
  bin_output->at_eventstart(*particles, event_id);

  /* Interaction smashon + smashon -> smashon */
  ParticleList initial_particles = particles->copy_to_vector();
  ParticleList final_state = {Test::smashon_random()};
  particles->replace(initial_particles, final_state);
  ParticleList final_particles = particles->copy_to_vector();
  Clock clock;

  DensityParameters dens_par(Test::default_parameters());
  bin_output->at_intermediate_time(*particles, clock, dens_par);

  /* Final state output */
  const double impact_parameter = 4.382;
  bin_output->at_eventend(*particles, event_id, impact_parameter);
  /*
   * Now we have an artificially generated binary output.
   * Let us try if we can read and understand it.
   */

  // Open file as a binary
  FILE *binF;
  const auto filename = (testoutputpath / "particles_binary.bin").native();
  binF = fopen(filename.c_str(), "rb");
  VERIFY(binF);
  // Header
  std::vector<char> buf(4);
  std::string magic, smash_version;
  int format_version_number;

  COMPARE(std::fread(&buf[0], 1, 4, binF), 4u);  // magic number
  magic.assign(&buf[0], 4);
  read_binary(format_version_number, binF);  // format version number
  read_binary(smash_version, binF);          // smash version

  COMPARE(magic, "SMSH");
  COMPARE(format_version_number, current_format_version);
  COMPARE(smash_version, VERSION_MAJOR);

  int npart;
  // particles at event start: expect two smashons
  npart = 2;
  VERIFY(compare_particles_block_header(npart, binF));
  VERIFY(compare_particle(initial_particles[0], binF));
  VERIFY(compare_particle(initial_particles[1], binF));

  // Periodic output: already after interaction. One smashon expected.
  npart = 1;
  VERIFY(compare_particles_block_header(npart, binF));
  VERIFY(compare_particle(final_particles[0], binF));

  // paricles at event end: nothing expected, because only_final option is off

  // after end of event
  VERIFY(compare_final_block_header(event_id, impact_parameter, binF));

  // remove file
  VERIFY(!std::fclose(binF));
  VERIFY(!std::remove(filename.c_str()));
}
