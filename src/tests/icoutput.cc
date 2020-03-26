/*
 *
 *    Copyright (c) 2019-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <vir/test.h>  // This include has to be first

#include "setup.h"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "../include/smash/icoutput.h"
#include "../include/smash/outputinterface.h"

using namespace smash;

static const bf::path testoutputpath = bf::absolute(SMASH_TEST_OUTPUT_PATH);

TEST(directory_is_created) {
  bf::create_directories(testoutputpath);
  VERIFY(bf::exists(testoutputpath));
}

TEST(init_particletypes) { Test::create_smashon_particletypes(); }

TEST(particlelist_format) {
  // Create 1 particle
  Particles particles;
  ParticleData p1 = particles.insert(Test::smashon_random());
  /*
  We need a little trick here to make sure the particle is actually written to
  the output. By construction, the ASCII IC output does not contain spectator
  particles. This is triggered by whether or not the particle has prior
  interactions (collisions_per_particle in the particle HistoryData). The test
  particle p1 has no prior interactions, so we manually have to change it's
  history. As collisions_per_particle is a private member of the HistoryData
  class, we can only change this property by setting the entire history. For
  this we also need a ParticleList usually containing the mother particles.
  Herein we are not interested in any actions so we set the mother particle
  to be the particle itself, only modifying the number of prior interactions.
  Physics-wise this has no meaning or interpretation.
  */
  // Create particle list with mother particles
  ParticleList mother_list = {ParticleData{p1.type()}};
  // Manually enforce that number of collisions = 1 (and therefore != 0)
  p1.set_history(1, 1, ProcessType::None, 0.01, mother_list);
  p1.set_4position(FourVector(2.3, 1.35722, 1.42223, 1.5));  // tau = 1.74356

  // Create and perform action ("hypersurface crossing")
  ActionPtr action = make_unique<HypersurfacecrossingAction>(p1, p1, 0.0);
  action->generate_final_state();
  action->perform(&particles, 1);

  const int event_id = 0;
  const bool empty_event = false;
  const double impact_parameter = 0.0;

  const bf::path outputfilepath = testoutputpath / "SMASH_IC.dat";
  bf::path outputfilepath_unfinished = outputfilepath;
  outputfilepath_unfinished += ".unfinished";

  {
    OutputParameters out_par = OutputParameters();
    auto IC_output =
        make_unique<ICOutput>(testoutputpath, "Initial_Conditions", out_par);

    VERIFY(bool(IC_output));
    VERIFY(bf::exists(outputfilepath_unfinished));

    /* Initial state output (write event number) */
    IC_output->at_eventstart(particles, event_id);

    /* Interaction Output (write particle's coordinates on hypersurface) */
    IC_output->at_interaction(*action, 0.);

    /* Final state output (event end line) */
    IC_output->at_eventend(particles, event_id, impact_parameter, empty_event);
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
          "# " VERSION_MAJOR
          " initial conditions: hypersurface of constant proper time\n"
          "# tau x y eta mt px py Rap pdg charge\n"
          "# fm fm fm none GeV GeV GeV none none e\n"
          "# event 1 start\n";
      int line_number = 0;
      do {
        line_number++;
        std::getline(outputfile, line);
        output_header += line + '\n';
      } while (line_number < 4);  // we expect the header to have 4 lines
      COMPARE(output_header, header);

      /* Check particle data */
      outputfile >> item;
      // Compare tau
      COMPARE_ABSOLUTE_ERROR(std::stod(item), p1.position().tau(), 1e-6);
      outputfile >> item;  // jump over x
      outputfile >> item;  // and also y
      outputfile >> item;
      // Compare eta
      COMPARE_ABSOLUTE_ERROR(std::stod(item), p1.position().eta(), 1e-6);
    }
    VERIFY(bf::remove(outputfilepath));
  }
}
