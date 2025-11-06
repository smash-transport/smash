/*
 *
 *    Copyright (c) 2019-2020,2022-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/icoutput.h"

#include <filesystem>

#include "setup.h"
#include "smash/fluidizationaction.h"
#include "smash/outputinterface.h"

using namespace smash;

static const std::filesystem::path testoutputpath =
    std::filesystem::absolute(SMASH_TEST_OUTPUT_PATH);

TEST(directory_is_created) {
  std::filesystem::create_directories(testoutputpath);
  VERIFY(std::filesystem::exists(testoutputpath));
}

TEST(init_particletypes) { Test::create_smashon_particletypes(); }

TEST(particlelist_format) {
  // Create 1 particle
  Particles particles;
  ParticleData p1 = particles.insert(Test::smashon_random());
  /*
  We need a little trick here to make sure the particle is actually written to
  the output. By construction, the IC output does not contain spectator
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
  p1.set_history(1, 0, ProcessType::None, 0.01, mother_list);
  p1.set_4position(FourVector(2.3, 1.35722, 1.42223, 1.5));  // tau = 1.74356

  // Create and perform action ("hypersurface crossing")
  ActionPtr action = std::make_unique<FluidizationAction>(p1, p1, 0.0);
  action->generate_final_state();
  action->perform(&particles, 1);

  const EventLabel event_id = {0, 0};
  const bool empty_event = false;
  const double impact_parameter = 0.0;
  EventInfo event = Test::default_event_info(impact_parameter, empty_event);

  const std::filesystem::path outputfilepath =
      testoutputpath / "SMASH_IC_For_vHLLE.dat";
  std::filesystem::path outputfilepath_unfinished = outputfilepath;
  outputfilepath_unfinished += ".unfinished";

  {
    OutputParameters out_par = OutputParameters();
    auto IC_output = std::make_unique<ICOutput>(testoutputpath,
                                                "Initial_Conditions", out_par);

    VERIFY(bool(IC_output));
    VERIFY(std::filesystem::exists(outputfilepath_unfinished));

    /* Initial state output (write event number) */
    IC_output->at_eventstart(particles, event_id, event);

    /* Interaction Output (write particle's coordinates on hypersurface) */
    IC_output->at_interaction(*action, 0.);

    /* Final state output (event end line) */
    IC_output->at_eventend(particles, event_id, event);
  }

  VERIFY(!std::filesystem::exists(outputfilepath_unfinished));
  VERIFY(std::filesystem::exists(outputfilepath));

  {
    std::fstream outputfile;
    outputfile.open(outputfilepath, std::ios_base::in);
    if (outputfile.good()) {
      std::string line, item;

      /* Check header */
      std::string output_header = "";
      std::string header =
          "# " SMASH_VERSION
          " initial conditions: hypersurface of constant proper time\n"
          "# tau x y eta mt px py Rap pdg charge "
          "baryon_number strangeness\n"
          "# fm fm fm none GeV GeV GeV none none e "
          "none none\n"
          "# event 0 ensemble 0 start\n";
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
      COMPARE_ABSOLUTE_ERROR(std::stod(item), p1.hyperbolic_time(), 1e-6);
      outputfile >> item;  // jump over x
      outputfile >> item;  // and also y
      outputfile >> item;
      // Compare eta
      COMPARE_ABSOLUTE_ERROR(std::stod(item), p1.spatial_rapidity(), 1e-6);
    }
    VERIFY(std::filesystem::remove(outputfilepath));
  }
}
