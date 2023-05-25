#include <filesystem>
#include <iostream>
#include <memory>
#include <string>

#include "smash/action.h"
#include "smash/collidermodus.h"
#include "smash/config.h"
#include "smash/experiment.h"
#include "smash/forwarddeclarations.h"
#include "smash/library.h"
#include "smash/particles.h"

void update_fluidization_background(
    smash::Experiment<smash::ColliderModus> &exp);

int main() {
  try {
    std::cout << "\nTest-run SMASH\n--------------" << '\n';

    // All the input that is needed
    const std::string config_file(SMASH_TOP_LEVEL_DIR
                                  "/build/config_dynIC.yaml");
    const std::filesystem::path output_path("./data");
    const std::string tabulations_path("./tabulations");
    const std::string particles_file(SMASH_TOP_LEVEL_DIR
                                     "/input/particles.txt");
    const std::string decaymodes_file(SMASH_TOP_LEVEL_DIR
                                      "/input/decaymodes.txt");

    ////////////////////////////
    // Setup SMASH            //
    ////////////////////////////

    // 1) Set-up config
    auto config = smash::setup_config_and_logging(config_file, particles_file,
                                                  decaymodes_file);

    // 2) Do addtional configurations e.g. set a custom end time by
    float new_end_time = 180.0;
    config.set_value({"General", "End_Time"}, new_end_time);

    std::string smash_version = SMASH_VERSION;

    // 3) Intialize decaymodes, particletypes, tabulations
    smash::initialize_particles_decays_and_tabulations(config, smash_version,
                                                       tabulations_path);

    // Create experiment
    auto experiment =
        smash::Experiment<smash::ColliderModus>(config, output_path);

    // Run manually in timesteps
    experiment.initialize_new_event();

    const double delta_time = 1.0; // time step for adding the background, we should make it =min(dt_vhlle,dt_smash)
    const double start_time = 0.0;
    double current_time = start_time;

    do {
      current_time += delta_time;
      std::cout << "Propagate to " << current_time << '\n';
      experiment.run_time_evolution(current_time);
      // the iterator for smash::Particles is not defined here somehow, so I
      // convert to vector
      smash::ParticleList particles =
          experiment.first_ensemble()->copy_to_vector();
      std::map<int32_t, double> fake_background;
      for (const smash::ParticleData &p : particles) {
       // get positions of SMASH particles
       smash::ThreeVector pos = p.position().threevec(); 
       // get e_den at each position from vhlle
       // add to corresponding particle index
       fake_background.insert_or_assign(p.id(), 0.2);
      }
      experiment.update_fluidization_background(fake_background);

    } while (current_time < new_end_time);  // maximum fluid/fluidization time

    std::cout << "Last time: " << current_time << '\n';
    experiment.do_final_decays();
    experiment.final_output();

  } catch (std::exception &e) {
    std::cout << "SMASH failed with the following error:\n" << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return 0;
}
