/*
Example of how other models can interface with SMASH, as for example done in
JETSCAPE.

Main steps are
- Setup smash configuration object
- Initialize particles, decay modes and tabulation for this config
- Create experiment
- Perform run time evolution

For more details have a look in the SMASH source for the main function in
smash.cc, Experiment::run() and the functions library.h.
*/

#include <filesystem>
#include <iostream>
#include <memory>
#include <string>

#include "smash/action.h"
#include "smash/collidermodus.h"
#include "smash/config.h"
#include "smash/experiment.h"
#include "smash/forwarddeclarations.h"
#include "smash/input_keys.h"
#include "smash/library.h"

int main() {
  try {
    std::cout << "\nTest-run SMASH\n--------------" << '\n';

    // All the input that is needed
    const std::string config_file(SMASH_INPUT_DIR "/config.yaml");
    const std::filesystem::path output_path("./data");
    const std::string tabulations_path("./tabulations");
    const std::string particles_file(SMASH_INPUT_DIR "/particles.txt");
    const std::string decaymodes_file(SMASH_INPUT_DIR "/decaymodes.txt");

    // Ensure output_path exists
    std::filesystem::create_directories(output_path);

    ////////////////////////////
    // Setup SMASH            //
    ////////////////////////////

    // 1) Set-up config
    auto config = smash::setup_config_and_logging(config_file, particles_file,
                                                  decaymodes_file);

    // 2) Do additional configurations e.g. set a custom end time by
    double new_end_time = 180.0;
    config.set_value(smash::InputKeys::gen_endTime, new_end_time);
    // ...

    std::string smash_version = SMASH_VERSION;

    // 3) Initialize decay modes, particle types
    const auto hash = smash::initialize_particles_decays_and_return_hash(
        config, smash_version);

    // Create experiment
    auto experiment =
        smash::Experiment<smash::ColliderModus>(config, output_path);

    // Tabulate integrals
    smash::tabulate_resonance_integrals(hash, tabulations_path);

    ////////////////////////////////////////////////////////////////////////////
    // Run the experiment in a special way here, mimicking new JETSCAPE features
    // i.e. feed additional particles at different times and perform the runtime
    // evolution manually.
    // If you want to just run SMASH as usual just call  experiment.run()
    // instead.
    ////////////////////////////////////////////////////////////////////////////

    // Run manually in time steps
    experiment.initialize_new_event();

    const double delta_time = 20.0;
    const double start_time = 0.0;
    double current_time = start_time;

    do {
      current_time += delta_time;
      std::cout << "Propagate to " << current_time << '\n';
      experiment.run_time_evolution(current_time);
    } while (current_time < new_end_time);

    std::cout << "Last time: " << current_time << '\n';
    experiment.do_final_decays();
    experiment.final_output();

  } catch (std::exception &e) {
    std::cout << "SMASH failed with the following error:\n" << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return 0;
}
