/*
 *
 *    Copyright (c) 2019-2020,2022-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "setup.h"
#include "smash/experiment.h"
#include "smash/fluidizationaction.h"
#include "smash/hypersurfacecrossingfinder.h"

using namespace smash;
using smash::Test::Momentum;
using smash::Test::Position;

TEST(init_particle_types) {
  Test::create_actual_particletypes();
  Test::create_actual_decaymodes();
}

TEST(hypersurface_crossing_action) {
  // Create 2 particles
  Particles particles;

  // Supposed to not cross hypersurface in next time step
  ParticleData a{ParticleType::find(0x2212)};
  a.set_4position(Position{0.5, 0.1, 0., 0.45});
  a.set_4momentum(Momentum{0.943, 0., 0., 0.1});

  // Supposed to cross hypersurface in next time step
  ParticleData b{ParticleType::find(0x2212)};
  b.set_4position(Position{0.5, 0., 0.1, 0.1});
  b.set_4momentum(Momentum{1.386, -1.0, 0., 0.2});

  particles.insert(a);
  particles.insert(b);

  ParticleList part_list = particles.copy_to_vector();

  // create finder at tau = 0.5fm without rapidity or pT cut
  double proper_time = 0.5;
  HyperSurfaceCrossActionsFinder finder(proper_time, 0.0, 0.0);
  FluidizationAction::remove_particle_ = true;

  // no grid means no grid cell volume
  const double grid_cell_vol = 0.0;

  // We do not test with frozen Fermi motion, so beam_mom vector is empty
  const std::vector<FourVector> beam_mom = {};

  // Find actions
  constexpr double time_step = 0.1;
  ActionList actions = finder.find_actions_in_cell(part_list, time_step,
                                                   grid_cell_vol, beam_mom);

  // Action list should only contain one element since one particle a crosses
  // the hypersurface in the given time step
  // Implicit test of HyperSurfaceCrossActionsFinder::crosses_hypersurface
  COMPARE(actions.size(), 1u);

  for (auto &action : actions) {
    // perform action
    uint32_t id_process = 1;

    // propagate particle to the action's time of execution
    const double t_until_action =
        action->time_of_execution() - b.position().x0();
    const ThreeVector &v = b.velocity();
    const FourVector distance = FourVector(0.0, v * t_until_action);
    FourVector position = b.position() + distance;
    position.set_x0(action->time_of_execution());
    b.set_4position(position);

    COMPARE_ABSOLUTE_ERROR(b.position().tau(), proper_time, 1e-7);

    // Perform action
    action->generate_final_state();
    action->perform(&particles, id_process);
    std::cout << action->get_type() << " " << ProcessType::Fluidization;
    COMPARE(action->get_type(), ProcessType::Fluidization);
    // 1 incoming, no outgoing particles expected (particle removed)
    COMPARE(action->incoming_particles().size(), 1u);
    COMPARE(action->outgoing_particles().size(), 0u);
  }
}
