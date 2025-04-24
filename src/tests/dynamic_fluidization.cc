/*
 *
 *    Copyright (c) 2023-2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "setup.h"
#include "smash/dynamicfluidfinder.h"
#include "smash/experiment.h"
#include "smash/random.h"

using namespace smash;
using smash::Test::Momentum;
using smash::Test::Position;

static const FourVector unit{1, 0, 0, 0};
static const DensityParameters dens_par(Test::default_parameters());
static const InitialConditionParameters ic_par(
    Test::default_dynamic_IC_parameters());

static ParticleData create_particle_at(PdgCode pdg, Position pos) {
  ParticleData particle{ParticleType::find(pdg)};
  particle.set_4position(pos);
  particle.set_4momentum(Momentum{particle.pole_mass(), 0, 0, 0});
  return particle;
}

// Functional tests for the fluidization finder
TEST(fluidization_finder) {
  Test::create_actual_particletypes();
  Test::create_actual_decaymodes();
  std::vector<Particles> ensembles(1);
  Particles &particles = ensembles[0];
  ParticleData a =
      particles.insert(create_particle_at(0x2212, Position{0, -5, -5, -5}));
  ParticleData b =
      particles.insert(create_particle_at(0x111, Position{0, 5, 5, 5}));

  a.set_4momentum(unit);  // set rest mass to 1 GeV
  particles.update_particle(a, a);

  // create large and coarse lattice
  const std::array<double, 3> length{20, 20, 20}, origin{-10, -10, -10};
  const std::array<int, 3> cell_array{2, 2, 2};
  auto lat = std::make_unique<RectangularLattice<EnergyMomentumTensor>>(
      length, cell_array, origin, false, LatticeUpdate::EveryTimestep);
  update_lattice_accumulating_ensembles(lat.get(), LatticeUpdate::EveryTimestep,
                                        DensityType::Hadron, dens_par,
                                        ensembles, false);

  // create fake background
  auto background = std::make_unique<std::map<int32_t, double>>();
  background->emplace(a.id(), 1.0);
  background->emplace(b.id(), 0.0);

  // create finder with 1 GeV threshold
  DynamicFluidizationFinder finder(*lat.get(), *background.get(), ic_par);

  // check that background works
  VERIFY(finder.above_threshold(a));
  VERIFY(!finder.above_threshold(b));

  EnergyMomentumTensor Tmunu;
  lat->value_at(a.position().threevec(), Tmunu);
  // with zero 3-momentum the LRF is the lab frame
  COMPARE(unit, Tmunu.landau_frame_4velocity());
  // for unit mass/momentum, T^{00} evaluates to the normalization
  COMPARE(Tmunu[0], dens_par.norm_factor_sf());

  // initial particles are not fluidizable and there is no action
  VERIFY(!finder.is_process_fluidizable(a.get_history()));
  ActionList actions =
      finder.find_actions_in_cell(particles.copy_to_vector(), 1, 0, {});
  COMPARE(actions.size(), 0u);

  ParticleList mother{ParticleData(ParticleType::find(0x2124), -1)};
  a.set_history(1, a.id_process(), ProcessType::Decay, 0, mother);
  particles.update_particle(a, a);
  b.set_history(1, b.id_process(), ProcessType::Decay, 0, mother);
  particles.update_particle(b, b);

  // decays are fluidizable by default
  VERIFY(finder.is_process_fluidizable(a.get_history()));
  actions = finder.find_actions_in_cell(particles.copy_to_vector(), 1, 0, {});
  COMPARE(actions.size(), 1u);

  // action found only after formation time
  a.set_formation_time(0.9);
  particles.update_particle(a, a);
  actions = finder.find_actions_in_cell(particles.copy_to_vector(), 0.8, 0, {});
  COMPARE(actions.size(), 0u);
  actions = finder.find_actions_in_cell(particles.copy_to_vector(), 1, 0, {});
  COMPARE(actions.size(), 1u);

  // modifying background has an effect even after the finder is constructed
  background->at(b.id()) = 1.0;
  VERIFY(finder.above_threshold(b));
  actions = finder.find_actions_in_cell(particles.copy_to_vector(), 1, 0, {});
  COMPARE(actions.size(), 2u);

  // the finder is ignored after max time
  const double max_time = ic_par.max_time.value();
  a.set_4position(Position{max_time + 1, 5, 5, 5});
  particles.update_particle(a, a);
  b.set_4position(Position{max_time + 1, -5, -5, -5});
  particles.update_particle(b, b);
  actions =
      finder.find_actions_in_cell(particles.copy_to_vector(), max_time, 0, {});
  COMPARE(actions.size(), 0u);
}

/*
 * The logic here is that a small region with 50 protons
 * (from appropriate processes) will fluidize completely
 */
TEST(dense_region) {
  ParticleList mother{ParticleData(ParticleType::find(0x2124), -1)};
  const int large_particle_number = 50;
  const double range = 0.5;
  std::vector<Particles> ensembles(1);
  Particles &particles = ensembles[0];
  random::set_seed(random::generate_63bit_seed());
  for (int i = 0; i < large_particle_number; i++) {
    auto [x, y, z] = std::make_tuple<double, double, double>(
        random::uniform(-range, range), random::uniform(-range, range),
        random::uniform(-range, range));
    ensembles[0].insert(create_particle_at(0x2212, Position{0, x, y, z}));
  }
  // create small and fine lattice
  const std::array<double, 3> length{2, 2, 2}, origin{-1, -1, -1};
  const std::array<int, 3> cell_array{20, 20, 20};
  auto lat = std::make_unique<RectangularLattice<EnergyMomentumTensor>>(
      length, cell_array, origin, false, LatticeUpdate::EveryTimestep);

  update_lattice_accumulating_ensembles(lat.get(), LatticeUpdate::EveryTimestep,
                                        DensityType::Hadron, dens_par,
                                        ensembles, false);
  auto background = std::make_unique<std::map<int32_t, double>>();
  DynamicFluidizationFinder finder(*lat.get(), *background.get(), ic_par);

  for (ParticleData &p : particles) {
    background->emplace(p.id(), 0.0);
    p.set_history(1, p.id_process(), ProcessType::Decay, 0, mother);
  }

  ActionList actions =
      finder.find_actions_in_cell(particles.copy_to_vector(), 1, 0, {});
  COMPARE(actions.size(), large_particle_number);
}
