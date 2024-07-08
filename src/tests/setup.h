/*
 *
 *    Copyright (c) 2015,2017-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_TESTS_SETUP_H_
#define SRC_TESTS_SETUP_H_

#include <filesystem>

#include "smash/decaymodes.h"
#include "smash/experiment.h"
#include "smash/outputinterface.h"
#include "smash/particledata.h"
#include "smash/particles.h"
#include "smash/particletype.h"
#include "smash/random.h"
#include "smash/scatteractionsfinderparameters.h"

namespace smash {
namespace Test {

/**
 * \addtogroup unittest
 * @{
 */

/**
 * Creates the ParticleType list containing the actual particles that SMASH
 * uses.
 */
inline void create_actual_particletypes() {
#ifndef DOXYGEN
/// not visible to doxygen, but compiled
#include <particles.txt.h>
#endif
  ParticleType::create_type_list(data);
}

/**
 * Creates the DecayModes list containing the actual decay modes that SMASH
 * uses.
 */
inline void create_actual_decaymodes() {
#ifndef DOXYGEN
/// not visible to doxygen, but compiled
#include <decaymodes.txt.h>
#endif
  DecayModes::load_decaymodes(data);
}

/// The mass of the smashon particle.
static constexpr double smashon_mass = 0.123;
/// The decay width of the smashon particle.
static constexpr double smashon_width = 1.2;
/// The PDG code of the smashon particle.
static constexpr const char smashon_pdg_string[] = "661";

/**
 * Creates a ParticleType list containing only the smashon test particle.
 */
inline void create_smashon_particletypes() {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
      "σ " +
      std::to_string(smashon_mass) + " " + std::to_string(smashon_width) +
      " + 661\n");
}

/**
 * Creates a ParticleType list containing only the smashon test particle with
 * width 0 (stable).
 */
inline void create_stable_smashon_particletypes() {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
      "σ " +
      std::to_string(smashon_mass) + " 0.0 + 661\n");
}

/// A FourVector that is marked as a position vector.
struct Position : public FourVector {
  using FourVector::FourVector;
};
/// A FourVector that is marked as a momentum vector.
struct Momentum : public FourVector {
  using FourVector::FourVector;
};

/**
 * Create a particle with 0 position and momentum vectors and optionally a given
 * \p id.
 */
inline ParticleData smashon(int id = -1) {
  ParticleData p{ParticleType::find(0x661), id};
  return p;
}
/**
 * Create a particle with 0 momentum vector, the given \p position, and
 * optionally a given \p id.
 */
inline ParticleData smashon(const Position &position, int id = -1) {
  ParticleData p{ParticleType::find(0x661), id};
  p.set_4position(position);
  p.set_formation_time(position[0]);
  return p;
}
/**
 * Create a particle with 0 position vector, the given \p momentum, and
 * optionally a given \p id.
 */
inline ParticleData smashon(const Momentum &momentum, int id = -1) {
  ParticleData p{ParticleType::find(0x661), id};
  p.set_4momentum(momentum);
  return p;
}
/**
 * Create a particle with the given \p position and \p momentum vectors, and
 * optionally a given \p id.
 */
inline ParticleData smashon(const Position &position, const Momentum &momentum,
                            int id = -1) {
  ParticleData p{ParticleType::find(0x661), id};
  p.set_4position(position);
  p.set_4momentum(momentum);
  p.set_formation_time(position[0]);
  p.set_unpolarized_spin_vector();
  return p;
}
/**
 * Create a particle with the given \p position and \p momentum vectors, and
 * optionally a given \p id.
 *
 * Convenience overload of the above to allow arbitrary order of momentum and
 * position.
 */
inline ParticleData smashon(const Momentum &momentum, const Position &position,
                            int id = -1) {
  ParticleData p{ParticleType::find(0x661), id};
  p.set_4position(position);
  p.set_4momentum(momentum);
  p.set_formation_time(position[0]);
  p.set_unpolarized_spin_vector();
  return p;
}
/**
 * Create a particle with random position and momentum vectors and optionally a
 * given \p id.
 */
inline ParticleData smashon_random(int id = -1) {
  auto random_value = random::make_uniform_distribution(-15.0, +15.0);
  ParticleData p{ParticleType::find(0x661), id};
  auto random_time = random_value();
  p.set_formation_time(random_time);
  p.set_4position(
      {random_time, random_value(), random_value(), random_value()});
  p.set_4momentum(smashon_mass,
                  {random_value(), random_value(), random_value()});
  return p;
}

/**
 * Create an experiment given an input configuration.
 */
inline std::unique_ptr<ExperimentBase> experiment(Configuration c) {
  return ExperimentBase::create(c, ".");
}

/**
 * Generate a list of particles from the given generator function.
 */
template <typename G>
inline ParticleList create_particle_list(std::size_t n, G &&generator) {
  ParticleList list;
  list.reserve(n);
  for (auto i = n; i; --i) {
    list.emplace_back(generator());
  }
  return list;
}

/// A type alias for a unique_ptr of Particles.
using ParticlesPtr = std::unique_ptr<Particles>;

/**
 * Creates a Particles object and fills it with \p n particles generated by the
 * \p generator function.
 */
template <typename G>
inline ParticlesPtr create_particles(int n, G &&generator) {
  ParticlesPtr p = std::make_unique<Particles>();
  for (auto i = n; i; --i) {
    p->insert(generator());
  }
  return p;
}

/**
 * Creates a Particles object and fills it with the particles passed as
 * initializer_list to this function.
 */
inline ParticlesPtr create_particles(
    const std::initializer_list<ParticleData> &init) {
  ParticlesPtr p = std::make_unique<Particles>();
  for (const auto &data : init) {
    p->insert(data);
  }
  return p;
}

/// returns BitSet of 2->2 reactions, where everything is on
inline ReactionsBitSet all_reactions_included() {
  return ReactionsBitSet().set();
}

/// returns BitSet for multi-particle reactions, where everything is off
inline MultiParticleReactionsBitSet no_multiparticle_reactions() {
  return MultiParticleReactionsBitSet().reset();
}

/**
 * Creates a standard ExperimentParameters object which works for almost all
 * testing purposes.
 *
 * If needed you can set the testparticles parameter to a different value than
 * 1.
 */
inline ExperimentParameters default_parameters(
    int testparticles = 1, double dt = 0.1,
    CollisionCriterion criterion = CollisionCriterion::Geometric,
    bool strings = false,
    NNbarTreatment nnbar_treatment = NNbarTreatment::NoAnnihilation,
    ReactionsBitSet included_2to2 = all_reactions_included()) {
  return ExperimentParameters{
      std::make_unique<UniformClock>(0., dt, 300.0),  // labclock
      std::make_unique<UniformClock>(0., 1., 300.0),  // outputclock
      1,                                              // ensembles
      testparticles,                                  // testparticles
      DerivativesMode::CovariantGaussian,             // derivatives mode
      RestFrameDensityDerivativesMode::Off,  // rest frame derivatives mode
      FieldDerivativesMode::ChainRule,       // field derivatives mode
      SmearingMode::CovariantGaussian,       // smearing mode
      1.0,                                   // Gaussian smearing width
      4.0,                                   // Gaussian smearing cut-off
      0.333333,                              // discrete smearing weight
      2.0,                                   // triangular smearing range
      criterion,                             // collision criterion
      true,                                  // two_to_one
      included_2to2,
      no_multiparticle_reactions(),
      strings,
      1.0,
      nnbar_treatment,
      0.,           // low energy sigma_NN cut-off
      false,        // potential_affect_threshold
      -1.0,         // box_length
      200.0,        // max. cross section
      2.5,          // fixed min. cell length
      1.0,          // cross section scaling
      false,        // in thermodynamics outputs spectators are included
      false,        // do non-strong decays
      true,         // decay initial particles
      std::nullopt  // use monash tune, not known
  };
}

/**
 * Creates a standard ScatterActionsFinderParameters object which works for
 * almost all testing purposes.
 *
 * The selected arguments are changed between different tests, which requires
 * setting the key by hand. This is not directly possible for enums, so one must
 * do it case by case.
 */
inline ScatterActionsFinderParameters default_finder_parameters(
    double elastic_parameter = 10,
    NNbarTreatment nnbar_treatment = NNbarTreatment::NoAnnihilation,
    ReactionsBitSet included_2to2 = all_reactions_included(),
    bool strings_switch = true, bool use_AQM = false,
    bool strings_with_probability = false,
    TotalCrossSectionStrategy xs_strategy =
        TotalCrossSectionStrategy::BottomUp) {
  Configuration config{
      R"(
  Collision_Term:
    Only_Warn_For_High_Probability: true
    Pseudoresonance: None
  )"};
  config.set_value(InputKeys::collTerm_elasticCrossSection, elastic_parameter);
  config.set_value(InputKeys::collTerm_useAQM, use_AQM);
  config.set_value(InputKeys::collTerm_stringsWithProbability,
                   strings_with_probability);
  if (xs_strategy == TotalCrossSectionStrategy::BottomUp) {
    config.merge_yaml(InputKeys::collTerm_totXsStrategy.as_yaml("BottomUp"));
  } else if (xs_strategy == TotalCrossSectionStrategy::TopDown) {
    config.merge_yaml(InputKeys::collTerm_totXsStrategy.as_yaml("TopDown"));
  } else if (xs_strategy == TotalCrossSectionStrategy::TopDownMeasured) {
    config.merge_yaml(
        InputKeys::collTerm_totXsStrategy.as_yaml("TopDownMeasured"));
  }
  return ScatterActionsFinderParameters(
      config,
      default_parameters(1, 0.1, CollisionCriterion::Geometric, strings_switch,
                         nnbar_treatment, included_2to2));
}

/// Creates default EventInfo object for testing purposes
inline EventInfo default_event_info(double impact_parameter = 0.0,
                                    bool empty_event = false) {
  return EventInfo{impact_parameter, 0.0,  0.0, 0.0, 0.0, 0.0, 1, 1,
                   empty_event,      false};
}

/// Creates a default StringProcessInterface object for testing
inline std::unique_ptr<StringProcess> default_string_process_interface() {
  return std::make_unique<StringProcess>(
      1.0,      // String_Tension
      1.0,      // String_Formation_Time
      0.5,      // Gluon_Beta
      0.001,    // Gluon_Pmin
      2.0,      // Quark_Alpha
      7.0,      // Quark_Beta
      0.16,     // Strange_Supp
      0.036,    // Diquark_Supp
      0.42,     // Sigma_Perp
      0.2,      // StringZ_A_Leading
      2.0,      // StringZ_B_Leading
      2.0,      // StringZ_A
      0.55,     // StringZ_B
      0.5,      // String_Sigma_T
      1.0,      // Form_Time_Factor
      false,    // Mass_Dependent_Formation_Times
      1. / 3.,  // Prob_proton_to_d_uu
      true,     // Separate_Fragment_Baryon
      0.15,     // Popcorn_Rate
      false);   // Use_Monash_Tune
}

/// Creates default parameters for dynamic IC
inline InitialConditionParameters default_dynamic_IC_parameters() {
  InitialConditionParameters parameters{};
  parameters.type = FluidizationType::Dynamic;
  parameters.fluidizable_processes = FluidizableProcessesBitSet{}.set();
  parameters.energy_density_threshold = 0.5;
  parameters.min_time = 0;
  parameters.max_time = 100;
  parameters.num_fluid_cells = 50;
  parameters.formation_time_fraction = 1;
  parameters.smearing_kernel_at_0 = std::pow(2 * M_PI, -1.5);
  parameters.delay_initial_elastic = false;
  return parameters;
}

/**
 * @}
 */
}  // namespace Test
}  // namespace smash

#endif  // SRC_TESTS_SETUP_H_
