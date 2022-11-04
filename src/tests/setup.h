/*
 *
 *    Copyright (c) 2015,2017-2022
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
  return p;
}
/**
 * Create a particle with random position and momentum vectors and optionally a
 * given \p id.
 */
inline ParticleData smashon_random(int id = -1) {
  auto random_value = random::make_uniform_distribution(-15.0, +15.0);
  ParticleData p{ParticleType::find(0x661), id};
  p.set_4position(
      {random_value(), random_value(), random_value(), random_value()});
  p.set_4momentum(smashon_mass,
                  {random_value(), random_value(), random_value()});
  return p;
}

/**
 * Return a configuration object filled with data from input/config.yaml. Note
 * that a change to that file may affect test results if you use it.
 *
 * If you want specific values in the config for testing simply overwrite the
 * relevant settings e.g. with:
 * \code
 * auto config = Test::configuration(
 *   "General:\n"
 *   "  Modus: Box\n"
 *   "  Testparticles: 100\n"
 * );
 * \endcode
 */
inline Configuration configuration(std::string overrides = {}) {
  Configuration c{std::filesystem::path{TEST_CONFIG_PATH} / "input"};
  if (!overrides.empty()) {
    c.merge_yaml(overrides);
  }
  return c;
}

/**
 * Create an experiment.
 *
 * If you want a specific configuration you can pass it as parameter, otherwise
 * it will use the result from configuration above.
 */
inline std::unique_ptr<ExperimentBase> experiment(
    Configuration &&c = configuration()) {
  return ExperimentBase::create(c, ".");
}

/**
 * Creates an experiment using the default config and the specified \p
 * configOverrides.
 */
inline std::unique_ptr<ExperimentBase> experiment(const char *configOverrides) {
  Configuration c = configuration(configOverrides);
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
    CollisionCriterion crit = CollisionCriterion::Geometric) {
  return ExperimentParameters{
      std::make_unique<UniformClock>(0., dt),  // labclock
      std::make_unique<UniformClock>(0., 1.),  // outputclock
      1,                                       // ensembles
      testparticles,                           // testparticles
      DerivativesMode::CovariantGaussian,      // derivatives mode
      RestFrameDensityDerivativesMode::Off,    // rest frame derivatives mode
      FieldDerivativesMode::ChainRule,         // field derivatives mode
      SmearingMode::CovariantGaussian,         // smearing mode
      1.0,                                     // Gaussian smearing width
      4.0,                                     // Gaussian smearing cut-off
      0.333333,                                // discrete smearing weight
      2.0,                                     // triangular smearing range
      crit,
      true,  // two_to_one
      all_reactions_included(),
      no_multiparticle_reactions(),
      false,  // strings switch
      1.0,
      NNbarTreatment::NoAnnihilation,
      0.,     // low energy sigma_NN cut-off
      false,  // potential_affect_threshold
      -1.0,   // box_length
      200.0,  // max. cross section
      2.5,    // fixed min. cell length
      1.0,    // cross section scaling
      false,  // in thermodynamics outputs spectators are included
      false   // do weak decays
  };
}

/**
 * Creates a standard ScatterActionsFinderParameters object which works for
 * almost all testing purposes.
 *
 * The selected arguments are changed between different tests.
 */
inline ScatterActionsFinderParameters default_finder_parameters(
    double elastic_parameter = 10,
    NNbarTreatment nnbar_treatment = NNbarTreatment::NoAnnihilation,
    ReactionsBitSet included_2to2 = all_reactions_included(),
    bool strings_switch = true, bool use_AQM = false,
    bool strings_with_probability = false) {
  return {
      elastic_parameter,
      0.,    // low_snn_cut
      1.,    // scale_xs
      0.,    // additional_el_xs
      200.,  // maximum_cross_section
      CollisionCriterion::Geometric,
      nnbar_treatment,
      included_2to2,
      no_multiparticle_reactions(),
      1,      // testparticles
      true,   // two_to_one
      false,  // allow_first_collisions_within_nucleus
      strings_switch,
      use_AQM,
      strings_with_probability,
      true  // only_warn_for_high_prob
  };
}

/// Creates default EventInfo object for testing purposes
inline EventInfo default_event_info(double impact_parameter = 0.0,
                                    bool empty_event = false) {
  return EventInfo{impact_parameter, 0.0,  0.0, 0.0, 0.0, 0.0, 1, 1,
                   empty_event,      false};
}

/**
 * @}
 */
}  // namespace Test
}  // namespace smash

#endif  // SRC_TESTS_SETUP_H_
