/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_TESTS_SETUP_H_
#define SRC_TESTS_SETUP_H_

#include "../include/cxx14compat.h"
#include "../include/experiment.h"
#include "../include/particles.h"
#include "../include/particletype.h"
#include "../include/particledata.h"
#include "../include/random.h"

#include <boost/filesystem.hpp>

namespace Smash {
namespace Test {

inline void create_actual_particletypes() {
#include <particles.txt.h>
  ParticleType::create_type_list(data);
}

static constexpr float smashon_mass = 0.123f;
static constexpr float smashon_width = 1.2f;

inline void create_smashon_particletypes() {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
      "smashon " +
      std::to_string(smashon_mass) + ' ' + std::to_string(smashon_width) +
      " 661\n");
}

auto random_value = Random::make_uniform_distribution(-15.0, +15.0);

inline ParticleData smashon(const FourVector &position = {
                                random_value(), random_value(), random_value(),
                                random_value()},
                            const FourVector &momentum = {
                                smashon_mass, random_value(), random_value(),
                                random_value()},
                            int id = -1) {
  ParticleData p{ParticleType::find(0x661), id};
  p.set_4position(position);
  p.set_4momentum(momentum);
  return p;
}

/**
 * Return a configuration object filled with data from src/config.yaml. Note
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
  Configuration c{TEST_CONFIG_PATH};
  if (!overrides.empty()) {
    c.merge_yaml(overrides);
  }
  return c;
}

inline std::unique_ptr<ExperimentBase> experiment(
    const Configuration &c = configuration()) {
  return ExperimentBase::create(c);
}

inline std::unique_ptr<ExperimentBase> experiment(const char *configOverrides) {
  return ExperimentBase::create(configuration(configOverrides));
}

using ParticlesPtr = std::unique_ptr<Particles>;
template <typename G>
inline ParticlesPtr create_particles(int n, G &&generator) {
  ParticlesPtr p = make_unique<Particles>();
  for (auto i = n; i; --i) {
    p->add_data(generator());
  }
  return p;
}

inline ParticlesPtr create_particles(
    const std::initializer_list<ParticleData> &init) {
  ParticlesPtr p = make_unique<Particles>();
  for (const auto &data : init) {
    p->add_data(data);
  }
  return p;
}

}  // namespace Test
}  // namespace Smash

#endif  // SRC_TESTS_SETUP_H_
