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

#include "../include/experiment.h"
#include "../include/particles.h"
#include "../include/particletype.h"
#include "../include/particledata.h"

#include <boost/filesystem.hpp>

namespace Smash {
namespace Test {

inline void create_actual_particletypes() {
#include <particles.txt.h>
  ParticleType::create_type_list(data);
}

inline void create_smashon_particletypes() {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
      "smashon 0.123 1.2 661\n");
}

inline ParticleData smashon(const FourVector &position,
                            const FourVector &momentum, int id = -1) {
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

}  // namespace Test
}  // namespace Smash

#endif  // SRC_TESTS_SETUP_H_
