/*
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_BOXMODUS_H_
#define SRC_INCLUDE_BOXMODUS_H_

#include <stdint.h>
#include <cmath>
#include <list>

#include "include/crosssections.h"
#include "include/modusdefault.h"
#include "include/particles.h"

namespace Smash {

class BoxModus;
class Configuration;
struct ExperimentParameters;

class BoxModus : public ModusDefault {
 public:
  explicit BoxModus(Configuration modus_config,
           const ExperimentParameters &parameters);

  void print_startup();  // TODO(mkretz): needs to be discoverable from an
                         // outside "printer"

  float initial_conditions(Particles *particles,
                          const ExperimentParameters &parameters);

  int sanity_check(Particles *particles);

  void propagate(Particles *particles, const ExperimentParameters &parameters, const OutputsList &);

 private:
  /* initial condition */
  const int initial_condition_;
  /* Cube edge length */
  const float length_;
  /* Temperature of the Boltzmann distribution for thermal initialization */
  const float temperature_;
  /* initial number density of the box */
  float number_density_initial_ = 0.f;
  /// initial time of the box
  const float start_time_ = 0.0f;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_BOXMODUS_H_
