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
#include "include/parameters.h"

namespace Smash {

class BoxModus;
class Configuration;
struct ExperimentParameters;

class BoxModus : public ModusDefault {
 public:
  explicit BoxModus(Configuration modus_config);

  void print_startup();  // TODO(mkretz): needs to be discoverable from an
                         // outside "printer"

  void initial_conditions(Particles *particles,
                          const ExperimentParameters &parameters);

  int sanity_check(Particles *particles);

  void check_collision_geometry(Particles *particles,
                                CrossSections *cross_sections,
                                std::list<int> *collision_list,
                                size_t *rejection_conflict,
                                const ExperimentParameters &parameters);

  void propagate(Particles *particles, const ExperimentParameters &parameters);

 private:
  /* initial condition */
  int initial_condition_ = 1;
  /* Cube edge length */
  float length_ = 10.f;
  /* Temperature of the Boltzmann distribution for thermal initialization */
  float temperature_ = 0.1f;
  /* initial number density of the box */
  float number_density_initial_ = 0.f;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_BOXMODUS_H_
