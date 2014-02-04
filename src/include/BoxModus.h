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

#include "../include/CrossSections.h"
#include "../include/ModusDefault.h"
#include "../include/Particles.h"
#include "../include/Parameters.h"

class BoxModus;
class ExperimentParameters;

class BoxModus : public ModusDefault {
 public:
    BoxModus() = default;

    /* special class funtions */
    // XXX: -> ctor
    void assign_params(std::list<Parameters> *configuration);
    // XXX: needs to be discoverable from an outside "printer"
    void print_startup();

    void initial_conditions(Particles *particles,
                            const ExperimentParameters &parameters);

    int sanity_check(Particles *particles);

    void check_collision_geometry(Particles *particles,
                                  CrossSections *cross_sections,
                                  std::list<int> *collision_list,
                                  size_t *rejection_conflict,
                                  const ExperimentParameters &parameters);

    void propagate(Particles *particles,
                   const ExperimentParameters &parameters);

 private:
    FourVector boundary_condition(FourVector position,
                                  bool *boundary_hit);

    /* initial condition */
    int initial_condition_ = 1;
    /* Cube edge length */
    float length_ = 10.f;
    /* Temperature of the Boltzmann distribution for thermal initialization */
    float temperature_ = 0.1f;
    /* initial number density of the box */
    float number_density_initial_ = 0.f;
};

#endif  // SRC_INCLUDE_BOXMODUS_H_
