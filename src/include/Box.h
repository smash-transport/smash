/*
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_BOX_H_
#define SRC_INCLUDE_BOX_H_

#include <stdint.h>
#include <time.h>
#include <cmath>
#include <list>

#include "../include/CrossSections.h"
#include "../include/BoundaryConditions.h"
#include "../include/Particles.h"
#include "../include/time.h"
#include "../include/Parameters.h"


class BoxBoundaryConditions : public BoundaryConditions
{
public:
  /* default constructor with probable values */
    BoxBoundaryConditions(): initial_condition(1), length(10.0f), temperature(0.1f),energy_initial(0.0f), number_density_initial(0.0f),
      time_start(set_timer_start()) {}
    /* special class funtions */
    virtual int evolve(Particles *particles, CrossSections *cross_sections);
    virtual void assign_params_specific(std::list<Parameters> *configuration);
    virtual void initial_conditions(Particles *particles);
    virtual void check_collision_geometry(Particles *particles, CrossSections *cross_sections, std::list<int> *collision_list, BoxBoundaryConditions const &box, size_t *rejection_conflict);
    inline timespec set_timer_start();
  private:
    /* initial condition */
    int initial_condition;
    /* Cube edge length */
    float length;
    /* Temperature of the Boltzmann distribution for thermal initialization */
    float temperature;
    /* initial total energy of the box */
    float energy_initial;
    /* initial number density of the box */
    float number_density_initial;
    /* starting time of the simulation */
    timespec time_start;
};

    /* set the timer to the actual time in nanoseconds precision */
    timespec inline BoxBoundaryConditions::set_timer_start(void) {
        timespec time;
        clock_gettime(&time);
        return time;
    }



#endif  // SRC_INCLUDE_BOX_H_
