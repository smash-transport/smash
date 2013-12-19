/*
 *
 *    Copyright (c) 2012-2013
 *      Hannah Petersen <petersen@fias.uni-frankfurt.de>
 *      
 *    GNU General Public License (GPLv3)
 *
 */

#include <list>

#include "include/Experiment.h"
#include "include/BoundaryConditions.h"
#include "include/Parameters.h"
#include "include/param-reader.h"
#include "include/Box.h"
#include "include/outputroutines.h"
#include "include/input-particles.h"
#include "include/input-decaymodes.h"
#include "include/CrossSections.h"

//#include "include/Sphere.h"

std::unique_ptr<Experiment> Experiment::create(char *modus)
{
  typedef std::unique_ptr<Experiment> ExperimentPointer;
  if (strcmp(modus,"Box") == 0) {
    return ExperimentPointer{new ExperimentImplementation<BoxBoundaryConditions>};
  }
//  else if (modus == 2) {
//    return ExperimentPointer{new ExperimentImplementation<SphereBoundaryConditions>};
//  }
    else {
    throw std::string("Invalid boundaryCondition requested from Experiment::create.");
  }
}



template <typename BoundaryConditions>
void ExperimentImplementation<BoundaryConditions>::config(std::list<Parameters> configuration)
{
    bc.assign_params(&configuration);
    warn_wrong_params(&configuration);
    bc.print_startup();
/* reducing cross section according to number of test particle */
    if (bc.testparticles > 1) {
      printf("IC test particle: %i\n", bc.testparticles);
      bc.cross_section = bc.cross_section / bc.testparticles;
      printf("Elastic cross section: %g mb\n", bc.cross_section);
    }
}

template <typename BoundaryConditions>
void ExperimentImplementation<BoundaryConditions>::initialize(char *path)
{
    srand48(bc.seed);
    input_particles(particles, path);
    input_decaymodes(particles, path);
    cross_sections->add_elastic_parameter(bc.cross_section);
    //
//      bc.initial_conditions;
//    bc.print_initial();
    
}

