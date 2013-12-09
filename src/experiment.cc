/*
 *
 *    Copyright (c) 2012-2013
 *      Hannah Petersen <petersen@fias.uni-frankfurt.de>
 *      
 *    GNU General Public License (GPLv3)
 *
 */


#include "include/Experiment.h"
#include "include/BoundaryConditions.h"
#include "include/Parameters.h"
#include "include/param-reader.h"

std::unique_ptr<Experiment> Experiment::create(const std::int &modus)
{
  typedef std::unique_ptr<Experiment> ExperimentPointer;
  if (modus == 1) {
    return ExperimentPointer{new ExperimentImplementation<BoxBoundaryConditions>};
  } else if (modus == 2) {
    return ExperimentPointer{new ExperimentImplementation<SphereBoundaryConditions>};
  } else {
    throw std::string("Invalid boundaryCondition requested from Experiment::create.");
  }
}

template <typename BoundaryConditions>
void ExperimentImplementation<BoundaryConditions>::config()
{
    process_config(Parameters, char *path);
    new BoundaryConditions bc;
    bc->assign_params_general(Parameters);
    bc->assign_params_specific(Parameters);
    
}

