/*
 *
 *    Copyright (c) 2012-2013
 *      Hannah Petersen <petersen@fias.uni-frankfurt.de>
 *      
 *    GNU General Public License (GPLv3)
 *
 */


#include "Experiment.h"
#include "BoundaryConditions.h"
#include "Parameters.h"
#include "param-reader.h"

std::unique_ptr<Experiment> Experiment::create(const std::int &BoundaryConditions)
{
  typedef std::unique_ptr<Experiment> ExperimentPointer;
  if (boundaryConditions == 1) {
    return ExperimentPointer{new ExperimentImplementation<BoxBoundaryConditions>};
  } else if (boundaryConditions == 2) {
    return ExperimentPointer{new ExperimentImplementation<SphereBoundaryConditions>};
  } else {
    throw std::string("Invalid boundaryCondition requested from Experiment::create.");
  }
}

template <typename Modus>
void ExperimentImplementation<Modus>::config()
{
    process_config(Parameters, char *path);
    new BoundaryConditions bc(Parameters);
    assign_params_general;
    assign_params_specific;
    
}

