/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_EXPERIMENT_H_
#define SRC_INCLUDE_EXPERIMENT_H_

#include <memory>
#include <list>
#include "../include/Parameters.h"
#include "../include/BoundaryConditions.h"
#include "../include/Particles.h"
#include "../include/CrossSections.h"

class Experiment
{
public:
    static std::unique_ptr<Experiment> create(char *modus);
    virtual void config(std::list<Parameters> configuration) = 0;
    virtual void initialize(char *path)=0;
};

template <typename BoundaryConditions> class ExperimentImplementation : public Experiment
{
public:
    virtual void config(std::list<Parameters> configuration);
    virtual void initialize(char *path);
private:
    BoundaryConditions bc;
    Particles *particles = new Particles;
    CrossSections *cross_sections = new CrossSections;
};

#endif  // SRC_INCLUDE_Experiment_H_

