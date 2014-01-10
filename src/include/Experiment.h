/*
 *    Copyright (c) 2013-2014
 *      Hannah Petersen <petersen@fias.uni-frankfurt.de>
 *      
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_EXPERIMENT_H_
#define SRC_INCLUDE_EXPERIMENT_H_

#include <memory>
#include <list>

#include "../include/CrossSections.h"
#include "../include/Modus.h"
#include "../include/Parameters.h"
#include "../include/Particles.h"

class Experiment
{
public:
    static std::unique_ptr<Experiment> create(char *modus_chooser);
    virtual void configure(std::list<Parameters> configuration) = 0;
    virtual void initialize(char *path)=0;
    virtual void run_time_evolution()=0;
    virtual void end()=0;
};

template <typename Modus> class ExperimentImplementation : public Experiment
{
public:
    virtual void configure(std::list<Parameters> configuration);
    virtual void initialize(char *path);
    virtual void run_time_evolution();
    virtual void end();
    
private:
    Modus bc;
    Particles *particles = new Particles;
    CrossSections *cross_sections = new CrossSections;
};

#endif  // SRC_INCLUDE_Experiment_H_

