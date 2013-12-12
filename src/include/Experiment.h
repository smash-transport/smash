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
#include "../include/Parameters.h"

class Experiment
{
public:
    static std::unique_ptr<Experiment> create(const int &modus);
    virtual void config(char *path) = 0;
};

template <typename Modus> class ExperimentImplementation : public Experiment
{
public:
    virtual void config(char *path);
//private:
//    BoundaryConditions bc(Parameters);
//    Particles,...

};

#endif  // SRC_INCLUDE_Experiment_H_

