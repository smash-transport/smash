/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_EXPERIMENT_H_
#define SRC_INCLUDE_EXPERIMENT_H_



class Experiment
{
public:
    virtual void run();
};

template <typename Modus> class ExperimentImplementation : public Experiment
{
public:
  virtual void run();
private:
    BoundaryConditions bc;
    Particles,...

};

#endif  // SRC_INCLUDE_Experiment_H_

