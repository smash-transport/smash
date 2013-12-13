/*
 *
 *    Copyright (c) 2012-2013
 *      Hannah Petersen <petersen@fias.uni-frankfurt.de>
 *      
 *    GNU General Public License (GPLv3)
 *
 */
 
 
#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include <iostream>
#include <stdint.h>
#include <cmath>
#include <list>
#include "../include/Parameters.h"

/* forward declartions */
class Particles;
class CrossSections;

class BoundaryConditions
{
public:
 /* default constructor with probable values */
   BoundaryConditions(): steps(10000), output_interval(100), testparticles(1),
   eps(0.001f), cross_section(10.0f), seed(1) {}
    /* special funtion should be called by specific subclass */
    virtual int evolve(Particles *p __attribute__((unused)),
      CrossSections *c __attribute__((unused))) { return -1; }
    virtual void initial_conditions(Particles *p __attribute__((unused))) { return; }
    virtual void assign_params(std::list<Parameters> *configuration);
    virtual void print_startup();
  public:
    /* number of steps */
    int steps;
    /* number of steps before giving measurables */
    int output_interval;
    /* number of test particle */
    int testparticles;
    /* temporal time step */
    float eps;
    /* cross section of the elastic scattering */
    float cross_section;
    /* initial seed for random generator */
    int64_t seed;

 };



#endif // BOUNDARYCONDITIONS_H
