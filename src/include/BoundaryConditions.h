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

/* forward declartions */
class Particles;
class CrossSections;

class BoundaryConditions
{
public:
 /* default constructor with probable values */
   BoundaryConditions(): steps_(10000), output_interval_(100), testparticles_(1), 
   eps_(0.001f), cross_section_(10.0f), seed_(1) {}
    /* member access funtions */
    inline int steps() const;
    inline void set_steps(int STEPS);
    inline int output_interval() const;
    inline void set_output_interval(int UPDATE);
    inline int testparticles() const;
    inline void set_testparticles(int TESTPARTICLES);
    inline float eps() const;
    inline void set_eps(float EPS);
    inline float cross_section() const;
    inline void set_cross_section(float sigma);
    inline int64_t seed() const;
    inline void set_seed(const int64_t &RANDOMSEED);
    /* special funtion should be called by specific subclass */
    int evolve(Particles *p __attribute__((unused)),
      CrossSections *c __attribute__((unused))) { return -1; }
    void process_config_general(char *path __attribute__((unused))) { return ; }
    void initial_conditions(Particles *p __attribute__((unused))) { return; }
    
  private:
    /* number of steps */
    int steps_;
    /* number of steps before giving measurables */
    int output_interval_;
    /* number of test particle */
    int testparticles_;
    /* temporal time step */
    float eps_;
    /* cross section of the elastic scattering */
    float cross_section_;
    /* initial seed for random generator */
    int64_t seed_;

 };

/* return the number of steps */
inline int BoundaryConditions::steps(void) const {
  return steps_;
}

/* set the number of steps */
inline void BoundaryConditions::set_steps(int STEPS) {
  steps_ = STEPS;
}

/* return the number on which interval output will be shown */
inline int BoundaryConditions::output_interval(void) const {
  return output_interval_;
}

/* set when to output physics */
inline void BoundaryConditions::set_output_interval(int update) {
  output_interval_ = update;
}

/* return the number of testparticles:
 * if equal to one a "testparticle" corresponds to a real particle
 */
inline int BoundaryConditions::testparticles(void) const {
  return testparticles_;
}

/* set the number of test particles */
inline void BoundaryConditions::set_testparticles(int TESTPARTICLES) {
  testparticles_ = TESTPARTICLES;
}

/* return the time step in use */
inline float BoundaryConditions::eps(void) const {
  return eps_;
}

/* set the time step */
inline void BoundaryConditions::set_eps(float EPS) {
  eps_ = EPS;
}

/* return the cross section */
inline float BoundaryConditions::cross_section(void) const {
  return cross_section_;
}

/* set the cross section */
inline void BoundaryConditions::set_cross_section(float sigma) {
  cross_section_ = sigma;
}

/* return the random seed */
inline int64_t BoundaryConditions::seed(void) const {
  return seed_;
}

/* set the seed */
inline void BoundaryConditions::set_seed(const int64_t &randomseed) {
  seed_ = randomseed;
}

 

#endif // BOUNDARYCONDITIONS_H
