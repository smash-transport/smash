/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_LABORATORY_H_
#define SRC_INCLUDE_LABORATORY_H_

class Laboratory {
  public:
    /* default constructor with probable values */
    Laboratory(): steps_(10000), output_interval_(100), initial_condition_(1),
      testparticles_(1), eps_(0.001f), cross_section_(10.0f), seed_(1) {}
    /* member funtions */
    inline int steps() const;
    inline void set_steps(const int &STEPS);
    inline int output_interval() const;
    inline void set_output_interval(const int &UPDATE);
    inline int modus() const;
    inline void set_modus(const int &MODUS);
    inline int initial_condition() const;
    inline void set_initial_condition(const int &INITIAL_CONDITION);
    inline int testparticles() const;
    inline void set_testparticles(const int &TESTPARTICLES);
    inline float eps() const;
    inline void set_eps(const float &EPS);
    inline float cross_section() const;
    inline void set_cross_section(const float &sigma);
    inline int64_t seed() const;
    inline void set_seed(const int64_t &RANDOMSEED);

  private:
    /* number of steps */
    int steps_;
    /* number of steps before giving measurables */
    int output_interval_;
    /* box or sphere or .. */
    int modus_;
    /* initial condition */
    int initial_condition_;
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
inline int Laboratory::steps(void) const {
  return steps_;
}

/* set the number of steps */
inline void Laboratory::set_steps(const int &STEPS) {
  steps_ = STEPS;
}

/* return the number on which interval output will be shown */
inline int Laboratory::output_interval(void) const {
  return output_interval_;
}

/* set when to output physics */
inline void Laboratory::set_output_interval(const int &update) {
  output_interval_ = update;
}

/* return the modus in use */
inline int Laboratory::modus(void) const {
  return modus_;
}

/* set the modus */
inline void Laboratory::set_modus(const int &MODUS) {
  modus_ = MODUS;
}

/* return the used initial condition */
inline int Laboratory::initial_condition(void) const {
  return initial_condition_;
}

/* set the initial condition */
inline void Laboratory::set_initial_condition(const int &INITIAL_CONDITION) {
  initial_condition_ = INITIAL_CONDITION;
}


/* return the number of testparticles:
 * if equal to one a "testparticle" corresponds to a real particle
 */
inline int Laboratory::testparticles(void) const {
  return testparticles_;
}

/* set the number of test particles */
inline void Laboratory::set_testparticles(const int &TESTPARTICLES) {
  testparticles_ = TESTPARTICLES;
}

/* return the time step in use */
inline float Laboratory::eps(void) const {
  return eps_;
}

/* set the time step */
inline void Laboratory::set_eps(const float &EPS) {
  eps_ = EPS;
}

/* return the cross section */
inline float Laboratory::cross_section(void) const {
  return cross_section_;
}

/* set the cross section */
inline void Laboratory::set_cross_section(const float &sigma) {
  cross_section_ = sigma;
}

/* return the random seed */
inline int64_t Laboratory::seed(void) const {
  return seed_;
}

/* set the seed */
inline void Laboratory::set_seed(const int64_t &randomseed) {
  seed_ = randomseed;
}

#endif  // SRC_INCLUDE_LABORATORY_H_
