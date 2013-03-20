/*
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_BOX_H_
#define SRC_INCLUDE_BOX_H_

class box {
  public:
    void inline set(const int STEPS, const int UPDATE,
      const float A, const float EPS, const float TEMPERATURE);
    float inline a();
    void inline set_a(const float &A);
    float inline eps();
    void inline set_eps(const float &EPS);
    int inline steps();
    void inline set_steps(const int &STEPS);
    int inline update();
    void inline set_update(const int &UPDATE);
    float inline temperature();
    void inline set_temperature(const float &T);

  private:
    /* number of steps */
    int steps_;
    /* number of steps before giving measurables */
    int update_;
    /* Cube edge length */
    float a_;
    /* temporal time step */
    float eps_;
    /* Temperature of the Boltzmann distribution for thermal initialization */
    float temperature_;
};

void inline box::set(const int STEPS, const int UPDATE, const float A,
      const float EPS, const float TEMPERATURE) {
  steps_ = STEPS;
  update_ = UPDATE;
  a_ = A;
  eps_ = EPS;
  temperature_ = TEMPERATURE;
}

float inline box::a(void) {
  return a_;
}

void inline box::set_a(const float &A) {
  a_ = A;
}

float inline box::eps(void) {
  return eps_;
}

void inline box::set_eps(const float &EPS) {
  eps_ = EPS;
}

int inline box::steps(void) {
  return steps_;
}

void inline box::set_steps(const int &STEPS) {
  steps_ = STEPS;
}

int inline box::update(void) {
  return update_;
}

void inline box::set_update(const int &UPDATE) {
  update_ = UPDATE;
}

float inline box::temperature(void) {
  return temperature_;
}

void inline box::set_temperature(const float &T) {
  temperature_ = T;
}

/* support for gcc branch prediction */
#ifdef __GNUC__
#define likely(x)       __builtin_expect((x), 1)
#define unlikely(x)     __builtin_expect((x), 0)
#else
#define likely(x)       (x)
#define unlikely(x)     (x)
#endif

/*XXX: get rid of those globals */
/* Default random seed */
extern unsigned int seedp;

#endif  // SRC_INCLUDE_BOX_H_
