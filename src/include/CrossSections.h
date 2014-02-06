/*
 *    Copyright (c) 2013
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_CROSSSECTIONS_H_
#define SRC_INCLUDE_CROSSSECTIONS_H_

#include "include/Particles.h"

class CrossSections {
 public:
  /* Default constructor */
  CrossSections() : elastic_parameter_(0.0), squared_mass_a_(-1.0),
                    squared_mass_b_(-1.0), mandelstam_s_(-1.0), p_lab_(-1.0) {}
  /* Fixed elastic cross section value */
  inline void add_elastic_parameter(float parameter);
  /* Compute kinematics */
  void compute_kinematics(Particles *particles, int id_a, int id_b);
  /* Return cross sections */
  float elastic(Particles *particles, int id_a, int id_b) const;
  float total(Particles *particles, int id_a, int id_b) const;

 private:
  /* Elastic cross section parameter */
  float elastic_parameter_;
  /* Mass of the first particle */
  float squared_mass_a_;
  /* Mass of the second particle */
  float squared_mass_b_;
  /* Mandelstam s of the collision (= CMS energy squared) */
  double mandelstam_s_;
  /* "Beam" momentum */
  double p_lab_;
};

/* Fixed elastic cross section value */
inline void CrossSections::add_elastic_parameter(float parameter) {
  elastic_parameter_ = parameter;
}


#endif  // SRC_INCLUDE_CROSSSECTIONS_H_
