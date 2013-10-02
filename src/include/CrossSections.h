/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_CROSSSECTIONS_H_
#define SRC_INCLUDE_CROSSSECTIONS_H_

#include <vector>

#include "../include/Particles.h"

class CrossSections {
 public:
  /* Add the values to be used in parametrizations */
  inline void add_elastic_parameter(const float parameter);
  inline void add_pp_elastic(std::vector<float> parameters);
  inline void add_pp_total(std::vector<float> parameters);
  inline void add_pn_elastic(std::vector<float> parameters);
  inline void add_pn_total(std::vector<float> parameters);
  inline void add_ppbar_elastic(std::vector<float> parameters);
  inline void add_ppbar_annihilation(std::vector<float> parameters);
  inline void add_ppbar_total(std::vector<float> parameters);
  /* Compute kinematics */
  void compute_kinematics(Particles *particles, int id_a, int id_b);
  /* Return cross sections */
  float elastic(Particles *particles, int id_a, int id_b) const;
  float annihilation(Particles *particles, int id_a, int id_b) const;
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
  /* Member containing the parametrizations */
  float parametrization_(std::vector<float> parameters) const;
  /* Values to be used in parametrizations */
  std::vector< std::vector<float> > pp_elastic_;
  std::vector< std::vector<float> > pp_total_;
  std::vector< std::vector<float> > pn_elastic_;
  std::vector< std::vector<float> > pn_total_;
  std::vector< std::vector<float> > ppbar_elastic_;
  std::vector< std::vector<float> > ppbar_annihilation_;
  std::vector< std::vector<float> > ppbar_total_;
};

/* Add the values to be used in parametrizations */
inline void CrossSections::add_elastic_parameter(const float parameter) {
  elastic_parameter_ = parameter;
}

inline void CrossSections::add_pp_elastic(std::vector<float> parameters) {
  pp_elastic_.push_back(parameters);
}

inline void CrossSections::add_pp_total(std::vector<float> parameters) {
  pp_total_.push_back(parameters);
}

inline void CrossSections::add_pn_elastic(std::vector<float> parameters) {
  pn_elastic_.push_back(parameters);
}

inline void CrossSections::add_pn_total(std::vector<float> parameters) {
  pn_total_.push_back(parameters);
}
inline void CrossSections::add_ppbar_elastic(
            std::vector<float> parameters) {
  ppbar_elastic_.push_back(parameters);
}

inline void CrossSections::add_ppbar_annihilation(
            std::vector<float> parameters) {
  ppbar_annihilation_.push_back(parameters);
}

inline void CrossSections::add_ppbar_total(std::vector<float> parameters) {
  ppbar_total_.push_back(parameters);
}

#endif  // SRC_INCLUDE_CROSSSECTIONS_H_
