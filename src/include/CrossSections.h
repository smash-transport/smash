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

class CrossSections {
 public:
  /* Add new parametrizations */
  inline void add_pp_elastic(std::vector<float> parametrization);
  inline void add_pp_total(std::vector<float> parametrization);
  inline void add_pn_elastic(std::vector<float> parametrization);
  inline void add_pn_total(std::vector<float> parametrization);
  inline void add_ppbar_elastic(std::vector<float> parametrization);
  inline void add_ppbar_annihilation(std::vector<float> parametrization);
  inline void add_ppbar_total(std::vector<float> parametrization);
 private:
  std::vector< std::vector<float> > pp_elastic_;
  std::vector< std::vector<float> > pp_total_;
  std::vector< std::vector<float> > pn_elastic_;
  std::vector< std::vector<float> > pn_total_;
  std::vector< std::vector<float> > ppbar_elastic_;
  std::vector< std::vector<float> > ppbar_annihilation_;
  std::vector< std::vector<float> > ppbar_total_;
};

/* Add parametrizations */
inline void CrossSections::add_pp_elastic(std::vector<float> parametrization) {
  pp_elastic_.push_back(parametrization);
}
inline void CrossSections::add_pp_total(std::vector<float> parametrization) {
  pp_total_.push_back(parametrization);
}
inline void CrossSections::add_pn_elastic(std::vector<float> parametrization) {
  pn_elastic_.push_back(parametrization);
}
inline void CrossSections::add_pn_total(std::vector<float> parametrization) {
  pn_total_.push_back(parametrization);
}
inline void CrossSections::add_ppbar_elastic(
            std::vector<float> parametrization) {
  ppbar_elastic_.push_back(parametrization);
}

inline void CrossSections::add_ppbar_annihilation(
            std::vector<float> parametrization) {
  ppbar_annihilation_.push_back(parametrization);
}

inline void CrossSections::add_ppbar_total(std::vector<float> parametrization) {
  ppbar_total_.push_back(parametrization);
}


#endif  // SRC_INCLUDE_CROSSSECTIONS_H_
