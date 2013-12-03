/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_PROCESSBRANCH_H_
#define SRC_INCLUDE_PROCESSBRANCH_H_

#include <vector>

class ProcessBranch {
 public:
  /* Add a particle to the list */
  inline void add_particle(int particle_pdg);
  /* Add the complete particle list */
  inline void add_particles(std::vector<int> particle_pdgs);
  /* Add branch ratio */
  inline void add_ratio(double ratio);
  /* Remove all modes */
  inline void clear(void);
  /* Pass the particle list */
  inline std::vector<int> particle_list(void);
  /* Pass the branch ratio */
  inline double ratio(void);
 private:
  std::vector<int> particle_list_;
  double branch_ratio_;
};

/* Add a particle to the list */
inline void ProcessBranch::add_particle(int particle_pdg) {
  particle_list_.push_back(particle_pdg);
}

/* Add the complete particle list */
inline void ProcessBranch::add_particles(std::vector<int> particle_pdgs) {
  particle_list_ = particle_pdgs;
}

/* Add the ratio of this branch */
inline void ProcessBranch::add_ratio(double ratio) {
  branch_ratio_ = ratio;
}

/* Remove all modes */
inline void ProcessBranch::clear(void) {
  particle_list_.clear();
  branch_ratio_ = -1.0;
}

/* Pass the particle list */
inline std::vector<int> ProcessBranch::particle_list(void) {
  return particle_list_;
}

/* Pass the branch ratio */
inline double ProcessBranch::ratio(void) {
  return branch_ratio_;
}

#endif  // SRC_INCLUDE_PROCESSBRANCH_H_
