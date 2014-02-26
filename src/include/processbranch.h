/*
 *    Copyright (c) 2013
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_PROCESSBRANCH_H_
#define SRC_INCLUDE_PROCESSBRANCH_H_

#include <vector>

class ProcessBranch {
 public:
  /* Default constructor */
  ProcessBranch() : branch_weight_(-1.0) {}
  /* Add a particle to the list */
  inline void add_particle(int particle_pdg);
  /* Add the complete particle list */
  inline void set_particles(std::vector<int> particle_pdgs);
  /* Add branch ratio */
  inline void set_weight(double process_weight);
  /* Add to the ratio of this branch */
  inline void change_weight(double additional_weight);
  /* Remove all modes */
  inline void clear(void);
  /* Pass the particle list */
  inline std::vector<int> particle_list(void) const;
  /* Pass the branch ratio */
  inline double weight(void) const;

 private:
  /**
   * \todo Is there a maximum to the number of branches? If not, a std::vector
   * is fine, otherwise (I assume 4 may be a useful limit?) switch to
   * std::array<int, 4>. std::vector<int> stores one size_t and one pointer,
   * which is as big as 4 ints. And then there's still the data of the vector
   * which is somewhere on the heap. Also the alignment of ints is only half
   * that of size_t/void*. (I was obviously talking about 64bit here...)
   */
  std::vector<int> particle_list_;
  /**
   * \todo This value is initialized from a single-precision float. Why is it a
   * double then?
   */
  double branch_weight_;
};

/* Add a particle to the list */
inline void ProcessBranch::add_particle(int particle_pdg) {
  particle_list_.push_back(particle_pdg);
}

/* Add the complete particle list */
inline void ProcessBranch::set_particles(std::vector<int> particle_pdgs) {
  particle_list_ = std::move(particle_pdgs);
}

/* Add the ratio of this branch */
inline void ProcessBranch::set_weight(double process_weight) {
  branch_weight_ = process_weight;
}

/* Add to the ratio of this branch */
inline void ProcessBranch::change_weight(double additional_weight) {
  branch_weight_ += additional_weight;
}

/* Remove all modes */
inline void ProcessBranch::clear(void) {
  particle_list_.clear();
  branch_weight_ = -1.0;
}

/* Pass the particle list */
inline std::vector<int> ProcessBranch::particle_list(void) const {
  return particle_list_;
}

/* Pass the branch ratio */
inline double ProcessBranch::weight(void) const {
  return branch_weight_;
}

#endif  // SRC_INCLUDE_PROCESSBRANCH_H_
