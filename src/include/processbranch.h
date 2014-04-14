/*
 *    Copyright (c) 2013
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_PROCESSBRANCH_H_
#define SRC_INCLUDE_PROCESSBRANCH_H_

#include <vector>

namespace Smash {

/**
 * ProcessBranch represents one possible final state
 * of an interaction process.
 *
 * Each final state has two components;
 * 1. The list of particle types present in this state.
 * 2. The weight of this state, i.e. how probable this outcome is
 * compared to other possible outcomes. Depending on context,
 * this can be either a cross section or a branching ratio.
 */
class ProcessBranch {
 public:
  /// Default constructor
  ProcessBranch() : branch_weight_(-1.0) {}
  /// Add one particle to the list
  inline void add_particle(PdgCode particle_pdg);
  /// Add a complete list of particles
  inline void set_particles(std::vector<PdgCode> particle_pdgs);
  /// Set the weight of the branch,
  /// i.e. how probable it is compared to other branches
  inline void set_weight(float process_weight);
  /// Modify the weight of the branch
  inline void change_weight(float additional_weight);
  /// Clear all information from the branch
  inline void clear(void);
  /// Return the particle list
  inline std::vector<PdgCode> particle_list(void) const;
  /// Return the branch weight
  inline float weight(void) const;

 private:
  /**
   * List of particles appearing in this process outcome.
   *
   * \todo Is there a maximum to the number of particles? If not, a std::vector
   * is fine, otherwise (I assume 4 may be a useful limit?) switch to
   * std::array<int, 4>. std::vector<int> stores one size_t and one pointer,
   * which is as big as 4 ints. And then there's still the data of the vector
   * which is somewhere on the heap. Also the alignment of ints is only half
   * that of size_t/void*. (I was obviously talking about 64bit here...)
   */
  std::vector<PdgCode> particle_list_;
  /// Weight of the branch, typically a cross section or a branching ratio
  float branch_weight_;
};

/// Add one particle to the list
inline void ProcessBranch::add_particle(PdgCode particle_pdg) {
  particle_list_.push_back(particle_pdg);
}

/// Add a complete list of particles
inline void ProcessBranch::set_particles(std::vector<PdgCode> particle_pdgs) {
  particle_list_ = std::move(particle_pdgs);
}

/// Set the weight of the branch,
/// i.e. how probable it is compared to other branches
inline void ProcessBranch::set_weight(float process_weight) {
  branch_weight_ = process_weight;
}

/// Modify the weight of the branch
inline void ProcessBranch::change_weight(float additional_weight) {
  branch_weight_ += additional_weight;
}

/// Clear all information from the branch
inline void ProcessBranch::clear(void) {
  particle_list_.clear();
  branch_weight_ = -1.0;
}

/// Return the particle list
inline std::vector<PdgCode> ProcessBranch::particle_list(void) const {
  return particle_list_;
}

/// Return the branch weight
inline float ProcessBranch::weight(void) const {
  return branch_weight_;
}

}  // namespace Smash

#endif  // SRC_INCLUDE_PROCESSBRANCH_H_
