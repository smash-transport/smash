/*
 *    Copyright (c) 2013-2014
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
 *
 * For example, create a list of decay modes for \f$\Delta^+\f$ resonance:
 * \code
 * std::vector<ProcessBranch> deltaplus_decay_modes;
 * ProcessBranch branch;
 * // Adding particle codes one by one
 * branch.add_particles(0x2212);
 * branch.add_particles(0x111)
 * branch.set_weight(2);
 * deltaplus_decay_modes.push_back(branch);
 * // Using already existing particle list
 * std::vector<PdgCode> particle_list(0x2112, 0x211);
 * // set_particles erases the previous list
 * branch.set_particles(particle_list);
 * // set_weight erases the previous weight
 * branch.set_weight(1);
 * deltaplus_decay_modes.push_back(branch);
 * \endcode
 */
class ProcessBranch {
 public:
  /// Default constructor
  ProcessBranch() : branch_weight_(-1.0) {}
  // Constructor with 1 particle
  inline ProcessBranch(PdgCode p, float w, int t);
  // Constructor with 2 particles
  inline ProcessBranch(PdgCode p1, PdgCode p2, float w, int t);
  /// Add one particle to the list
  inline void add_particle(PdgCode particle_pdg);
  /**
   * Create a particle list.
   * This will remove any previously added particles,
   * but more members can be added to list afterwards
   * with add_particle(int)
   */
  inline void set_particles(std::vector<PdgCode> particle_pdgs);
  /**
   * Set the weight of the branch.
   * In other words, how probable this branch is
   * compared to other branches
   */
  inline void set_weight(float process_weight);
  /// Set the type of interaction
  inline void set_type(int t);
  /// Clear all information from the branch
  inline void clear(void);
  /// Return the particle list
  inline std::vector<PdgCode> particle_list(void) const;
  /// Return the branch weight
  inline float weight(void) const;
  /// Return the type of interaction
  inline int type(void) const;

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
  /// Type of interaction
  int interaction_type_;
};

// Constructor with 1 particle
ProcessBranch::ProcessBranch (PdgCode p, float w, int t)
                        : branch_weight_(w), interaction_type_(t) {
  this->add_particle(p);
}

// Constructor with 2 particles
ProcessBranch::ProcessBranch (PdgCode p1, PdgCode p2, float w, int t)
                          : branch_weight_(w), interaction_type_(t) {
  this->add_particle(p1);
  this->add_particle(p2);
}

/// Add one particle to the list
inline void ProcessBranch::add_particle(PdgCode particle_pdg) {
  particle_list_.push_back(particle_pdg);
}

/**
 * Create a particle list.
 * This will remove any previously added particles,
 * but more members can be added to list afterwards
 * with add_particle(int)
 */
inline void ProcessBranch::set_particles(std::vector<PdgCode> particle_pdgs) {
  particle_list_ = std::move(particle_pdgs);
}

/**
 * Set the weight of the branch.
 * In other words, how probable this branch is
 * compared to other branches
 */
inline void ProcessBranch::set_weight(float process_weight) {
  branch_weight_ = process_weight;
}

/// Set the type of interaction.
inline void ProcessBranch::set_type (int t) {
  interaction_type_ = t;
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

/// Return the type of interaction
inline int ProcessBranch::type (void) const {
  return interaction_type_;
}

}  // namespace Smash

#endif  // SRC_INCLUDE_PROCESSBRANCH_H_
