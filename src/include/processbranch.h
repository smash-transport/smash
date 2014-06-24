/*
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_PROCESSBRANCH_H_
#define SRC_INCLUDE_PROCESSBRANCH_H_

#include "forwarddeclarations.h"
#include "pdgcode.h"

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
 * std::vector<PdgCode> pdg_list(0x2112, 0x211);
 * // set_particles erases the previous list
 * branch.set_particles(pdg_list);
 * // set_weight erases the previous weight
 * branch.set_weight(1);
 * deltaplus_decay_modes.push_back(branch);
 * \endcode
 */
class ProcessBranch {
 public:
  /// Default constructor
  ProcessBranch() : branch_weight_(-1.0) {}
  /// Constructor with 1 particle
  inline ProcessBranch(PdgCode p, float w);
  /// Constructor with 2 particles
  inline ProcessBranch(PdgCode p1, PdgCode p2, float w);
  /// Constructor with particle vector
  inline ProcessBranch(std::vector<PdgCode> pdg_list, float w);
  /// Add one particle to the list
  inline void add_particle(PdgCode particle_pdg);
  /**
   * Create a particle list.
   * This will remove any previously added particles,
   * but more members can be added to list afterwards
   * with add_particle(int)
   */
  inline void set_particles(std::vector<PdgCode> pdg_list);
  /**
   * Set the weight of the branch.
   * In other words, how probable this branch is
   * compared to other branches
   */
  inline void set_weight(float process_weight);
  /// Clear all information from the branch
  inline void clear(void);
  /// Return the particle list
  inline std::vector<PdgCode> pdg_list(void) const;

  /**
   * Return a list of ParticleData initialized with the ParticleType for the PDG
   * codes from pdg_list.
   */
  ParticleList particle_list() const;

  /// Return the branch weight
  inline float weight(void) const;

  /**
   * Determine the threshold for this branch, i.e. the minimum energy that is
   * required to produce all final-state particles.
   */
  float threshold() const;

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
  std::vector<PdgCode> pdg_list_;
  /// Weight of the branch, typically a cross section or a branching ratio
  float branch_weight_;
};

// Constructor with 1 particle
ProcessBranch::ProcessBranch (PdgCode p, float w) : branch_weight_(w) {
  add_particle(p);
}

// Constructor with 2 particles
ProcessBranch::ProcessBranch (PdgCode p1, PdgCode p2, float w)
                          : branch_weight_(w) {
  add_particle(p1);
  add_particle(p2);
}

// Constructor with particle vector
ProcessBranch::ProcessBranch (std::vector<PdgCode> new_pdg_list, float w)
                          : branch_weight_(w) {
  set_particles(std::move(new_pdg_list));
}

/// Add one particle to the list
inline void ProcessBranch::add_particle(PdgCode particle_pdg) {
  pdg_list_.push_back(particle_pdg);
}

/**
 * Create a particle list.
 * This will remove any previously added particles,
 * but more members can be added to list afterwards
 * with add_particle(int)
 */
inline void ProcessBranch::set_particles(std::vector<PdgCode> new_pdg_list) {
  pdg_list_ = std::move(new_pdg_list);
}

/**
 * Set the weight of the branch.
 * In other words, how probable this branch is
 * compared to other branches
 */
inline void ProcessBranch::set_weight(float process_weight) {
  branch_weight_ = process_weight;
}

/// Clear all information from the branch
inline void ProcessBranch::clear(void) {
  pdg_list_.clear();
  branch_weight_ = -1.0;
}

/// Return the particle list
inline std::vector<PdgCode> ProcessBranch::pdg_list(void) const {
  return pdg_list_;
}

/// Return the branch weight
inline float ProcessBranch::weight(void) const {
  return branch_weight_;
}

/**
 * DecayBranch is a derivative of ProcessBranch,
 * which is used to represent decay channels.
 * It contains additional information like the angular momentum.
 */
class DecayBranch : public ProcessBranch {
 public:
  inline int angular_momentum() const;
  inline void set_angular_momentum(const int L);
 private:
  int angular_momentum_;
};

inline int DecayBranch::angular_momentum() const {
  return angular_momentum_;
};

inline void DecayBranch::set_angular_momentum(const int L)
{
  angular_momentum_ = L;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_PROCESSBRANCH_H_
