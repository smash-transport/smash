/*
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_PROCESSBRANCH_H_
#define SRC_INCLUDE_PROCESSBRANCH_H_

#include "forwarddeclarations.h"
#include "particletype.h"
#include "decaytype.h"

#include <vector>
#include <memory>

namespace Smash {

/**
 * \ingroup data
 *
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
 * branch.set_weight(2);
 * deltaplus_decay_modes.push_back(branch);
 * // set_weight erases the previous weight
 * branch.set_weight(1);
 * deltaplus_decay_modes.push_back(branch);
 * \endcode
 */
class ProcessBranch {
 public:
  /// Create a ProcessBranch without final states
  ProcessBranch() : branch_weight_(-1.0) {}
  ProcessBranch(float w) : branch_weight_(w) {}

  /// Copying is disabled. Use std::move or create a new object.
  ProcessBranch(const ProcessBranch &) = delete;

  /** Virtual Destructor.
   * The declaration of the destructor is necessary to make it virtual.
   */
  virtual ~ProcessBranch() = default;

  /**
   * Set the weight of the branch.
   * In other words, how probable this branch is
   * compared to other branches
   */
  inline void set_weight(float process_weight);
  /// Clear all information from the branch
  inline void clear(void);

  /// Return the particle types associated with this branch.
  virtual const ParticleTypePtrList &particle_types() const = 0;

  /**
   * Return a list of ParticleData initialized with the stored ParticleType
   * objects.
   */
  ParticleList particle_list() const;

  /// Return the branch weight
  inline float weight(void) const;

  /**
   * Determine the threshold for this branch, i.e. the minimum energy that is
   * required to produce all final-state particles.
   */
  float threshold() const;

  virtual unsigned int particle_number() const = 0;
 protected:
  /// Weight of the branch, typically a cross section or a branching ratio
  float branch_weight_;
};

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
  branch_weight_ = -1.0;
}

/// Return the branch weight
inline float ProcessBranch::weight(void) const {
  return branch_weight_;
}

/** \relates ProcessBranch
 * Calculates the total weight by summing all weights of the ProcessBranch
 * objects in the list \p l.
 */
inline float total_weight(const ProcessBranchList &l) {
  float sum = 0.f;
  for (const auto &p : l) {
    sum += p->weight();
  }
  return sum;
}


/**
 * \ingroup data
 *
 * CollisionBranch is a derivative of ProcessBranch,
 * which is used to represent particular final-state channels in a collision.
 */
class CollisionBranch : public ProcessBranch {
 public:
   CollisionBranch() : ProcessBranch() {}
  /// Constructor with 1 particle
  CollisionBranch(const ParticleType &type, float w) : ProcessBranch(w) {
    particle_types_.reserve(1);
    particle_types_.push_back(&type);
  }
  /// Constructor with 2 particles
  CollisionBranch(const ParticleType &type_a, const ParticleType &type_b, float w)
      : ProcessBranch(w) {
    particle_types_.reserve(2);
    particle_types_.push_back(&type_a);
    particle_types_.push_back(&type_b);
  }
  /// Constructor with a list of particles
  CollisionBranch(ParticleTypePtrList new_types, float w)
      : ProcessBranch(w), particle_types_(std::move(new_types)) {}
  /// The move constructor efficiently moves the particle-type list member.
  CollisionBranch(CollisionBranch &&rhs)
      : ProcessBranch(rhs.branch_weight_),
        particle_types_(std::move(rhs.particle_types_)) {}
  /// Return the particle types associated with this branch.
  const ParticleTypePtrList &particle_types() const override {
    return particle_types_;
  }
  /// Clear all information from the branch
  inline void clear(void) {
    particle_types_.clear();
    ProcessBranch::clear();
  }
  unsigned int particle_number() const override {
    return particle_types_.size();
  }
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
  ParticleTypePtrList particle_types_;
};


/**
 * \ingroup data
 *
 * DecayBranch is a derivative of ProcessBranch,
 * which is used to represent decay channels.
 * It contains additional information like the angular momentum.
 */
class DecayBranch : public ProcessBranch {
 public:
   DecayBranch(const DecayType *t, float w) : ProcessBranch(w), type_(t) {}
  /// The move constructor efficiently moves the particle-type list member.
  DecayBranch(DecayBranch &&rhs) : ProcessBranch(rhs.branch_weight_),
                                   type_(rhs.type_) {}
  /// Get the angular momentum of this branch.
  inline int angular_momentum() const;
  /// Return the particle types associated with this branch.
  const ParticleTypePtrList &particle_types() const override {
    return type_->particle_types();
  }
  unsigned int particle_number() const override {
    return type_->particle_number();
  }
  inline const DecayType& type() const {
    return *type_;
  }
  /// Clear all information from the branch
  inline void clear(void) {
    delete type_;
    ProcessBranch::clear();
  }
 private:
  // decay type (including final-state particles and angular momentum
  const DecayType *type_;
};

inline int DecayBranch::angular_momentum() const {
  return type_->angular_momentum();
};


}  // namespace Smash

#endif  // SRC_INCLUDE_PROCESSBRANCH_H_
