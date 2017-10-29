/*
 *    Copyright (c) 2013-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_PROCESSBRANCH_H_
#define SRC_INCLUDE_PROCESSBRANCH_H_

#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include "decaytype.h"
#include "forwarddeclarations.h"
#include "particletype.h"

namespace Smash {

/** Process Types are used to identify the type of the process,
 * currently we have 12 of these:
 * (0) nothing (None)
 * (1) elastic (Elastic)
 * (2) resonance formation (2->1) (TwoToOne)
 * (3) 2->2 (inelastic) (TwoToTwo)
 * (4) string excitation by PYTHIA (String)
 * (5) resonance decays (Decay)
 * (6) Wall transition (Wall)
 * (7) Forces thermalization
 * (41) Soft string excitation
 * (42) Hard string process involving 2->2 QCD process
 */
enum class ProcessType {
  None = 0,
  Elastic = 1,
  TwoToOne = 2,
  TwoToTwo = 3,
  String = 4,
  Decay = 5,
  Wall = 6,
  Thermalization = 7,
  StringSoft = 41,
  StringHard = 42
};

std::ostream &operator<<(std::ostream &os, ProcessType process_type);

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
 * 3. The process id that identifies a certain class of processes.
 * If the outgoing particles are not known yet, e.g. for strings
 * there will be only the weight and the process id.
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
  ProcessBranch() : branch_weight_(0.) {}
  explicit ProcessBranch(double w) : branch_weight_(w) {}

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
  inline void set_weight(double process_weight);
  /// Return the process type
  virtual ProcessType get_type() const = 0;

  /// Return the particle types associated with this branch.
  virtual const ParticleTypePtrList &particle_types() const = 0;

  /**
   * Return a list of ParticleData initialized with the stored ParticleType
   * objects.
   */
  ParticleList particle_list() const;

  /// Return the branch weight
  inline double weight() const;

  /**
   * Determine the threshold for this branch, i.e. the minimum energy that is
   * required to produce all final-state particles.
   */
  double threshold() const;

  virtual unsigned int particle_number() const = 0;

 protected:
  /// Weight of the branch, typically a cross section or a branching ratio
  double branch_weight_;
  /// Threshold of the branch
  mutable double threshold_ = -1.;
};

/**
 * Set the weight of the branch.
 * In other words, how probable this branch is
 * compared to other branches
 */
inline void ProcessBranch::set_weight(double process_weight) {
  branch_weight_ = process_weight;
}

/// Return the branch weight
inline double ProcessBranch::weight() const { return branch_weight_; }

/** \relates ProcessBranch
 * Calculates the total weight by summing all weights of the ProcessBranch
 * objects in the list \p l.
 */
template <typename Branch>
inline double total_weight(const ProcessBranchList<Branch> &l) {
  double sum = 0.;
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
  CollisionBranch(double w, ProcessType p_type)
      : ProcessBranch(w), process_type_(p_type) {}
  /// Constructor with 1 particle
  CollisionBranch(const ParticleType &type, double w, ProcessType p_type)
      : ProcessBranch(w), process_type_(p_type) {
    particle_types_.reserve(1);
    particle_types_.push_back(&type);
  }
  /// Constructor with 2 particles
  CollisionBranch(const ParticleType &type_a, const ParticleType &type_b,
                  double w, ProcessType p_type)
      : ProcessBranch(w), process_type_(p_type) {
    particle_types_.reserve(2);
    particle_types_.push_back(&type_a);
    particle_types_.push_back(&type_b);
  }
  /// Constructor with a list of particles
  CollisionBranch(ParticleTypePtrList new_types, double w, ProcessType p_type)
      : ProcessBranch(w),
        particle_types_(std::move(new_types)),
        process_type_(p_type) {}
  /// The move constructor efficiently moves the particle-type list member.
  CollisionBranch(CollisionBranch &&rhs)
      : ProcessBranch(rhs.branch_weight_),
        particle_types_(std::move(rhs.particle_types_)),
        process_type_(rhs.process_type_) {}
  /// Return the particle types associated with this branch.
  const ParticleTypePtrList &particle_types() const override {
    return particle_types_;
  }
  /// Set the process type
  inline void set_type(ProcessType p_type) { process_type_ = p_type; }
  /// Return the process type
  inline ProcessType get_type() const override { return process_type_; }
  unsigned int particle_number() const override {
    return particle_types_.size();
  }

 private:
  /**
   * List of particles appearing in this process outcome.
   *
   * \note This currently uses a std::vector and thus works for any number of
   * particles. But this number is bounded (4?) and a std::array may therefore
   * be more efficient.
   */
  ParticleTypePtrList particle_types_;
  /**
   * Process type are used to distinguish different types of processes,
   * e.g. string formation, resonance formation, elastic scattering and so on.
   *
   */
  ProcessType process_type_;
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
  DecayBranch(const DecayType &t, double w) : ProcessBranch(w), type_(t) {}
  /// The move constructor efficiently moves the particle-type list member.
  DecayBranch(DecayBranch &&rhs)
      : ProcessBranch(rhs.branch_weight_), type_(rhs.type_) {}
  /// Get the angular momentum of this branch.
  inline int angular_momentum() const { return type_.angular_momentum(); }
  /// Return the particle types associated with this branch.
  const ParticleTypePtrList &particle_types() const override {
    return type_.particle_types();
  }
  unsigned int particle_number() const override {
    return type_.particle_number();
  }
  inline const DecayType &type() const { return type_; }
  /// Return the process type
  inline ProcessType get_type() const override { return ProcessType::Decay; }

 private:
  // decay type (including final-state particles and angular momentum
  const DecayType &type_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_PROCESSBRANCH_H_
