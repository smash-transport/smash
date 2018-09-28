/*
 *    Copyright (c) 2013-2018
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

namespace smash {

/**
 * Process Types are used to identify the type of the process. Corresponding
 * integer numbers are given explicitly, because they appear in the output.
 */
enum class ProcessType {
  /// nothing
  None = 0,
  /// elastic scattering: particles remain the same, only momenta change
  Elastic = 1,
  /// resonance formation (2->1)
  TwoToOne = 2,
  /// 2->2 inelastic scattering
  TwoToTwo = 3,
  /// resonance decays
  Decay = 5,
  /// box wall crossing
  Wall = 6,
  /**
   *  forced thermalization, many particles are replaced
   *  by a thermalized ensemble
   */
  Thermalization = 7,
  /**
   * (41-45) soft string excitations. Here "soft" means that the process does
   *         not involve quark or gluon scattering. A string is formed by quark
   *         and antiquark, or quark and diquark, in its ends. Then this
   *         string decays. Depending on which quark and anti- (or di-)quarks
   *         are selected for string formation, the process has one of the
   *         following types.
   */
  /// single diffractive AB->AX. Both quark and anti-/di-quark taken from B.
  StringSoftSingleDiffractiveAX = 41,
  /**
   *  single diffractive AB->XB. Both quark and anti-/di-quark taken from A.
   *  It makes sense to distinguish it from AB->AX, because A and B can be
   *  particles of different types, for example, a pion and a proton.
   *  It matters then, whether the pion creates a string or the proton.
   */
  StringSoftSingleDiffractiveXB = 42,
  /// double diffractive. Two strings are formed, one from A and one from B.
  StringSoftDoubleDiffractive = 43,
  /**
   * a special case of baryon-antibaryon annihilation. One pair qqbar
   * annihilates immediately and then two strings are formed.
   */
  StringSoftAnnihilation = 44,
  /// non-diffractive. Two strings are formed both have ends in A and B.
  StringSoftNonDiffractive = 45,
  /**
   *  hard string process involving 2->2 QCD process by PYTHIA. Here quarks
   *  do not simply form a string. They actually scatter on parton level first.
   */
  StringHard = 46
};

/**
 * Check if a given process type is a soft string excitation
 * \param[in] p The process type
 */
bool is_string_soft_process(ProcessType p);

/**
 * \ingroup logging
 * Writes the textual representation of the \p process_type
 * to the output stream \p os.
 */
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
  /// Create a ProcessBranch without final states and weight.
  ProcessBranch() : branch_weight_(0.) {}

  /**
   * Create a ProcessBranch with with weight but without final states.
   * \param[in] w Weight of new branch.
   */
  explicit ProcessBranch(double w) : branch_weight_(w) {}

  /// Copying is disabled. Use std::move or create a new object.
  ProcessBranch(const ProcessBranch &) = delete;

  /**
   * Virtual Destructor.
   * The declaration of the destructor is necessary to make it virtual.
   */
  virtual ~ProcessBranch() = default;

  /**
   * Set the weight of the branch.
   * In other words, how probable this branch is
   * compared to other branches.
   * \param[in] process_weight Weight of the process.
   */
  inline void set_weight(double process_weight);

  /// \return the process type
  virtual ProcessType get_type() const = 0;

  /// \return the particle types associated with this branch.
  virtual const ParticleTypePtrList &particle_types() const = 0;

  /**
   * \return a list of ParticleData initialized with the stored ParticleType
   * objects.
   */
  ParticleList particle_list() const;

  /// \return the branch weight.
  inline double weight() const;

  /**
   * \return the threshold for this branch, i.e. the minimum energy that is
   * required to produce all final-state particles.
   */
  double threshold() const;

  /// \return the number of particles in the final state.
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
 * compared to other branches.
 * \param[in] process_weight weight of the branch
 */
inline void ProcessBranch::set_weight(double process_weight) {
  branch_weight_ = process_weight;
}

/// \return the branch weight
inline double ProcessBranch::weight() const { return branch_weight_; }

/**
 * \relates ProcessBranch
 * \param[in] l The list of all the processes that would be summed.
 * \return the total weight calculated by summing all weights of the
 * ProcessBranch objects in the list \p l.
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
  /**
   * Construct collision branch with empty final state.
   * \param[in] w Weight of created branch.
   * \param[in] p_type Process type of created branch.
   */
  CollisionBranch(double w, ProcessType p_type)
      : ProcessBranch(w), process_type_(p_type) {}
  /**
   * Construct collision branch with 1 particle in final state.
   * \param[in] type Particle type of final state particle.
   * \param[in] w Weight of created branch.
   * \param[in] p_type Process type of created branch.
   */
  CollisionBranch(const ParticleType &type, double w, ProcessType p_type)
      : ProcessBranch(w), process_type_(p_type) {
    particle_types_.reserve(1);
    particle_types_.push_back(&type);
  }
  /**
   * Construct collision branch with 2 particles in final state.
   * \param[in] type_a Particle types of one final state particle.
   * \param[in] type_b Particle types of other final state particle.
   * \param[in] w Weight of created branch.
   * \param[in] p_type Process type of created branch.
   */
  CollisionBranch(const ParticleType &type_a, const ParticleType &type_b,
                  double w, ProcessType p_type)
      : ProcessBranch(w), process_type_(p_type) {
    particle_types_.reserve(2);
    particle_types_.push_back(&type_a);
    particle_types_.push_back(&type_b);
  }
  /**
   * Construct collision branch with a list of particles in final state.
   * \param[in] new_types List of particle types of final state particles.
   * \param[in] w Weight of created branch.
   * \param[in] p_type Process type of created branch.
   */
  CollisionBranch(ParticleTypePtrList new_types, double w, ProcessType p_type)
      : ProcessBranch(w),
        particle_types_(std::move(new_types)),
        process_type_(p_type) {}
  /// The move constructor efficiently moves the particle-type list member.
  CollisionBranch(CollisionBranch &&rhs)
      : ProcessBranch(rhs.branch_weight_),
        particle_types_(std::move(rhs.particle_types_)),
        process_type_(rhs.process_type_) {}
  const ParticleTypePtrList &particle_types() const override {
    return particle_types_;
  }
  /**
   * Set the process type
   *
   * \param[in] p_type The new value of the process type
   */
  inline void set_type(ProcessType p_type) { process_type_ = p_type; }
  /// \return type of the process
  inline ProcessType get_type() const override { return process_type_; }
  /// \return number of particles involved in the process
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
   */
  ProcessType process_type_;
};

/**
 * \ingroup logging
 * Writes the textual representation of the Collision Branch \p cbranch
 * to the output stream \p os.
 */
std::ostream &operator<<(std::ostream &os, const CollisionBranch &cbranch);

/**
 * \ingroup data
 *
 * DecayBranch is a derivative of ProcessBranch,
 * which is used to represent decay channels.
 * It contains additional information like the angular momentum.
 */
class DecayBranch : public ProcessBranch {
 public:
  /**
   * Construct decay branch.
   * \param[in] t DecayType of branch.
   * \param[in] w Weight of created branch.
   */
  DecayBranch(const DecayType &t, double w) : ProcessBranch(w), type_(t) {}
  /// The move constructor efficiently moves the particle-type list member.
  DecayBranch(DecayBranch &&rhs)
      : ProcessBranch(rhs.branch_weight_), type_(rhs.type_) {}
  /// \return the quantized angular momentum of this branch.
  inline int angular_momentum() const { return type_.angular_momentum(); }
  const ParticleTypePtrList &particle_types() const override {
    return type_.particle_types();
  }
  unsigned int particle_number() const override {
    return type_.particle_number();
  }
  /// \return the DecayType of the branch.
  inline const DecayType &type() const { return type_; }
  /// \return "Decay" to the process type.
  inline ProcessType get_type() const override { return ProcessType::Decay; }

 private:
  /// Decay type (including final-state particles and angular momentum)
  const DecayType &type_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_PROCESSBRANCH_H_
