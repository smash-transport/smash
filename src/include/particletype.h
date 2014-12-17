/*
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_PARTICLETYPE_H_
#define SRC_INCLUDE_PARTICLETYPE_H_

#include "forwarddeclarations.h"
#include "macros.h"
#include "pdgcode.h"

#include <assert.h>
#include <string>
#include <vector>

namespace Smash {

/**
 * \ingroup data
 *
 * Particle type contains the static properties of a particle
 *
 * Before creation of Experiment, SMASH initializes the list of particles
 * (\ref list_all). After construction these values are immutable.
 *
 * The list of particles is stored in such a way that look up of a ParticleType
 * object (\ref find) for a given PDG code is as efficient as possible
 * (\f$\mathcal O(\log N)\f$). This is still not efficient enough to use PdgCode
 * as a substitute for storing information about a particle type, though. Use
 * ParticleTypePtr instead.
 */
class ParticleType {
 public:
  /**
   * Creates a fully initialized ParticleType object.
   *
   * \param n The name of the particle (only used for debug output).
   * \param m The mass of the particle.
   * \param w The width of the particle.
   * \param id The PDG code of the particle.
   *
   * \note The remaining properties ParticleType provides are derived from the
   *       PDG code and therefore cannot be set explicitly (this avoids the
   *       chance of introducing inconsistencies).
   */
  ParticleType(std::string n, float m, float w, PdgCode id);

  /**
   * Copies are not allowed as they break intended use. Instead use a const-ref
   * or ParticleTypePtr (as returned from operator&).
   */
  ParticleType(const ParticleType &) = delete;
  /// assignment is not allowed, see copy constructor above
  ParticleType &operator=(const ParticleType &) = delete;

  // move ctors are needed for std::sort
  ParticleType(ParticleType &&) = default;
  ParticleType &operator=(ParticleType &&) = default;

  /// Returns the DecayModes object for this particle type.
  const DecayModes &decay_modes() const;

  /// Returns the name of the particle (for debug output only).
#ifdef NDEBUG
  std::string name() const { return {}; }
#else
  const std::string &name() const { return name_; }
#endif

  /// Returns the particle mass.
  float mass() const { return mass_; }

  /// Returns the squared particle mass.
  float mass_sqr() const { return mass_*mass_; }

  /// Returns the particle width (at the mass pole).
  float width_at_pole() const { return width_; }

  /// Returns the PDG code of the particle.
  PdgCode pdgcode() const { return pdgcode_; }

  /// \copydoc PdgCode::has_antiparticle
  bool has_antiparticle() const { return pdgcode_.has_antiparticle(); }

  /// Return a pointer to the corresponding antiparticle ParticleType object.
  ParticleTypePtr get_antiparticle() const;

  /// \copydoc PdgCode::isospin_total
  unsigned int isospin() const { return isospin_; }

  /// \copydoc PdgCode::isospin3
  int isospin3() const { return pdgcode_.isospin3(); }

  /// \copydoc PdgCode::charge
  int charge() const { return charge_; }

  /// \copydoc PdgCode::spin
  unsigned int spin() const { return pdgcode_.spin(); }

  /// \copydoc PdgCode::is_hadron
  bool is_hadron() const { return pdgcode_.is_hadron(); }

  /// \copydoc PdgCode::baryon_number
  int baryon_number() const { return pdgcode_.baryon_number(); }

  /// Check if the particle is stable
  inline bool is_stable() const;

  /**
   * The minimum mass of the resonance.
   *
   * Calculate the minimum rest energy the resonance must have
   * for any decay channel to be kinematically available.
   * (In other words, find the smallest sum of final-state particle masses.)
   *
   * \return The minimum mass that a particle of this type can assume.
   */
  float minimum_mass() const;

  /**
   * Get the mass-dependent partial decay width of a particle with mass m
   * in a particular decay mode.
   *
   * \param m Invariant mass of the decaying particle.
   * \param mode Decay mode to consider.
   */
  float partial_width(const float m, const DecayBranch &mode) const;

  /**
   * Get the mass-dependent total width of a particle with mass m.
   *
   * \param m Invariant mass of the decaying particle.
   */
  float total_width(const float m) const;

  /**
   * Get the mass-dependent partial decay widths of a particle with mass m.
   * Returns a list of process branches, whose weights correspond to the
   * actual partial widths.
   *
   * \param m Invariant mass of the decaying particle.
   */
  ProcessBranchList get_partial_widths(const float m) const;

  /**
   * Get the mass-dependent partial in-width of a resonance with mass m,
   * decaying into two given daughter particles. For stable daughter
   * particles, the in-width equals the 'normal' partial decay width
   * (i.e. the 'out-width').
   *
   * \param m Invariant mass of the decaying resonance.
   * \param p_a First daughter particle.
   * \param p_b Second daughter particle.
   */
  float get_partial_in_width(const float m, const ParticleData &p_a,
                                            const ParticleData &p_b) const;

  /**
   * Returns a list of all ParticleType objects.
   *
   * \note The order of the list may be sorted by PDG codes, but do not rely on
   *       this.
   */
  static const ParticleTypeList &list_all();

  /** Returns a list of all nucleons (i.e. proton and neutron). */
  static std::vector<ParticleTypePtr> list_nucleons();
  /** Returns a list of all baryon resonances, i.e. unstable baryons (not including antibaryons). */
  static std::vector<ParticleTypePtr> list_baryon_resonances();

  /**
   * Returns the ParticleType object for the given \p pdgcode.
   *
   * \note The complexity of the search is \f$\mathcal O(\log N)\f$. Therefore,
   * do not use this function except for user input that selects a particle
   * type. All other internal references for a particle type should use
   * ParticleTypePtr instead.
   */
  static const ParticleType &find(PdgCode pdgcode) SMASH_CONST;

  /**
   * Returns whether the ParticleType with the given \p pdgcode exists.
   *
   * \note The complexity of the search is \f$\mathcal O(\log N)\f$.
   */
  static bool exists(PdgCode pdgcode) SMASH_CONST;

  /**
   * Initialize the global ParticleType list (list_all) from the given input
   * data.
   *
   * \param particles A string that contains the definition of ParticleTypes to
   *                  be created.
   */
  static void create_type_list(const std::string &particles);

  /**
   * Returns an object that acts like a pointer, except that it requires only 2
   * bytes and inhibits pointer arithmetics.
   *
   * This is an optimization for creating references to ParticleType objects.
   * With a normal pointer you would require the original object to stay in its
   * place in memory, as otherwise the pointer would dangle. The ParticleTypePtr
   * does not have this problem as it only stores the index of the ParticleType
   * in the global vector of ParticleType objects.
   *
   * In addition to returning a more efficient reference type, the overload of
   * operator& effectively inhibits passing ParticleType objects by pointer.
   * This is an intended restriction since ParticleType objects should only be
   * passed by const-ref. (You can now pass by ParticleTypePtr instead, but
   * please prefer not to. It might be useful when you want to store a reference
   * inside the function anyway, but that can just as well be done with a
   * const-ref parameter. A ParticleTypePtr can be invalid, a const-ref
   * is always a valid reference semantically.)
   *
   * \par Pre-condition:
   * The operator expects that the ParticleType object is stored in the vector
   * returned by ParticleType::list_all. Therefore, never create new ParticleType
   * objects (that includes copies and moves)!
   *
   * \par Note on distributed execution:
   * At some point we might want to have several copies of the ParticleType
   * vector - on different machines or NUMA nodes. In that case
   * ParticleType::list_all will return the local vector. This operator will
   * continue to work as expected as long as the ParticleType object is an entry
   * of this local vector. The ParticleTypePtr can then be used to communicate a
   * ParticleType over node / NUMA boundaries (an actual pointer would not work,
   * though).
   *
   * \return A pointer-like object referencing this ParticleType object.
   *
   * \see ParticleTypePtr
   */
  ParticleTypePtr operator&() const;

 private:
#ifndef NDEBUG
  /// name of the particle
  /// This variable is only used for debug output. Non-debug builds save the
  /// memory to be more cache-efficient.
  std::string name_;
#endif
  /// mass of the particle
  float mass_;
  /// width of the particle
  float width_;
  /// PDG Code of the particle
  PdgCode pdgcode_;
  /** twice the isospin of the particle
   *
   * This is filled automatically from pdgcode_.
   */
  unsigned int isospin_;
  /** charge of the particle
   *
   * This is filled automatically from pdgcode_.
   */
  int charge_;

  /**\ingroup logging
   * Writes all information about the particle type to the output stream.
   */
  friend std::ostream &operator<<(std::ostream &out, const ParticleType &type);
};

inline bool ParticleType::is_stable() const {
  /* We currently regard a particle type as stable if its on-shell width is
   * less than 10 keV. */
  return width_ < 1E-5f;
}

/**
 * \ingroup data
 *
 * A pointer-like interface to global references to ParticleType objects.
 *
 * \see ParticleType::operator&
 */
class ParticleTypePtr {
 public:
  /// Dereferences the pointer and returns the ParticleType object.
  const ParticleType &operator*() const {
    return lookup();
  }

  /// Dereferences the pointer and returns the ParticleType object.
  const ParticleType *operator->() const {
    // this requires std::addressof because &lookup() would call
    // ParticleType::operator& and return ParticleTypePtr again
    return std::addressof(lookup());
  }

  /// Default construction initializes with an invalid index.
  ParticleTypePtr() = default;

  /// Initialization with \c nullptr constructs an object with an invalid index.
  ParticleTypePtr(std::nullptr_t) {}

  /// Returns whether the two objects reference the same ParticleType object.
  bool operator==(const ParticleTypePtr &rhs) const {
    return index_ == rhs.index_;
  }

  /// Returns whether the two objects reference different ParticleType objects.
  bool operator!=(const ParticleTypePtr &rhs) const {
    return index_ != rhs.index_;
  }

  /// Returns whether the objects stores a valid ParticleType reference.
  operator bool() const { return index_ != 0xffff; }

 private:
  /// ParticleType::operator& is a friend in order to call the constructor
  friend ParticleTypePtr ParticleType::operator&() const;

  /// Constructs a pointer to the ParticleType object at offset \p i.
  ParticleTypePtr(std::uint16_t i) : index_(i) {}

  /**
   * Helper function that does the ParticleType lookup from the stored index.
   *
   * \par Implementation:
   * In debug builds this function asserts that the index is valid.
   * It then asks for the vector of all ParticleType objects
   * (ParticleType::list_all) and uses vector::operator[] to return an lvalue
   * reference to the ParticleType object at the offset \p index_.
   */
  const ParticleType &lookup() const {
    assert(index_ != 0xffff);
    return ParticleType::list_all()[index_];
  }

  /**
   * Stores the index of the references ParticleType object in the global
   * vector. The value 0xffff is used to denote an invalid index (similar to a
   * null pointer).
   */
  std::uint16_t index_= 0xffff;
};

//#define SMASH_INLINE_LIST_ALL 1
#ifdef SMASH_INLINE_LIST_ALL
extern const ParticleTypeList *all_particle_types;
inline const ParticleTypeList &ParticleType::list_all() {
  assert(all_particle_types);
  return *all_particle_types;
}
inline ParticleTypePtr ParticleType::operator&() const {
  const auto offset = this - std::addressof(list_all()[0]);
  assert(offset >= 0 && offset < 0xffff);
  return {static_cast<uint16_t>(offset)};
}
#endif

inline ParticleTypePtr ParticleType::get_antiparticle() const {
  assert(has_antiparticle());
  return &find(pdgcode_.get_antiparticle());
}

}  // namespace Smash

#endif  // SRC_INCLUDE_PARTICLETYPE_H_
