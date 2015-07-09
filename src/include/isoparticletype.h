/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_ISOPARTICLETYPE_H_
#define SRC_INCLUDE_ISOPARTICLETYPE_H_

#include <string>

#include "particletype.h"

namespace Smash {

/**
 * \ingroup data
 *
 * IsoParticleType is a class to represent isospin multiplets.
 * It is similar to ParticleType, but refers to whole multiplets instead of
 * single particle states.
 */
class IsoParticleType {
 public:
  /**
   * Creates a fully initialized IsoParticleType object.
   *
   * \param n The name of the multiplet.
   * \param m The (average) mass of the multiplet.
   * \param w The (average) width of the multiplet.
   * \param i Twice the total isospin of the multiplet.
   */
  IsoParticleType(std::string n, float m, float w, int i);

  /**
   * Copies are not allowed as they break intended use. Instead use a const-ref
   * or ParticleTypePtr (as returned from operator&).
   */
  IsoParticleType(const IsoParticleType &) = delete;
  /// assignment is not allowed, see copy constructor above
  IsoParticleType &operator=(const IsoParticleType &) = delete;

  // move ctors are needed for std::sort
  IsoParticleType(IsoParticleType &&) = default;
  IsoParticleType &operator=(IsoParticleType &&) = default;

  /// Returns the name of the multiplet.
  const std::string &name() const { return name_; }

  /// Returns the (average) multiplet mass.
  float mass() const { return mass_; }

  /// Returns the (average) multiplet width.
  float width() const { return width_; }

  ParticleTypePtrList get_states() const { return states_; }

  /// Add a new state to an existing multiplet.
  void add_state(const ParticleType &type) { states_.push_back(&type); }

  /**
   * Returns the IsoParticleType object for the given \p name.
   *
   * \note The complexity of the search is \f$\mathcal O(\log N)\f$.
   */
  static const IsoParticleType& find(std::string name) SMASH_CONST;

  /**
   * Returns whether the ParticleType with the given \p pdgcode exists.
   *
   * \note The complexity of the search is \f$\mathcal O(\log N)\f$.
   */
  static bool exists(std::string name) SMASH_CONST;

  /**
   * Returns the ParticleType object for the given \p name, by first finding the
   * correct multiplet and then looking for the desired state.
   */
  static const ParticleTypePtr find_state(std::string name) SMASH_CONST;

  /**
   * Add a new multiplet to the global list of IsoParticleTypes, which contains
   * the type t. If the multiplet exists already, the type t will be added to it.
   */
  static void create_multiplet(const ParticleType &type);

 private:
  /// name of the multiplet
  std::string name_;
  /// (average) mass of the multiplet
  float mass_;
  /// (average) width of the multiplet
  float width_;
  /// twice the total isospin of the multiplet
  unsigned int isospin_;
  /// list of states that are contained in the multiplet
  ParticleTypePtrList states_;
  /**
   * Private version of the 'find' method that returns a non-const reference.
   */
  static IsoParticleType& find_private(std::string name);
};

}  // namespace Smash

#endif  // SRC_INCLUDE_ISOPARTICLETYPE_H_
