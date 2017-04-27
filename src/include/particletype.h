/*
 *    Copyright (c) 2012-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_PARTICLETYPE_H_
#define SRC_INCLUDE_PARTICLETYPE_H_

#include <cassert>
#include <string>
#include <utility>
#include <vector>

#include "forwarddeclarations.h"
#include "macros.h"
#include "pdgcode.h"

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
   * Decay width cutoff for considering a particle as stable.
   *
   * We currently regard a particle type as stable if its on-shell width is less
   * than 200 keV. The cutoff is chosen such that the η and the η' are stable.
   */
  static constexpr float width_cutoff = 1e-5f;

  /**
   * Creates a fully initialized ParticleType object.
   *
   * \param n The name of the particle.
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

  /// Returns the name of the particle.
  const std::string &name() const { return name_; }

  /// Returns the particle pole mass.
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

  /// \copydoc PdgCode::antiparticle_sign
  int antiparticle_sign() const { return pdgcode_.antiparticle_sign(); }

  /** Returns twice the isospin vector length \f$I\f$.
   *
   * This returns e.g. 1 for nucleons, 2 for pions and 3 for Deltas.
   * It is always positive.
   */
  int isospin() const;

  /// \copydoc PdgCode::isospin3
  int isospin3() const { return I3_; }

  /// Returns the isospin-3 component relative to the total isospin.
  float isospin3_rel() const {
    unsigned int I = isospin();
    return (I == 0) ? 0 : static_cast<float>(isospin3())/I;
  }

  /**
   * Returns a pointer to the Isospin-multiplet of this PDG Code.
   */
  IsoParticleType* iso_multiplet() const {
    return iso_multiplet_;
  }

  /// \copydoc PdgCode::charge
  int charge() const { return charge_; }

  /// \copydoc PdgCode::spin
  unsigned int spin() const { return pdgcode_.spin(); }

  /// \copydoc PdgCode::is_hadron
  bool is_hadron() const { return pdgcode_.is_hadron(); }

  /// \copydoc PdgCode::is_lepton
  bool is_lepton() const { return pdgcode_.is_lepton(); }

  /// \copydoc PdgCode::is_baryon
  bool is_baryon() const { return pdgcode_.is_baryon(); }

  /// \copydoc PdgCode::is_meson
  bool is_meson() const { return pdgcode_.is_meson(); }

  /// \copydoc PdgCode::baryon_number
  int baryon_number() const { return pdgcode_.baryon_number(); }

  /// \copydoc PdgCode::strangeness
  int strangeness() const {return pdgcode_.strangeness(); }

  /// \copydoc PdgCode::is_nucleon
  bool is_nucleon() const { return pdgcode_.is_nucleon(); }

  /// \copydoc PdgCode::is_Delta
  bool is_Delta() const { return pdgcode_.is_Delta(); }

  /// Is this a nucleon resonance (N*)?
  inline bool is_Nstar() const {
    return is_baryon() && isospin() == 1 && !pdgcode_.is_nucleon() &&
           pdgcode_.strangeness() == 0 && pdgcode_.charmness() == 0;
  }

  /// Is this a Delta resonance (Delta*)?
  inline bool is_Deltastar() const {
    return is_baryon() && isospin() == 3 && !pdgcode_.is_Delta() &&
           pdgcode_.strangeness() == 0 && pdgcode_.charmness() == 0;
  }

  /// Check if the particle is stable
  inline bool is_stable() const {
    return width_ < width_cutoff;
  }

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
  float partial_width(const float m, const DecayBranch *mode) const;

  /**
   * Get the mass-dependent total width of a particle with mass m.
   *
   * \param m Invariant mass of the decaying particle.
   */
  float total_width(const float m) const;

  /**
   * Get the mass-dependent partial decay widths of a particle with mass m.
   * Returns a list of process branches, whose weights correspond to the
   * actual partial widths. The list contains all branches.
   *
   * \param m Invariant mass of the decaying particle.
   */
  DecayBranchList get_partial_widths(const float m) const;

  /**
  * Get the mass-dependent partial decay widths of a particle with mass m.
  * Returns a list of process branches, whose weights correspond to the
  * actual partial widths. The list contains all but the dilepton branches.
  *
  * \param m Invariant mass of the decaying particle.
  */
  DecayBranchList get_partial_widths_hadronic(const float m) const;

  /**
  * Get the mass-dependent partial decay widths of a particle with mass m.
  * Returns a list of process branches, whose weights correspond to the
  * actual partial widths. The list contains only the dilepton branches.
  *
  * \param m Invariant mass of the decaying particle.
  */
  DecayBranchList get_partial_widths_dilepton(const float m) const;

  /**
   * Get the mass-dependent partial width of a resonance with mass m,
   * decaying into two given daughter particles.
   *
   * \param m Invariant mass of the decaying resonance.
   * \param t_a Type of first daughter particle.
   * \param t_b Type of second daughter particle.
   */
  float get_partial_width(const float m, const ParticleType &t_a,
                                         const ParticleType &t_b) const;

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
   * Full spectral function
   * \f$ A(m) = \frac{2}{\pi} N \frac{m^2\Gamma(m)}{(m^2-m_0^2)^2+(m\Gamma(m))^2} \f$
   * of the resonance (relativistic Breit-Wigner distribution with
   * mass-dependent width, where N is a normalization factor).
   * \param m Actual off-shell mass of the resonance, where the
   *          spectral function is supposed to be evaluated.
   * \note The normalization factor N ensures that the spectral function is
   *       normalized to unity.
   */
  float spectral_function(float m) const;

  /**
   * Full spectral function without normalization factor. */
  float spectral_function_no_norm(float m) const;

  /**
   * The spectral function with a constant width (= width at pole).
   * It is guaranteed to be normalized to 1, when integrated from 0 to inf. */
  float spectral_function_const_width(float m) const;

  /**
   * This one is the most simple form of the spectral function, using a
   * Cauchy distribution (non-relativistic Breit-Wigner with constant width).
   * It can be integrated analytically, and is normalized to 1 when integrated
   * from -inf to inf.
   */
  float spectral_function_simple(float m) const;

  /**
  * Resonance mass sampling for 2-particle final state with one resonance
  * (type given by 'this') and one stable particle.
  *
  * \param[in] mass_stable Mass of the stable particle.
  * \param[in] cms_energy center-of-mass energy of the 2-particle final state.
  * \param[in] L relative angular momentum of the final-state particles
  *
  * \return The mass of the resonance particle.
  */
  float sample_resonance_mass(const float mass_stable,
                              const float cms_energy, int L = 0) const;

  /**
  * Resonance mass sampling for 2-particle final state with two resonances.
  *
  * \param[in] t2 Type of the second resonance
  *               (the first resonance is given by 'this').
  * \param[in] cms_energy center-of-mass energy of the 2-particle final state.
  * \param[in] L relative angular momentum of the final-state particles
  *
  * \return The masses of the resonance particles.
  */
  std::pair<float, float> sample_resonance_masses(const ParticleType &t2,
                                                  const float cms_energy,
                                                  int L = 0) const;

  /**
   *  Prints out width and spectral function versus mass to the
   *  standard output. This is useful for debugging and analysis.
   */
  void dump_width_and_spectral_function() const;

  /**
   * Returns a list of all ParticleType objects.
   *
   * \note The order of the list may be sorted by PDG codes, but do not rely on
   *       this.
   */
  static const ParticleTypeList &list_all();

  /** Returns a list of all nucleons (i.e. proton and neutron). */
  static ParticleTypePtrList &list_nucleons();
  /** Returns a list of all anti-nucleons (i.e. anti-proton and anti-neutron).
    */
  static ParticleTypePtrList &list_anti_nucleons();
  /** Returns a list of the Delta(1232) baryons // oliiny: only 1232?!
   *  (i.e. all four charge states). */
  static ParticleTypePtrList &list_Deltas();
  /** Returns a list of the anti-Delta(1232) baryons
   *  (i.e. all four charge states). */
  static ParticleTypePtrList &list_anti_Deltas();
  /** Returns a list of all baryon resonances,
   * i.e. unstable baryons (not including antibaryons). */
  static ParticleTypePtrList &list_baryon_resonances();

  /**
   * Returns the ParticleTypePtr for the given \p pdgcode.
   * If the particle type is not found, an invalid ParticleTypePtr is returned.
   * You can convert a ParticleTypePtr to a bool to check whether it is valid.
   *
   * \note The complexity of the search is \f$\mathcal O(\log N)\f$. Therefore,
   * do not use this function except for user input that selects a particle
   * type. All other internal references for a particle type should use
   * ParticleTypePtr instead.
   */
  static const ParticleTypePtr try_find(PdgCode pdgcode);

  /**
   * Returns the ParticleType object for the given \p pdgcode.
   * If the particle is not found, a PdgNotFoundFailure is thrown.
   *
   * \note The complexity of the search is \f$\mathcal O(\log N)\f$. Therefore,
   * do not use this function except for user input that selects a particle
   * type. All other internal references for a particle type should use
   * ParticleTypePtr instead.
   */
  static const ParticleType &find(PdgCode pdgcode);
  /// \ingroup exception
  struct PdgNotFoundFailure : public std::runtime_error {
    using std::runtime_error::runtime_error;
  };


  /**
   * Returns whether the ParticleType with the given \p pdgcode exists.
   *
   * \note The complexity of the search is \f$\mathcal O(\log N)\f$.
   */
  static bool exists(PdgCode pdgcode);

  /**
   * Returns whether the ParticleType with the given \p name exists.
   *
   * \note The complexity of the search is \f$\mathcal O(N)\f$.
   */
  static bool exists(const std::string& name);

  /**
   * Initialize the global ParticleType list (list_all) from the given input
   * data. This function must only be called once (will fail on second
   * invocation).
   *
   * \param particles A string that contains the definition of ParticleTypes to
   *                  be created.
   */
  static void create_type_list(const std::string &particles);

  /// Returns whether the two ParticleType objects have the same PDG code.
  bool operator==(const ParticleType &rhs) const {
    return pdgcode() == rhs.pdgcode();
  }
  /// Returns whether the two ParticleType objects have different PDG codes.
  bool operator!=(const ParticleType &rhs) const {
    return pdgcode() != rhs.pdgcode();
  }
  /// "Less than" operator for sorting the ParticleType list (by PDG code)
  bool operator<(const ParticleType &rhs) const {
    return pdgcode() < rhs.pdgcode();
  }

  /**
   * Check if unstable particles have any decay modes and throw errors.
   */
  static void check_consistency();

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

  /// \ingroup exception
  struct LoadFailure : public std::runtime_error {
    using std::runtime_error::runtime_error;
  };

 private:
  /// name of the particle
  std::string name_;
  /// pole mass of the particle
  float mass_;
  /// width of the particle
  float width_;
  /// PDG Code of the particle
  PdgCode pdgcode_;
  /// minimum mass of the particle
  /* Mutable, because it is initialized at first call of minimum mass function,
     so it's logically const, but not physically const, which is a classical
     case for using mutable. */
  mutable float minimum_mass_;
  /** This normalization factor ensures that the spectral function is normalized
   * to unity, when integrated over its full domain. */
  mutable float norm_factor_ = -1.;
  /** charge, isospin and isospin projection of the particle
   *
   * This is filled automatically from pdgcode_.
   */
  int charge_;
  mutable int isospin_;
  int I3_;

  IsoParticleType *iso_multiplet_ = nullptr;

  // Maximum factor for single-res mass sampling, cf. sample_resonance_mass.
  mutable float max_factor1_ = 1.;
  // Maximum factor for double-res mass sampling, cf. sample_resonance_masses.
  mutable float max_factor2_ = 1.;

  /**\ingroup logging
   * Writes all information about the particle type to the output stream.
   */
  friend std::ostream &operator<<(std::ostream &out, const ParticleType &type);
};

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

  /// Returns whether the two objects reference the same ParticleType object.
  bool operator==(const ParticleTypePtr &rhs) const {
    return index_ == rhs.index_;
  }

  /// Returns whether the two objects reference different ParticleType objects.
  bool operator!=(const ParticleTypePtr &rhs) const {
    return index_ != rhs.index_;
  }
  // "Less than" operator
  bool operator<(const ParticleTypePtr &rhs) const {
    return index_ < rhs.index_;
  }

  /// Returns whether the objects stores a valid ParticleType reference.
  operator bool() const { return index_ != 0xffff; }

 private:
  /// ParticleType::operator& is a friend in order to call the constructor
  friend ParticleTypePtr ParticleType::operator&() const;

  /// Constructs a pointer to the ParticleType object at offset \p i.
  explicit ParticleTypePtr(std::uint16_t i) : index_(i) {}

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
  std::uint16_t index_ = 0xffff;
};

// #define SMASH_INLINE_LIST_ALL 1
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
