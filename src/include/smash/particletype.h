/*
 *    Copyright (c) 2012-2018
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

namespace smash {

/**
 * Represent the parity of a particle type.
 */
enum class Parity {
  /// Positive parity.
  Pos,
  /// Negative parity.
  Neg
};

/**
 * \param p Given parity.
 * \return Inverted parity.
 */
inline Parity invert_parity(Parity p) {
  switch (p) {
    case Parity::Pos:
      return Parity::Neg;
    case Parity::Neg:
      return Parity::Pos;
  }
  // This is unreachable and should be optimized away.
  // It is required to silence a compiler warning.
  throw std::runtime_error("unreachable");
}

/**
 * \ingroup data
 *
 * Particle type contains the static properties of a particle species
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
   * than 10 keV. The cutoff is chosen such that the η and the η' are stable.
   */
  // If this is changed, make sure to update the userguide in
  // `include/configuration.h`.
  static constexpr double width_cutoff = 1e-5;

  /**
   * Creates a fully initialized ParticleType object.
   *
   * \param[in] n The name of the particle.
   * \param[in] m The mass of the particle.
   * \param[in] w The width of the particle.
   * \param[in] p The parity of the particle.
   * \param[in] id The PDG code of the particle.
   *
   * \note The remaining properties ParticleType provides are derived from the
   *       PDG code and therefore cannot be set explicitly (this avoids the
   *       chance of introducing inconsistencies).
   */
  ParticleType(std::string n, double m, double w, Parity p, PdgCode id);

  /**
   * Copies are not allowed as they break intended use. Instead use a const-ref
   * or ParticleTypePtr (as returned from operator&).
   */
  ParticleType(const ParticleType &) = delete;
  /// assignment is not allowed, see copy constructor above
  ParticleType &operator=(const ParticleType &) = delete;

  /// move ctors are needed for std::sort
  ParticleType(ParticleType &&) = default;
  /// move ctors are needed for std::sort
  ParticleType &operator=(ParticleType &&) = default;

  /// \return the DecayModes object for this particle type.
  const DecayModes &decay_modes() const;

  /// \return the name of the particle.
  const std::string &name() const { return name_; }

  /// \return the particle pole mass.
  double mass() const { return mass_; }

  /// \return the squared particle mass.
  double mass_sqr() const { return mass_ * mass_; }

  /// \return the particle width (at the mass pole).
  double width_at_pole() const { return width_; }

  /// \return the parity of the particle.
  Parity parity() const { return parity_; }

  /// \return the PDG code of the particle.
  PdgCode pdgcode() const { return pdgcode_; }

  /// \copydoc PdgCode::has_antiparticle
  bool has_antiparticle() const { return pdgcode_.has_antiparticle(); }

  /// \return a pointer to the corresponding antiparticle ParticleType object.
  ParticleTypePtr get_antiparticle() const;

  /// \copydoc PdgCode::antiparticle_sign
  int antiparticle_sign() const { return pdgcode_.antiparticle_sign(); }

  /**
   * Returns twice the isospin vector length \f$I\f$.
   *
   * This returns e.g. 1 for nucleons, 2 for pions and 3 for Deltas.
   * It is always positive.
   */
  int isospin() const;

  /// \copydoc PdgCode::isospin3
  int isospin3() const { return I3_; }

  /// \return the isospin-3 component relative to the total isospin.
  double isospin3_rel() const {
    unsigned int I = isospin();
    return (I == 0) ? 0 : static_cast<double>(isospin3()) / I;
  }

  /// \return a pointer to the Isospin-multiplet of this PDG Code.
  IsoParticleType *iso_multiplet() const { return iso_multiplet_; }

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
  int strangeness() const { return pdgcode_.strangeness(); }

  /// \copydoc PdgCode::is_nucleon
  bool is_nucleon() const { return pdgcode_.is_nucleon(); }

  /// \copydoc PdgCode::is_Delta
  bool is_Delta() const { return pdgcode_.is_Delta(); }

  /// \copydoc PdgCode::is_rho
  bool is_rho() const { return pdgcode_.is_rho(); }

  /// \return Is this a nucleon resonance (N*)?
  inline bool is_Nstar() const {
    return is_baryon() && isospin() == 1 && !pdgcode_.is_nucleon() &&
           pdgcode_.strangeness() == 0 && pdgcode_.charmness() == 0;
  }

  /// \copydoc PdgCode::is_Nstar1535
  bool is_Nstar1535() const { return pdgcode_.is_Nstar1535(); }

  /// \return Is this a Delta resonance (Delta*)?
  inline bool is_Deltastar() const {
    return is_baryon() && isospin() == 3 && !pdgcode_.is_Delta() &&
           pdgcode_.strangeness() == 0 && pdgcode_.charmness() == 0;
  }

  /// \return whether the particle is stable
  inline bool is_stable() const { return width_ < width_cutoff; }

  /// \return whether the particle is a nucleus
  inline bool is_nucleus() const { return pdgcode_.is_nucleus(); }

  /// \return whether the particle is an (anti-)deuteron
  inline bool is_deuteron() const {
    return is_nucleus() && std::abs(pdgcode_.get_decimal()) == 1000010020;
  }

  /// \return whether the particle is an artificial d' resonance
  inline bool is_dprime() const {
    return is_nucleus() && std::abs(pdgcode_.get_decimal()) == 1000010021;
  }

  /**
   * The minimum mass of the resonance that is kinematically allowed.
   *
   * Calculate the minimum rest energy the resonance must have
   * for any decay channel to be kinematically available.
   * (In other words, find the smallest sum of final-state particle masses.)
   *
   * \return The minimum mass that a particle of this type can assume, where at
   *         least one decay is possible.
   */
  double min_mass_kinematic() const;

  /**
   * The minimum mass of the resonance, where the spectral function is non-zero.
   *
   * Calculate the the smallest mass where the spectral function still has a
   * contribution. This value can be different from min_mass_kinematic,
   * if the spectral function becomes zero at masses higher than
   * min_mass_kinematic.
   *
   * \return The minimum mass that a particle of this type can assume, where the
   * spectral function still has a non-zero value.
   */
  double min_mass_spectral() const;

  /**
   * Get the mass-dependent partial decay width of a particle with mass m
   * in a particular decay mode.
   *
   * \param[in] m Invariant mass of the decaying particle.
   * \param[in] mode Decay mode to consider.
   * \return the partial width of this specific mode for this mass
   */
  double partial_width(const double m, const DecayBranch *mode) const;

  /**
   * Get the mass-dependent total width of a particle with mass m.
   *
   * \param[in] m Invariant mass of the decaying particle.
   * \return the total width for all modes for this mass
   */
  double total_width(const double m) const;

  /**
   * Get all the mass-dependent partial decay widths of a particle with mass m.
   * \todo lots of code duplication in general in these partial width functions
   *
   * \param[in] m Invariant mass of the decaying particle.
   * \return a list of process branches, whose weights correspond to the
   *         actual partial widths. The list contains all branches.
   */
  DecayBranchList get_partial_widths(const double m) const;

  /**
   * Get the mass-dependent hadronic partial decay widths
   * of a particle with mass m.
   * \todo lots of code duplication with get_partial_widths_dilepton
   *
   * \param[in] m Invariant mass of the decaying particle.
   * \return a list of process branches, whose weights correspond to the
   *         actual partial widths. The list contains only hadronic branches.
   * \throw runtime_error if a decay has less than 2 or more than 3 products
   */
  DecayBranchList get_partial_widths_hadronic(const double m) const;

  /**
   * Get the mass-dependent dilepton partial decay widths
   * of a particle with mass m.
   * \todo lots of code duplication with get_partial_widths_hadronic
   *
   * \param[in] m Invariant mass of the decaying particle.
   * \return a list of process branches, whose weights correspond to the
   *         actual partial widths. The list contains only dilepton branches.
   * \throw runtime_error if a decay has less than 2 or more than 3 products
   */
  DecayBranchList get_partial_widths_dilepton(const double m) const;

  /**
   * Get the mass-dependent partial width of a resonance with mass m,
   * decaying into two given daughter particles.
   *
   * \param[in] m Invariant mass of the decaying resonance.
   * \param[in] t_a Type of first daughter particle.
   * \param[in] t_b Type of second daughter particle.
   * \return the partial width for this mass and this specific decay channel
   */
  double get_partial_width(const double m, const ParticleType &t_a,
                           const ParticleType &t_b) const;

  /**
   * Get the mass-dependent partial in-width of a resonance with mass m,
   * decaying into two given daughter particles. For stable daughter
   * particles, the in-width equals the 'normal' partial decay width
   * (i.e. the 'out-width').
   *
   * \param[in] m Invariant mass of the decaying resonance.
   * \param[in] p_a First daughter particle.
   * \param[in] p_b Second daughter particle.
   * \return the partial in-width for this mass and this specific decay channel
   */
  double get_partial_in_width(const double m, const ParticleData &p_a,
                              const ParticleData &p_b) const;

  /**
   * Full spectral function
   * \f$ A(m) = \frac{2}{\pi} N
   * \frac{m^2\Gamma(m)}{(m^2-m_0^2)^2+(m\Gamma(m))^2} \f$
   * of the resonance (relativistic Breit-Wigner distribution with
   * mass-dependent width, where N is a normalization factor).
   *
   * \param[in] m Actual off-shell mass of the resonance, where the
   *              spectral function is to be evaluated.
   * \return the value of the spectral function for this mass
   *
   * \note The normalization factor N ensures that the spectral function is
   *       normalized to unity.
   */
  double spectral_function(double m) const;

  /**
   * Full spectral function without normalization factor.
   * \see spectral_function
   *
   * \param[in] m Actual off-shell mass of the resonance, where the
   *              spectral function is to be evaluated.
   * \return the value of the non-normalized spectral function for this mass
   */
  double spectral_function_no_norm(double m) const;

  /**
   * \todo unused
   * The spectral function with a constant width (= width at pole).
   * It is guaranteed to be normalized to 1, when integrated from 0 to inf.
   */
  double spectral_function_const_width(double m) const;

  /**
   * This one is the most simple form of the spectral function, using a
   * Cauchy distribution (non-relativistic Breit-Wigner with constant width).
   * It can be integrated analytically, and is normalized to 1 when integrated
   * from -inf to inf.
   *
   * \param[in] m Actual off-shell mass of the resonance, where the
   *              spectral function is to be evaluated.
   * \return the Cauchy spectral function at mass m
   */
  double spectral_function_simple(double m) const;

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
  double sample_resonance_mass(const double mass_stable,
                               const double cms_energy, int L = 0) const;

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
  std::pair<double, double> sample_resonance_masses(const ParticleType &t2,
                                                    const double cms_energy,
                                                    int L = 0) const;

  /**
   * Prints out width and spectral function versus mass to the
   * standard output. This is useful for debugging and analysis.
   *
   * \throw if the particle type is stable
   */
  void dump_width_and_spectral_function() const;

  /**
   * \return a list of all ParticleType objects.
   *
   * \note This list is currently sorted, but do not rely on it.
   */
  static const ParticleTypeList &list_all();

  /// \return a list of all nucleons (i.e. proton and neutron).
  static ParticleTypePtrList &list_nucleons();
  /// \return a list of all anti-nucleons (i.e. anti-proton and anti-neutron).
  static ParticleTypePtrList &list_anti_nucleons();
  /**
   * \return a list of the Delta(1232) baryons
   *         (i.e. all four charge states).
   */
  static ParticleTypePtrList &list_Deltas();
  /**
   * \return a list of the anti-Delta(1232) baryons
   *         (i.e. all four charge states).
   */
  static ParticleTypePtrList &list_anti_Deltas();
  /**
   * \return a list of all baryon resonances,
   *         i.e. unstable baryons (not including antibaryons).
   */
  static ParticleTypePtrList &list_baryon_resonances();
  /**
   * \return a list of all light nuclei from SMASH particle list.
   *         Nucleons are not included into light nuclei by convention.
   */
  static ParticleTypePtrList &list_light_nuclei();

  /**
   * Returns the ParticleTypePtr for the given \p pdgcode.
   * If the particle type is not found, an invalid ParticleTypePtr is returned.
   * You can convert a ParticleTypePtr to a bool to check whether it is valid.
   *
   * \note The complexity of the search is \f$\mathcal O(\log N)\f$. Therefore,
   * do not use this function except for user input that selects a particle
   * type. All other internal references for a particle type should use
   * ParticleTypePtr instead.
   *
   * \param[in] pdgcode the unique pdg code to try to find
   * \return the ParticleTypePtr that corresponds to this pdg code, or an
             invalid pointer
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
   *
   * \param[in] pdgcode the unique pdg code to try to find
   * \return the ParticleTypePtr that corresponds to this pdg code
   * \throw PdgNotFoundFailure pdgcode not found in available particle types
   */
  static const ParticleType &find(PdgCode pdgcode);

  /// \ingroup exception
  struct PdgNotFoundFailure : public std::runtime_error {
    using std::runtime_error::runtime_error;
  };

  /**
   * \param[in] pdgcode the PdgCode to look for
   * \return whether the ParticleType with the given \p pdgcode exists.
   *
   * \note The complexity of the search is \f$\mathcal O(\log N)\f$.
   */
  static bool exists(PdgCode pdgcode);

  /**
   * \param[in] name the name to look for
   * \return whether the ParticleType with the given \p name exists.
   *
   * \note The complexity of the search is \f$\mathcal O(N)\f$.
   */
  static bool exists(const std::string &name);

  /**
   * Initialize the global ParticleType list (list_all) from the given input
   * data. This function must only be called once (will fail on second
   * invocation).
   *
   * \param[in] particles A string that contains the definition of ParticleTypes
   *                      to be created.
   * \throw LoadFailure if a line in the particle file could not be read, or if
   *                    there are duplicates in it
   * \throw runtime_error if the mass of of nucleons, kaons and deltas are
   *                      different from the hardcoded masses, or if this
   *                      function is called more than once
   */
  static void create_type_list(const std::string &particles);

  /**
   * \param[in] rhs another ParticleType to compare to
   * \return whether the two ParticleType objects have the same PDG code.
   */
  bool operator==(const ParticleType &rhs) const {
    return pdgcode() == rhs.pdgcode();
  }
  /**
   * \param[in] rhs another ParticleType to compare to
   * \return whether the two ParticleType objects have different PDG codes.
   */
  bool operator!=(const ParticleType &rhs) const {
    return pdgcode() != rhs.pdgcode();
  }
  /**
   * "Less than" operator for sorting the ParticleType list (by PDG code)
   *
   * \param[in] rhs another ParticleType to compare to
   * \return whether the PDG code of rhs is larger than the one from *this
   */
  bool operator<(const ParticleType &rhs) const {
    return pdgcode() < rhs.pdgcode();
  }

  /// \throw runtime_error if unstable particles have no decay modes
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
   * returned by ParticleType::list_all. Therefore, never create new
   * ParticleType
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
  double mass_;
  /// width of the particle
  double width_;
  /// Parity of the particle
  Parity parity_;
  /// PDG Code of the particle
  PdgCode pdgcode_;
  /**
   * minimum kinematically allowed mass of the particle
   * Mutable, because it is initialized at first call of minimum mass function,
   * so it's logically const, but not physically const, which is a classical
   * case for using mutable.
   */
  mutable double min_mass_kinematic_;
  /**
   * minimum mass, where the spectral function is non-zero
   * Mutable, because it is initialized at first call of minimum mass function,
   * so it's logically const, but not physically const, which is a classical
   * case for using mutable.
   */
  mutable double min_mass_spectral_;
  /** This normalization factor ensures that the spectral function is normalized
   * to unity, when integrated over its full domain. */
  mutable double norm_factor_ = -1.;
  /// Charge of the particle; filled automatically from pdgcode_.
  int charge_;
  /// Isospin of the particle; filled automatically from pdgcode_.
  mutable int isospin_;
  /// Isospin projection of the particle; filled automatically from pdgcode_.
  int I3_;

  /// Container for the isospin multiplet information
  IsoParticleType *iso_multiplet_ = nullptr;

  /// Maximum factor for single-res mass sampling, cf. sample_resonance_mass.
  mutable double max_factor1_ = 1.;
  /// Maximum factor for double-res mass sampling, cf. sample_resonance_masses.
  mutable double max_factor2_ = 1.;

  /**\ingroup logging
   * Writes all information about the particle type to the output stream.
   *
   * \param[out] out The ostream into which to output
   * \param[in] type The ParticleType object to write into out
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
  /// \return Dereferences the pointer and returns the ParticleType object.
  const ParticleType &operator*() const { return lookup(); }

  /// \return Dereferences the pointer and returns the ParticleType object.
  const ParticleType *operator->() const {
    // this requires std::addressof because &lookup() would call
    // ParticleType::operator& and return ParticleTypePtr again
    return std::addressof(lookup());
  }

  /// Default construction initializes with an invalid index.
  ParticleTypePtr() = default;

  /**
   * \param[in] rhs the ParticleTypePtr to compare to
   * \return whether the two objects reference the same ParticleType object.
   */
  bool operator==(const ParticleTypePtr &rhs) const {
    return index_ == rhs.index_;
  }

  /**
   * \param[in] rhs the ParticleTypePtr to compare to
   * \return whether the two objects reference different ParticleType objects.
   */
  bool operator!=(const ParticleTypePtr &rhs) const {
    return index_ != rhs.index_;
  }
  /**
   * "Less than" operator
   *
   * \param[in] rhs the ParticleTypePtr to compare to
   * \return whether the index is smaller than rhs' index.
   */
  bool operator<(const ParticleTypePtr &rhs) const {
    return index_ < rhs.index_;
  }

  /// \return whether the objects stores a valid ParticleType reference.
  operator bool() const { return index_ != 0xffff; }

 private:
  /**
   * ParticleType::operator& is a friend in order to call the constructor
   *
   * \return the pointer to a ParticleType
   */
  friend ParticleTypePtr ParticleType::operator&() const;

  /** Constructs a pointer to the ParticleType object at offset \p i.
   *
   * \param[in] i the offset where to create.
   */
  explicit ParticleTypePtr(std::uint16_t i) : index_(i) {}

  /**
   * Helper function that does the ParticleType lookup from the stored index.
   *
   * \par Implementation:
   * In debug builds this function asserts that the index is valid.
   * It then asks for the vector of all ParticleType objects
   * (ParticleType::list_all) and uses vector::operator[] to return an lvalue
   * reference to the ParticleType object at the offset \p index_.
   *
   * \return a reference to the ParticleType object at offset \p index_.
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

// #define some global variables and functions
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

}  // namespace smash

#endif  // SRC_INCLUDE_PARTICLETYPE_H_
