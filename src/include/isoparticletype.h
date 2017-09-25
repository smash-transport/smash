/*
 *    Copyright (c) 2015-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_ISOPARTICLETYPE_H_
#define SRC_INCLUDE_ISOPARTICLETYPE_H_

#include <string>
#include <unordered_map>

#include "particletype.h"
#include "tabulation.h"

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
   * \param s Twice the spin of the multiplet.
   */
  IsoParticleType(const std::string &n, double m, double w, unsigned int s);

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

  /**
   * Returns whether the two IsoParticleType objects have the same PDG code for
   * their first state; if it is, it is the same iso multiplet.
   */
  bool operator==(const IsoParticleType &rhs) const {
    return states_[0]->pdgcode() == rhs.states_[0]->pdgcode();
  }

  /// Returns the name of the multiplet.
  const std::string &name() const { return name_; }

  /// Returns the (average) multiplet mass.
  double mass() const { return mass_; }

  /// Returns the (average) multiplet width.
  double width() const { return width_; }

  /// Returns twice the total isospin of the multiplet
  int isospin() const { return states_.size() - 1; }

  /**
   * Returns twice the spin of the multiplet. All particles in the multiplet
   * are required to have the same spin.
   */
  unsigned int spin() const { return spin_; }

  ParticleTypePtrList get_states() const { return states_; }

  /**
   * Add a new state to an existing multiplet
   * (and check if isospin symmetry is fulfilled).
   */
  void add_state(const ParticleType &type);

  /**
   * Check if there is a multiplet of antiparticles, which is different from
   * the original multiplet.
   */
  bool has_anti_multiplet() const;

  /**
   * Returns a list of all IsoParticleTypes
   */
  static const IsoParticleTypeList &list_all();

  /**
   * Returns the IsoParticleType pointer for the given \p name.
   * If the particle type is not found, an invalid pointer is returned.
   * You can convert the pointer to a bool to check whether it is valid.
   *
   * \note The complexity of the search is \f$\mathcal O(\log N)\f$.
   */
  static const IsoParticleType *try_find(const std::string &name);

  /**
   * Returns the IsoParticleType object for the given \p name.
   *
   * \note The complexity of the search is \f$\mathcal O(\log N)\f$.
   */
  static const IsoParticleType &find(const std::string &name);

  /**
   * Returns the IsoParticleType object for the given \p type.
   *
   * \note The complexity of the search is \f$\mathcal O(\log N)\f$.
   */
  static IsoParticleType *find(const ParticleType &type);

  /// \ingroup exception
  struct ParticleNotFoundFailure : public std::runtime_error {
    using std::runtime_error::runtime_error;
  };

  /**
   * Returns whether the ParticleType with the given \p pdgcode exists.
   *
   * \note The complexity of the search is \f$\mathcal O(\log N)\f$.
   */
  static bool exists(const std::string &name);

  /**
   * Returns the ParticleType object for the given \p name, by first finding the
   * correct multiplet and then looking for the desired state.
   */
  static const ParticleTypePtr find_state(const std::string &name);

  /**
   * Add a new multiplet to the global list of IsoParticleTypes, which contains
   * the type t. If the multiplet exists already, the type t will be added to
   * it.
   */
  static void create_multiplet(const ParticleType &type);

  /// Look up the tabulated resonance integral for the XX -> NR cross section.
  double get_integral_NR(double sqrts);

  /// Look up the tabulated resonance integral for the XX -> RR cross section.
  double get_integral_RR(const ParticleType &type_res_2, double sqrts);

  /// Utility function to help compute various XX->RR spectral integrals
  TabulationPtr integrate_RR(ParticleTypePtr &type_res_2);

  /// Look up the tabulated resonance integral for the XX -> RK cross section.
  double get_integral_RK(double sqrts);

 private:
  /// name of the multiplet
  std::string name_;
  /// (average) mass of the multiplet
  double mass_;
  /// (average) width of the multiplet
  double width_;
  /// twice the spin of the multiplet
  unsigned int spin_;
  /// list of states that are contained in the multiplet
  ParticleTypePtrList states_;

  /// A tabulation of the spectral integral for the NK -> RK cross sections.
  TabulationPtr XS_RK_tabulation_;
  /* A tabulation for the NN -> NR and NN -> DR cross sections,
   * where R is a resonance from this multiplet. */
  TabulationPtr XS_NR_tabulation_, XS_DR_tabulation_;
  /* A tabulation list for the NN -> RR' cross sections,
   * where R is this multiplet and R' is a baryon resonance, associated
   * with a list of resonances R' for the NN -> RR' cross sections;
   * used to calculate every multiplet spectral function only once*/
  std::unordered_map<IsoParticleType *, TabulationPtr> XS_RR_tabulations;

  /**
   * Private version of the 'find' method that returns a non-const reference.
   */
  static IsoParticleType &find_private(const std::string &name);
};

}  // namespace Smash

#endif  // SRC_INCLUDE_ISOPARTICLETYPE_H_
