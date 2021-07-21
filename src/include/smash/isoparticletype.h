/*
 *    Copyright (c) 2015-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_SMASH_ISOPARTICLETYPE_H_
#define SRC_INCLUDE_SMASH_ISOPARTICLETYPE_H_

#include <string>
#include <unordered_map>
#include <vector>

#include "particletype.h"
#include "sha256.h"
#include "tabulation.h"

namespace smash {

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
   * \param p Parity of the multiplet.
   */
  IsoParticleType(const std::string &n, double m, double w, unsigned int s,
                  Parity p);

  /**
   * Copies are not allowed as they break intended use. Instead use a const-ref
   * or ParticleTypePtr (as returned from operator&).
   */
  IsoParticleType(const IsoParticleType &) = delete;
  /// Assignment is not allowed, see copy constructor above
  IsoParticleType &operator=(const IsoParticleType &) = delete;

  /// Move constructor of IsoParticleType (needed for std::sort)
  IsoParticleType(IsoParticleType &&) = default;
  /// Move constructor of IsoParticleType "="-operator (needed for std::sort)
  IsoParticleType &operator=(IsoParticleType &&) = default;

  /**
   * Returns whether the two IsoParticleType objects have the same PDG code for
   * their first state; if so, it is the same iso multiplet.
   *
   * \param rhs The other multiplet.
   */
  bool operator==(const IsoParticleType &rhs) const {
    return states_[0]->pdgcode() == rhs.states_[0]->pdgcode();
  }

  /// Returns the name of the multiplet.
  const std::string &name() const { return name_; }

  /// Returns the name of the multiplet, after replacing "'" with "_prime"
  const std::string name_filtered_prime() const {
    std::string tmp_s = name_;
    std::size_t found_position = tmp_s.find("'");
    if (found_position != std::string::npos) {
      tmp_s.erase(found_position, 1);
      tmp_s.insert(found_position, "_prime");
    }
    return tmp_s;
  }

  /// Returns the (average) multiplet mass.
  double mass() const { return mass_; }

  /// Returns the (average) multiplet width.
  double width() const { return width_; }

  /// Returns twice the total isospin of the multiplet.
  int isospin() const { return states_.size() - 1; }

  /**
   * Returns twice the spin of the multiplet. All particles in the multiplet
   * are required to have the same spin.
   */
  unsigned int spin() const { return spin_; }

  /**
   * \return The parity of the multiplet.
   */
  Parity parity() const { return parity_; }

  /**
   * \return Is this a hadron multiplet?
   */
  bool is_hadron() const { return states_[0]->is_hadron(); }

  /// Returns list of states that form part of the multiplet.
  ParticleTypePtrList get_states() const { return states_; }

  /**
   * Add a new state to an existing multiplet
   * (and check if isospin symmetry is fulfilled).
   *
   * \param type The particle state to be added.
   */
  void add_state(const ParticleType &type);

  /**
   * Return a multiplet of antiparticles, if it is different from
   * the original multiplet. Otherwise, return a nullptr.
   */
  const IsoParticleType *anti_multiplet() const;

  /**
   * Check if there is a multiplet of antiparticles, which is different from
   * the original multiplet.
   */
  bool has_anti_multiplet() const;

  /// Returns a list of all IsoParticleTypes.
  static const IsoParticleTypeList &list_all();

  /// Returns a list of all IsoParticleTypes that are baryon
  /// resonances.
  static const std::vector<const IsoParticleType *> list_baryon_resonances();

  /**
   * Returns the IsoParticleType pointer for the given \p name.
   * If the particle type is not found, an invalid pointer is returned.
   * You can convert the pointer to a bool to check whether it is valid.
   *
   * \param name The name of the particle type.
   * \note The complexity of the search is \f$\mathcal O(\log N)\f$.
   */
  static const IsoParticleType *try_find(const std::string &name);

  /**
   * Returns the IsoParticleType object for the given \p name.
   *
   * \param name The name of the of the particle type to be found.
   * \throw ParticleNotFoundFailure if \p name not found.
   * \note The complexity of the search is \f$\mathcal O(\log N)\f$.
   */
  static const IsoParticleType &find(const std::string &name);

  /**
   * Returns the IsoParticleType object for the given \p type.
   *
   * \param type The particle type to be found.
   * \throw ParticleNotFoundFailure if \p type not found.
   * \note The complexity of the search is \f$\mathcal O(\log N)\f$.
   */
  static IsoParticleType *find(const ParticleType &type);

  /** \ingroup exception
   *
   * Throw when requested particle could not be found.
   */
  struct ParticleNotFoundFailure : public std::runtime_error {
    using std::runtime_error::runtime_error;
  };

  /**
   * Returns whether the ParticleType with the given \p pdgcode exists.
   *
   * \param name The name of the particle type to be found.
   * \note The complexity of the search is \f$\mathcal O(\log N)\f$.
   */
  static bool exists(const std::string &name);

  /**
   * Returns the ParticleType object for the given \p name, by first finding the
   * correct multiplet and then looking for the desired state.
   *
   * \param name The name of the particle state to be found.
   * \throw std::runtime_error if \p name is not found.
   */
  static const ParticleTypePtr find_state(const std::string &name);

  /**
   * Add a new multiplet to the global list of IsoParticleTypes, which contains
   * \p type. If the multiplet exists already, the \p type will be added to
   * it.
   *
   * \param type The multiplet to be created.
   */
  static void create_multiplet(const ParticleType &type);

  /**
   * Tabulate all relevant integrals.
   *
   * \param hash The hash of the particle properties.
   *             This is used to determine whether a cached tabulation can be
   *             reused or not.
   * \param tabulations_path The path to the directory where the tabulations are
   * cached.
   */
  static void tabulate_integrals(sha256::Hash hash,
                                 const bf::path &tabulations_path);

  /**
   * Look up the tabulated resonance integral for the XX -> NR cross section.
   *
   * \param sqrts The center-of-mass energy.
   */
  double get_integral_NR(double sqrts);

  /**
   * Look up the tabulated resonance integral for the XX -> RR cross section.
   *
   * \param type_res_2 Type of the two resonances in the final state.
   * \param sqrts The center-of-mass energy.
   */
  double get_integral_RR(IsoParticleType *type_res_2, double sqrts);

  /**
   * Look up the tabulated resonance integral for the XX -> RK cross section.
   *
   * \param sqrts The center-of-mass energy.
   */
  double get_integral_RK(double sqrts);

  /**
   * Look up the tabulated resonance integral for the XX -> piR cross section.
   *
   * \param sqrts The center-of-mass energy.
   */
  double get_integral_piR(double sqrts);

  /**
   * Look up the tabulated resonance integral for the XX -> rhoR cross section.
   *
   * \param sqrts The center-of-mass energy.
   */
  double get_integral_rhoR(double sqrts);

 private:
  /// name of the multiplet
  std::string name_;
  /// (average) mass of the multiplet
  double mass_;
  /// (average) width of the multiplet
  double width_;
  /// twice the spin of the multiplet
  unsigned int spin_;
  /// parity of the multiplet
  Parity parity_;
  /// list of states that are contained in the multiplet
  ParticleTypePtrList states_;

  /// A tabulation of the spectral integral for the dpi -> d'pi cross sections.
  Tabulation *XS_piR_tabulation_ = nullptr;
  /// A tabulation of the spectral integral for the NK -> RK cross sections.
  Tabulation *XS_RK_tabulation_ = nullptr;
  /**
   * A tabulation for the NN -> NR cross sections,
   * where R is a resonance from this multiplet.
   */
  Tabulation *XS_NR_tabulation_ = nullptr;
  /**
   * A tabulation for the NN -> RΔ cross sections,
   * where R is a resonance from this multiplet.
   */
  Tabulation *XS_DeltaR_tabulation_ = nullptr;
  /**
   * A tabulation for the ρρ integrals.
   */
  Tabulation *XS_rhoR_tabulation_ = nullptr;

  /**
   * Private version of the 'find' method that returns a non-const reference.
   *
   * \param[in] name The name of the of the particle type to be found.
   * \throw ParticleNotFoundFailure if \p name not found.
   */
  static IsoParticleType &find_private(const std::string &name);
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_ISOPARTICLETYPE_H_
