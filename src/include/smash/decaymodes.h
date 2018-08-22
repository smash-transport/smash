/*
 *    Copyright (c) 2013-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_DECAYMODES_H_
#define SRC_INCLUDE_DECAYMODES_H_

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "processbranch.h"

namespace smash {

/**
 * \ingroup data
 *
 * The DecayModes class is used to store and update information about decay
 * branches (i.e. the possible children a mother particle can decay into),
 * including their relative weights.
 *
 * If you want to find a DecayModes object for a specific particle type use
 * ParticleType::decay_modes().
 */
class DecayModes {
 public:
  /**
   * Add a decay mode using all necessary information
   *
   * \param[in] mother the particle which decays
   * \param[in] ratio the weight to add to the current mode
   * \param[in] L angular momentum
   * \param[in] particle_types a list of the products of the decay
   */
  void add_mode(ParticleTypePtr mother, double ratio, int L,
                ParticleTypePtrList particle_types);

  /**
   * Add a decay mode from an already existing decay branch
   *
   * \param[in] branch the decay branch to add
   */
  void add_mode(DecayBranchPtr branch) {
    decay_modes_.push_back(std::move(branch));
  }

  /**
   * Renormalize the branching ratios to add up to 1.
   *
   * \param[in] name the name of the decaying particle
   */
  void renormalize(const std::string &name);

  /// \return true if empty (i.e. no decay modes)
  bool is_empty() const { return decay_modes_.empty(); }

  /// \return pass out the decay modes list
  const DecayBranchList &decay_mode_list() const { return decay_modes_; }

  /**
   * Loads the DecayModes map as described in the \p input string.
   *
   * It does sanity checking - that the particles it talks about are in the
   * ParticleType map.
   *
   * \param[in] input the full decaymodes input file, as a string
   * \throw MissingDecays if there are no decays specified for an unstable
   *                      particle
   * \throw LoadFailure if there are duplicate entries detailing decaymodes
   *                    for the same particle, or if the angular momentum is
   *                    smaller than 0 or larger than 4
   * \throw InvalidDecay if product of the decay does not exist in the
   *                     particle list or as an isospin multiplet, or if
   *                     decay is forbidden by isospin or charge conservation,
   *                     or if sum of the minimum mass of products is less
   *                     than the pole mass of the mother particle
   * \throw runtime_error if there are less than 2 or more than 3 products
   *                      to a given decay branch, or if a branch cannot be
   *                      found to be part of an isospin multiplet
   */
  static void load_decaymodes(const std::string &input);

  /**
   * Retrieve a decay type.
   *
   * \param[in] mother the decaying particle
   * \param[in] particle_types the products of the decay
   * \param[in] L the angular momentum
   * \return the corresponding DecayType object
   * \throw InvalidDecay if there are less than 2 or more than 3 products
   */
  static DecayType *get_decay_type(ParticleTypePtr mother,
                                   ParticleTypePtrList particle_types, int L);

  /// \ingroup exception
  struct InvalidDecay : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };
  /// \ingroup exception
  struct LoadFailure : public std::runtime_error {
    using std::runtime_error::runtime_error;
  };
  /// \ingroup exception
  struct MissingDecays : public LoadFailure {
    using LoadFailure::LoadFailure;
  };
  /// \ingroup exception
  struct ReferencedParticleNotFound : public LoadFailure {
    using LoadFailure::LoadFailure;
  };

 private:
  /**
   * Vector of decay modes.
   * Each mode consists of a vector of the pdg codes of decay products
   * and a ratio of this decay mode compared to all possible modes
   */
  DecayBranchList decay_modes_;

  /// allow ParticleType::decay_modes to access all_decay_modes
  friend const DecayModes &ParticleType::decay_modes() const;

  /**
   * A list of all DecayModes objects using the same indexing as
   * all_particle_types.
   */
  static std::vector<DecayModes> *all_decay_modes;
};

}  // namespace smash

#endif  // SRC_INCLUDE_DECAYMODES_H_
