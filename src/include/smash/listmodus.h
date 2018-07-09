/*
 *    Copyright (c) 2013-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_LISTMODUS_H_
#define SRC_INCLUDE_LISTMODUS_H_

#include <cmath>
#include <cstdint>
#include <list>
#include <string>
#include <utility>

#include "forwarddeclarations.h"
#include "modusdefault.h"

namespace smash {

/**
 * \ingroup modus
 * ListModus: Provides a modus for running SMASH on an external particle list,
 * for example as an afterburner calculation.
 *
 * To use this modus, choose
    Modus:         List
 * \code
 * General:
 *      Modus: List
 * \endcode
 * in the configuration file.
 *
 * Options for ListModus go in the "Modi"→"List" section of the
 * configuration:
 *
 * \code
 * Modi:
 *      List:
 *              # options here
 * \endcode

 * For configuring see \ref input_modi_list_.
 *
 *
 * Since SMASH is searching for collisions in computational frame time 't',
 * all particles need to be at the same time. If this is not the case in
 * the list provided, the particles will be propagated backwards on
 * straight lines ("anti-freestreaming"). To avoid unphysical interactions
 * of these particles, the back-propagated particles receive a
 * formation_time and zero cross_section_scaling_factor. The cross-sections
 * are set to zero during the time, where the particle will just propagate
 * on a straight line again to appear at the formation_time into the system.
 *
 */
class ListModus : public ModusDefault {
 public:
  /**
   * Constructor
   *
   * Gathers all configuration variables for the List.
   *
   * \param[in] modus_config The configuration object that sets all
   *                         initial conditions of the experiment.
   * \param[in] parameters Unused, but necessary because of templated
   *                       initialization
   */
  explicit ListModus(Configuration modus_config,
                     const ExperimentParameters &parameters);

  /// Construct an empty list. Useful for convenient JetScape connection.
  ListModus() :
    shift_id_(0) {}

  /**
   * Generates initial state of the particles in the system according to a list.
   *
   * \param[out] particles An empty list that gets filled up by this function
   * \param[in] parameters Unused, but necessary because of templated use of
   *                       this function
   * \return The starting time of the simulation
   * \throw runtime_error if an input list file could not be found
   * \throw LoadFailure if an input list file is not correctly formatted
   * \throw invalid_argument if the listed charge of a particle does not
   *                         correspond to its pdg charge
   */
  double initial_conditions(Particles *particles,
                            const ExperimentParameters &parameters);
  /**
   * Judge whether formation times are the same for all the particles;
   * Don't do anti-freestreaming if all particles start already at the same
   * time.
   *
   * If particles are at different times, calculate earliest formation
   * time as start_time_ and free-stream all particles back to this time.
   *
   * \param particles particles to be checked and possibly back-streamed
   */
  void backpropagate_to_same_time(Particles &particles);

  /**
   * Tries to add a new particle to particles and performs consistency checks:
   * (i) The PDG code is legal and exists in SMASH. If not, a warning is printed
   *     and the particle is ignored.
   * (ii) The mass matches the pole mass of `pdgcode` in SMASH. If it does not,
   *      then a warning is printed, the pole mass of the particle is set equal
   *      to the corresponding mass from SMASH particle table and it's energy
   *      is recomputed as \f$ E^2 = p^2 + m^2 \f$.
   * (iii) Any stable particle is on-shell, i.e. \f$ E^2 - p^2 = m^2 \f$. If it
   *       is not, then a warning is printed and the energy is set to
   *       \f$ E^2 = p^2 + m^2 \f$.
   * This very tolerant behaviour is justified by the practical usage of SMASH
   * as afterburner. Usually particles unknown to SMASH are rare resonances,
   * which do not play a large role. Mass mismatch is typically less than 1%
   * and comes from rounding and from SMASH enforcing isospin symmetry
   * (for example the mass of neutral pion is artificially forced to be the
   * same as charged pion). On-shellness violation typically comes from the
   * insufficient number of significant digits in the input file + rounding.
   *
   * \param[in] pdgcode pdg code of added particle
   * \param[in] pdgcode pdg code of added particle
   * \param[in] t       time of added particle
   * \param[in] x       x-coordinate of added particle
   * \param[in] y       y-coordinate of added particle
   * \param[in] z       z-coordinate of added particle
   * \param[in] mass    mass of added particle
   * \param[in] E       energy of added particle
   * \param[in] px      x-component of momentum of added particle
   * \param[in] py      y-component of momentum of added particle
   * \param[in] pz      z-component of momentum of added particle
   * \param[out] particles structure, to which the particle is added
   */
  void try_create_particle(Particles &particles, PdgCode pdgcode, double t,
                           double x, double y, double z, double mass, double E,
                           double px, double py, double pz);

  /** \ingroup exception
   * Used when external particle list cannot be found.
   */
  struct LoadFailure : public std::runtime_error {
    using std::runtime_error::runtime_error;
  };

 protected:
  /// Starting time for the List; changed to the earliest formation time
  double start_time_ = 0.;

 private:
  /// File directory of the particle list
  std::string particle_list_file_directory_;

  /// File prefix of the particle list
  std::string particle_list_file_prefix_;

  /// File name of current file
  std::string current_particle_list_file_;

  /// shift_id is the start number of file_id_
  const int shift_id_;

  /// event_id_ = the unique id of the current event
  int event_id_;

  /// file_id_ is the id of the current file
  int file_id_;

  /// Counter for mass-check warnings to avoid spamming
  int n_warns_precision_ = 0;
  /// Counter for energy-momentum conservation warnings to avoid spamming
  int n_warns_mass_consistency_ = 0;

  /** Check if the file given by filepath has events left after streampos
   * last_position
   *
   * \param[in] filepath Path to file to be checked.
   * \param[in] last_position Streamposition in file after which check is
   * performed
   * \return True if there is at least one event left, false otherwise
   * \throws runtime_error If file could not be read for whatever reason.
   */
  bool file_has_events_(bf::path filepath, std::streampos last_position);

  /// last read position in current file
  std::streampos last_read_position_;

  /** Return the absolute file path based on given integer. The filename
   * is assumed to have the form (particle_list_prefix)_(file_id)
   *
   * \param[in] file_id integer of wanted file
   * \return Absolute file path to file
   * \throws
   * runtime_error if file does not exist.
   */
  bf::path file_path_(const int file_id);

  /**  Read the next event. Either from the current file if it has more events
   * or from the next file (with file_id += 1)
   *
   * \returns
   *  One event as string.
   *  \throws runtime_error If file could not be read for whatever reason.
   */
  std::string next_event_();

  /**\ingroup logging
   * Writes the initial state for the List to the output stream.
   *
   * \param[in] out The ostream into which to output
   * \param[in] m The ListModus object to write into out
   */
  friend std::ostream &operator<<(std::ostream &, const ListModus &);
};

}  // namespace smash

#endif  // SRC_INCLUDE_LISTMODUS_H_
