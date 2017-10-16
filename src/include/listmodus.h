/*
 *    Copyright (c) 2013-2017
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

namespace Smash {

/**
 * \ingroup modus
 * ListModus: Provides a modus for running SMASH on any external particle list,
 * for example as an afterburner calculations for hybrid codes
 *
 * For configuring see \ref input_modi_list_.
*
*#!OSCAR2013 particle_lists SMASH-githash t x y z mass p0 px py pz pdg
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
  /* default constructor with probable values */
  explicit ListModus(Configuration modus_config,
                     const ExperimentParameters &parameters);

  /** creates initial conditions for the particles.
   */
  double initial_conditions(Particles *particles,
                            const ExperimentParameters &parameters);

  /// \ingroup exception
  struct LoadFailure : public std::runtime_error {
    using std::runtime_error::runtime_error;
  };

 private:
  /// File directory of the particle list
  std::string particle_list_file_directory_;

  /// File prefix of the particle list
  std::string particle_list_file_prefix_;

  /// Starting time for the List; changed to the earliest formation time
  double start_time_ = 0.;

  /// shift_id is the start number of event_id
  const int shift_id_;

  /// event_id_ = the unique id of the current even
  int event_id_;

  /// check whether anti-freestreaming is needed, if yes return
  /// earliest formation time as start_time_
  std::pair<bool, double> check_formation_time_(
      const std::string &particle_list);

  /**\ingroup logging
   * Writes the initial state for the List to the output stream.
   */
  friend std::ostream &operator<<(std::ostream &, const ListModus &);
};

}  // namespace Smash

#endif  // SRC_INCLUDE_LISTMODUS_H_
