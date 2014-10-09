/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_COLLIDERMODUS_H_
#define SRC_INCLUDE_COLLIDERMODUS_H_

#include "modusdefault.h"

#include "forwarddeclarations.h"
#include "pdgcode.h"

namespace Smash {

/**
 * \ingroup modus
 * ColliderModus: Provides a modus for collisions of single particles
 *
 * To use this modus, chose
 * \code
 * General:
 *      MODUS: Collider
 * \endcode
 * in the configuration file.
 *
 * Options for ColliderModus go in the "Modi"→"Collider" section of the
 * configuration:
 *
 * \code
 * Modi:
 *      Collider:
 *              # definitions here
 * \endcode
 *
 * The following configuration options are understood: \ref input_modi_collider_
 */
class ColliderModus : public ModusDefault {
 public:
  /** Constructor
   *
   * Takes all there is to take from the (truncated!) configuration
   * object (only contains configuration for this modus).
   **/
  explicit ColliderModus(Configuration modus_config,
           const ExperimentParameters &parameters);

  /** creates initial conditions from the particles.
   *
   * In particular, it initializes target and projectile.
   */
  float initial_conditions(Particles *particles,
                          const ExperimentParameters &parameters);

  // in ModusDefault:
  // * sanity_check
  // * check_collision_geometry
  // * propagate

 private:
  /// PdgCode of Projectile particle
  const PdgCode projectile_;
  /// PdgCode of Target particle
  const PdgCode target_;
  /// Center-of-mass energy of the collision in GeV
  const float sqrts_;

  /**\ingroup logging
   * Writes the initial state for the Collider to the output stream.
   */
  friend std::ostream &operator<<(std::ostream &, const ColliderModus &);
};

}  // namespace Smash

#endif  // SRC_INCLUDE_COLLIDERMODUS_H_
