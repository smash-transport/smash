/*
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_SPHEREMODUS_H_
#define SRC_INCLUDE_SPHEREMODUS_H_

#include <stdint.h>
#include <cmath>
#include <list>

#include "forwarddeclarations.h"
#include "modusdefault.h"

namespace Smash {

/**
 * \ingroup modus
 * SphereModus: Provides a modus for expanding matter calculations
*
*  Matter is put in a sphere of radius R with isotropic thermal momenta.
*
* To use this modus, chose
* \code
* General:
*      MODUS: Sphere
* \endcode
* in the configuration file.
*
* Options for SphereModus go in the "Modi"→"Sphere" section of the
* configuration:
*
* \code
* Modi:
*      Sphere:
*              # definitions here
* \endcode
*
* The following configuration options are understood: \ref input_modi_sphere_
*/
class SphereModus : public ModusDefault {
 public:
  /* default constructor with probable values */
  explicit SphereModus(Configuration modus_config,
                       const ExperimentParameters &parameters);

  /** creates initial conditions for the particles.
   */
  float initial_conditions(Particles *particles,
                          const ExperimentParameters &parameters);
 private:
  /// Sphere radius length
  float radius_;
  /// Temperature for momentum distribution
  float sphere_temperature_;
  /// Starting time for the Sphere
  const float start_time_ = 0.0f;
  /// particle multiplicities at initialization
  const std::map<PdgCode, int> init_multipl_;
  /**\ingroup logging
   * Writes the initial state for the Sphere to the output stream.
   */
  friend std::ostream &operator<<(std::ostream &, const SphereModus &);
};
}  // namespace Smash
#endif  // SRC_INCLUDE_SPHEREMODUS_H_
