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
#include <map>

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
  double initial_conditions(Particles *particles,
                          const ExperimentParameters &parameters);

 private:
  /// Sphere radius length
  double radius_;
  /// Temperature for momentum distribution
  double sphere_temperature_;
  /// Starting time for the Sphere
  const double start_time_ = 0.;
  /** whether to use a thermal initialization for all particles
   *  instead of specific numbers */
  const bool use_thermal_ = false;
  /// baryon chemical potential for thermal box
  const double mub_;
  /// strange chemical potential for thermal box
  const double mus_;
  /// particle multiplicities at initialization
  const std::map<PdgCode, int> init_multipl_;
  /** Initialization scheme for momenta in the sphere;
   *  used for expanding metric setup */
  const SphereInitialCondition init_distr_;
  /**\ingroup logging
   * Writes the initial state for the Sphere to the output stream.
   */
  friend std::ostream &operator<<(std::ostream &, const SphereModus &);
};
}  // namespace Smash
#endif  // SRC_INCLUDE_SPHEREMODUS_H_
