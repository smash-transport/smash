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

/** SphereModus: Provides a modus for expanding matter calculations
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
* The following directives are understood:
*
* Modi:Sphere:
*--------
*/
// Userguide {
/**
* \if user
* \page input_modi_box_ Input Section Modi:Box
* \endif
*
* `RADIUS`: Radius of the Sphere.
*
* `NUMBEROFPARTICLES`: Total number of particles in the Sphere.
*
* `SPHERETEMPERATURE`: Temperature for the momentum sampling in the sphere in GeV.
*
* `START_TIME`: Starting time of Sphere calculation.
*/
// } Userguide
class SphereModus : public ModusDefault {
 public:
  /* default constructor with probable values */
  explicit SphereModus(Configuration modus_config,
                       const ExperimentParameters &parameters);
  /** Prints some information about the initialization of SphereModus
   *
   * \see ModusDefalt::print_startup()
   */
  void print_startup();  // TODO(mkretz): needs to be discoverable from an
    // outside "printer"
    /** creates initial conditions for the particles.
     */
  float initial_conditions(Particles *particles,
                          const ExperimentParameters &parameters);
 private:
  /// Sphere radius length
  float radius_;
  /// Total number of particles in Sphere
  int number_of_particles_;
  /// Temperature for momentum distribution
  float sphere_temperature_;
  /// Starting time for the Sphere
  const float start_time_ = 0.0f;
};
}  // namespace Smash
#endif  // SRC_INCLUDE_SPHEREMODUS_H_
