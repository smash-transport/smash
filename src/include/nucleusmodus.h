/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_NUCLEUSMODUS_H_
#define SRC_INCLUDE_NUCLEUSMODUS_H_

#include "modusdefault.h"

#include "forwarddeclarations.h"
#include "nucleus.h"
#include "deformednucleus.h"
#include "pdgcode.h"

#include <cstring>
#include <memory>
#include <utility>

namespace Smash {

struct ExperimentParameters;

/**
 * \ingroup modus
 * NucleusModus: Provides a modus for colliding nuclei.
 *
 * To use this modus, chose
 * \code
 * General:
 *      MODUS: Nucleus
 * \endcode
 * in the configuration file.
 *
 * Options for NucleusModus go in the "Modi"→"Nucleus" section of the
 * configuration.
 *
 * The following configuration options are understood: \ref input_modi_nucleus_
 */
class NucleusModus : public ModusDefault {
 public:
  /** Constructor
   *
   * Takes all there is to take from the (truncated!) configuration
   * object (only contains configuration for this modus).
   **/
  explicit NucleusModus(Configuration modus_config,
           const ExperimentParameters &parameters);

  /** Creates initial conditions from the particles.
   *
   * In particular, it initializes the nuclei.
   */
  float initial_conditions(Particles *particles,
                          const ExperimentParameters &parameters);

  /// \ingroup exception
  /// Thrown when either \a projectile_ or \a target_ nuclei are empty.
  struct NucleusEmpty : public ModusDefault::BadInput {
    using ModusDefault::BadInput::BadInput;
  };

 private:
  /** Projectile.
   *
   * The object that comes from negative z-values at positive x-values
   * with positive velocity.
   **/
  std::unique_ptr<Nucleus> projectile_;
  /** Target.
   *
   * The object that comes from positive z-values at negative x-values
   * with negative velocity. In fixed target experiments, the target is
   * at rest.
   **/
  std::unique_ptr<Nucleus> target_;
  // Center-of-mass energy of the nucleus-nucleus collision.
  float total_s_;
  /** impact parameter
   *
   * The nuclei projectile_ and target_ will be shifted along the x axis
   * so that their centers move on antiparallel lines that are this
   * distance apart from each other.
   **/
  float impact_ = 0.f;
  /** sample impact parameter
   *
   * sets the impact parameter to a value between min and max.
   *
   * @param s if true, use quadratic sampling (probability for a given
   * impact parameter \f$dP(b)\f$ is proportional to \f$b\f$: \f$dP(b) =
   * b\cdot db\f$), else every \f$b\f$ has same probability.
   * @param min minimum value for impact parameter
   * @param max maximum value for impact parameter
   *
   * Note that max less than min also works fine.
   *
   **/
  void sample_impact(bool s, float min, float max);
  /** initial z displacement of nuclei
   *
   * each nucleus is shifted so that
   * the outermost particle on the side facing the other nucleus is at
   * \f$\pm\f$ this value.
   **/
  double initial_z_displacement_ = 1.0;
  // Reference frame for the system.
  // 1 = Center of velocity
  // 2 = Center of mass
  // 3 = Fixed target
  int frame_ = 1;
  // Get the frame dependent velocity for each nucleus, using
  // the current reference frame. \see frame_
  //
  // @param s The total mandelstam S of the system
  // @param m1 The mass of the projectile.
  // @param m2 The mass of the target.
  // @return < v1, v2 >
  std::pair<double, double> get_velocities(float mandelstam_s, float m1, float m2);

  /**\ingroup logging
   * Writes the initial state for the NucleusModus to the output stream.
   */
  friend std::ostream &operator<<(std::ostream &, const NucleusModus &);
};

}  // namespace Smash

#endif  // SRC_INCLUDE_NUCLEUSMODUS_H_
