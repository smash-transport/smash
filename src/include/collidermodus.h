/*
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_COLLIDERMODUS_H_
#define SRC_INCLUDE_COLLIDERMODUS_H_

#include <cstring>
#include <memory>
#include <utility>

#include "deformednucleus.h"
#include "forwarddeclarations.h"
#include "interpolation.h"
#include "modusdefault.h"
#include "nucleus.h"
#include "pdgcode.h"

namespace Smash {

struct ExperimentParameters;

enum class Sampling {
    UNIFORM,
    QUADRATIC,
    CUSTOM,
};

/**
 * \ingroup modus
 * ColliderModus: Provides a modus for colliding nuclei.
 *
 * To use this modus, chose
 * \code
 * General:
 *      MODUS: Collider
 * \endcode
 * in the configuration file.
 *
 * Options for ColliderModus go in the "Modi"→"Collider" section of the
 * configuration.
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

  /** Creates initial conditions from the particles.
   *
   * In particular, it initializes the nuclei.
   */
  float initial_conditions(Particles *particles,
                          const ExperimentParameters &parameters);

  /// \ingroup exception
  /// Thrown when either \a projectile_ or \a target_ nuclei are empty.
  struct ColliderEmpty : public ModusDefault::BadInput {
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
  /** Center-of-mass energy of the nucleus-nucleus collision.
   *
   **/
  float total_s_;
  /** Impact parameter.
   *
   * The nuclei projectile_ and target_ will be shifted along the x axis
   * so that their centers move on antiparallel lines that are this
   * distance apart from each other.
   **/
  float impact_ = 0.f;
  /// Method used for sampling of impact parameter.
  Sampling sampling_ = Sampling::QUADRATIC;
  /// Minimum value of impact parameter.
  float imp_min_ = 0.0;
  /// Maximum value of impact parameter.
  float imp_max_ = 0.0;
  /// Maximum value of yield. Needed for custom impact parameter sampling.
  float yield_max_ = 0.0;
  /// Pointer to the impact parameter interpolation.
  std::unique_ptr<InterpolateData<float>> impact_interpolation_ = nullptr;
  /** Sample impact parameter.
   *
   * Samples the impact parameter from values between imp_min_ and imp_max_, if
   * linear or quadratic sampling is used. By specifying impact parameters and
   * corresponding yields, custom sampling can be used.
   * This depends on the value of sampling_.
   *
   * Note that imp_max_ less than imp_min_ also works fine.
   *
   **/
  void sample_impact();
  /** Initial z displacement of nuclei.
   *
   * Each nucleus is shifted so that
   * the outermost particle on the side facing the other nucleus is at
   * \f$\pm\f$ this value.
   *
   * \fpPrecision Why \c double?
   **/
  double initial_z_displacement_ = 1.0;
  /** Reference frame for the system.
   *
   * 1 = Center of velocity<br>
   * 2 = Center of mass<br>
   * 3 = Fixed target<br>
   **/
  int frame_ = 1;
  /**
   * An option to include Fermi motion
   */
  bool fermi_motion_;
  /** Get the frame dependent velocity for each nucleus, using
   * the current reference frame. \see frame_
   *
   * @param mandelstam_s The total center-of-mass energy of the system.
   * @param m_a The mass of the projectile.
   * @param m_b The mass of the target.
   * @return < v_a, v_b > Velocities of the nuclei.
   *
   * \fpPrecision Why \c double?
   **/
  std::pair<double, double> get_velocities(float mandelstam_s,
                                           float m_a, float m_b);

  /**\ingroup logging
   * Writes the initial state for the ColliderModus to the output stream.
   */
  friend std::ostream &operator<<(std::ostream &, const ColliderModus &);
};

}  // namespace Smash

#endif  // SRC_INCLUDE_COLLIDERMODUS_H_
