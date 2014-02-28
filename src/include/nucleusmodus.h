/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_NUCLEUSMODUS_H_
#define SRC_INCLUDE_NUCLEUSMODUS_H_

#include <stdint.h>
#include <cmath>
#include <list>
#include <string>

#include "include/configuration.h"
#include "include/modusdefault.h"
#include "include/nucleus.h"
#include "include/particles.h"
#include "include/parameters.h"

struct ExperimentParameters;

/** NucleusModus: Provides a modus for colliding nuclei.
 *
 * To use this modus, chose
 * \code
 * General:
 *      MODUS: Nucleus
 * \endcode
 * in the configuration file.
 *
 * Options for NucleusModus go in the "Modi"→"Nucleus" section of the
 * configuration:
 *
 * \code
 * Modi:
 *      Nucleus:
 *              # definitions here
 * \endcode
 *
 * The following directives are understood:
 *
 * `SQRTSNN:` Defines the energy of the collision as center-of-mass
 * energy in the collision of one participant each from both nuclei.
 * Since not all participants have the same mass, and hence
 * \f$\sqrt{s_{\rm NN}}\f$ is different for \f$NN\f$ = proton+proton and
 * \f$NN\f$=neutron+neutron, you can specify which \f$NN\f$-pair you
 * want this to refer to with `SQRTS_N`
 *
 * `SQRTS_N:` The "N" that's used for `SQRTSNN`. Expects a vector of two
 * PDG Codes, e.g. `SQRTS_N: [2212, 2212]` for proton-proton (that's the
 * default behaviour). The important part here is which mass is used for
 * the calculation \f$\sqrt{s_{\rm NN}} \rightarrow \beta\f$.
 *
 * \todo It should be possible to say something like "take the average
 * mass of a particle in the nuclei". Maybe even "take the average mass
 * of a particle in this nucleus" (and specify the other one). This can
 * be done, e.g., by using a unknown or otherwise invalid PDG code.
 *
 * `Projectile:` Section for projectile nucleus. The projectile will
 * start at \f$z < 0\f$ and fly in positive \f$z\f$-direction, at \f$x
 * \ge 0\f$.
 *
 * `Target:` Section for target nucleus. The target will start at \f$z
 * > 0\f$ and fly in negative \f$z\f$-direction, at \f$x \le 0\f$.
 *
 * `Projectile:`/`Target:` options:
 * \li `PARTICLES:` A map in which the keys are PDG codes and the
 * values are number of particles with that PDG code that should be in
 * the current nucleus. E.g.\ `PARTICLES: {2212: 82, 2112: 126}` for a
 * lead-208 nucleus (82 protons and 126 neutrons = 208 nucleons), and
 * `PARTICLES: {2212: 1, 2112: 1, 3122: 1}` for Hyper-Triton (one
 * proton, one neutron and one Lambda).
 * \li `SOFTNESS:` The softness used in the Woods-Saxon distribution
 * for this nucleus.
 * 0 means a hard sphere.
 *
 * `Impact:` A section for the impact parameter (= distance of the two
 * straight lines that the center of masses of the nuclei travel on).
 *
 * \li `VALUE: `fixed value for the impact parameter. No other \a
 * Impact: directive is looked at.
 * \li `SAMPLE:` if `linear`, use linear sampling of the impact
 * parameter (\f$dP(b) = db\f$). If else, use quadratical input sampling
 * (the physical case, where the probability of an input parameter range
 * is proportional to the area corresponding to that range, \f$dP(b) =
 * b\cdot db\f$).
 * \li `RANGE:` A vector of minimal and maximal impact parameters
 * between which b should be chosen. (The order of these is not
 * important.)
 * \li `MAX:` Like `RANGE: [0.0, MAX]`. Note that if both `RANGE` and
 * `MAX` are specified, `MAX` takes precedence.
 *
 * `INITIAL_DISTANCE:` The initial distance of the two nuclei. That
 * means \f$z_{\rm min}^{\rm target} - z_{\rm max}^{\rm projectile}\f$.
 **/
class NucleusModus : public ModusDefault {
 public:
  /** Constructor
   *
   * Takes all there is to take from the (truncated!) configuration
   * object (only contains configuration for this modus).
   **/
  explicit NucleusModus(Configuration modus_config);

  void print_startup();

  void initial_conditions(Particles *particles,
                          const ExperimentParameters &parameters);

 private:
  /** Projectile.
   *
   * The object that comes from negative z-values at positive x-values
   * with positive velocity.
   **/
  Nucleus projectile_;
  /** Target.
   *
   * The object that comes from positive z-values at negative x-values
   * with negative velocity. In fixed target experiments, the target is
   * at rest.
   **/
  Nucleus target_;
  /** Center-of-mass energy of the individual nucleon-nucleon
   * collisions.
   *
   * Note that \f$\sqrt{s}\f$ is different for neutron-neutron and
   * proton-proton collisions (because of the different masses).
   * Therefore, pdg_sNN_1_ and pdg_sNN_2_ are needed to specify which
   * two particles' collisions have this \f$\sqrt{s}\f$. (They each
   * specify the PDG code of the particle species that we want to use
   * for the definition of \f$\sqrt{s}\f$).
   **/
  float sqrt_s_NN_;
  /// \see sqrt_s_NN_
  int pdg_sNN_1_ = 2212;
  /// \see sqrt_s_NN_
  int pdg_sNN_2_ = 2212;
  /** impact paramter
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
  void sample_impact(const bool s, const float min, const float max);
  /** initial z displacement of nuclei
   *
   * each nucleus is shifted so that
   * the outermost particle on the side facing the other nucleus is at
   * \f$\pm\f$ this value.
   **/
  double initial_z_displacement_ = 1.0;
};

#endif  // SRC_INCLUDE_NUCLEUSMODUS_H_
