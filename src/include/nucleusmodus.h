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
  double initial_z_displacement = 1.0;
};

#endif  // SRC_INCLUDE_NUCLEUSMODUS_H_
