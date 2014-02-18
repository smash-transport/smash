/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_NUCLEUS_H_
#define SRC_INCLUDE_NUCLEUS_H_

#include<vector>

#include "include/particles.h"

/// A Nucleus is a collection of Particles (ParticleData thingys) that
/// are initialized before the beginning of the simulation and all have
/// the same velocity (and spatial proximity).
/// This class inherits from Particles, which is the collection of all
/// particles in the simulation and contains special functions for the
/// initialization of nuclei.
class Nucleus : public Particles {
 public:
  /// returns the mass of the nucleus
  float mass() const;
  /// returns the radius of the nucleus
  float nuclear_radius() const;
  /** returns a length with a probability distribution according to a
   * Woods-Saxon distribution suitable for this nucleus.
   * $$\frac{dN}{dr} = \frac{r^2}{\exp\left(\frac{r-R}{d}\right) + 1}$$
   * where $d$ is the softness_ parameter and $R$ is nuclear_radius(). */
  float distribution_nucleons() const;
  /// sets the positions of the nuclei inside nucleus A.
  void arrange_nucleons();
  /// Boosts the nuclei (no shifting yet!)
  void boost(const double& beta_squared);
  /// Adds a particle to the nucleus
  void add_particle(const int pdgcode);

 private:
  /** softness of Woods-Saxon-distribution in this nucleus im fm
   * (for softness_ == 0, we obtain a hard sphere (but don't do that; we
   * don't want to divide by zero). */
  float softness_ = .545f;
  /** single-proton-radius */
  float proton_radius_ = 1.2f;
  /// z (beam direction-) coordinate of the outermost particle
  float z_max_ = 0.f;
  /// z (beam direction-) coordinate of the outermost particle
  float z_min_ = 0.f;
};

#endif  // SRC_INCLUDE_NUCLEUS_H_
