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

class Nucleus : public Particles {
 public:
  /** returns the mass of the nucleus */
  float get_mass() const;
  /** returns the radius of the nucleus */
  float nuclear_radius() const;
  /* returns a length r with a relative probability of $$\frac{dN}{dr} =
   * \frac{r^2}{\exp\left(\frac{r-R}{d}\right) + 1}$$ where $d$ is a class
   * parameter and $R$ is the function parameter ''radius''. */
  float distribution_nucleons() const;
  /// sets the positions of the nuclei inside nucleus A.
  void arrange_nucleons();
  /// Boosts the nuclei (no shifting yet!)
  void boost(const double& beta_squared);
 private:
  /** softness of Woods-Saxon-distribution in this nucleus im fm (this
   * is the parameter $d$ in the distribution; for softness_ == 0, we
   * obtain a hard sphere (but don't do that; we don't want to divide by
   * zero). */
  float softness_ = .545f;
  /** single-proton-radius */
  float proton_radius_ = 1.2f;
  /** z (beam direction-) coordinate of the outermost particle */
  float max_z_ = 0.0;
};

#endif  // SRC_INCLUDE_NUCLEUS_H_
