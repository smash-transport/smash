/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_DEFORMEDNUCLEUS_H_
#define SRC_INCLUDE_DEFORMEDNUCLEUS_H_

#include "configuration.h"
#include "forwarddeclarations.h"
#include "nucleus.h"
#include "threevector.h"

namespace Smash {

class DeformedNucleus : public Nucleus {
 public:
  DeformedNucleus();

  /** Return the deformed Woods-Saxon probability for
   * the given position.
   * 
   * @param r The sample radius
   * @param cosx The cosine of the sample polar angle
   * @return The woods-saxon probability
   **/
  double deformed_woods_saxon(double r, double cosx) const;

  /**  Deformed Woods-Saxon sampling routine.
   *
   * @return a spatial position from uniformly sampling 
   * the deformed woods-saxon distribution
   **/
  virtual ThreeVector distribute_nucleon() const;

  /** Sets the deformation parameters of the Woods-Saxon distribution
   * according to the current mass number.
   *
   * Ref. for deformation parameters is "Nuclear Ground-State Masses and Deformations" 
   * by P. Möller, J. R. Nix, W. D. Myers, and W. J. Swiatecki.
   * Corrections to deformation parameter beta2 in Uranium come from 
   * arxiv:nuclth/0506088 by A. Kuhlman and U. Heinz.
   * For finite nucleon size corrections to the nuclear density and radius, see ref.
   * arxiv:0904.4080 [nucl-th] by T. Hirano and Y. Nara for copper and gold, and 
   * arxiv:1010.6222 [nucl-th] by T. Hirano, P. Huovinen, and Y. Nara for uranium.
   */
  virtual void set_parameters_automatic();

  // Set parameters for Woods-Saxon by hand using the configuration file.
  // \see Nucleus::set_parameters_from_config
  virtual void set_parameters_from_config(bool is_projectile, Configuration &config);

  // Rotates the nucleus according to members nucleus_polar_angle_
  // and nucleus_azimuthal_angle_ and updates nucleon positions.
  virtual void rotate();

  // Spherical harmonics Y_2_0 and Y_4_0.
  double y_l_0(int l, double cosx) const;
  
  // Set deformation coefficient for Y_2_0.
  inline void set_beta_2(double b2) {
    beta2_ = b2;
  }
  // Set deformation coefficient for Y_4_0.
  inline void set_beta_4(double b4) {
    beta4_ = b4;
  }
  // Set the nucleus polar angle.
  inline void set_polar_angle(double theta) {
    polar_theta_ = theta;
  }
  // Set the nucleus azimuthal angle.
  inline void set_azimuthal_angle(double phi) {
    azimuthal_phi_ = phi;
  }

 private:
  // Deformation parameters.
  double beta2_ = 0.0;
  double beta4_ = 0.0;
  // Nucleus orientation (initial profile in xz plane).
  double polar_theta_ = 0.0;
  double azimuthal_phi_ = 0.0;
};

}

#endif /* SRC_INCLUDE_DEFORMEDNUCLEUS_H_ */