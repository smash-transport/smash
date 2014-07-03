/*
 *    Andy's Deformed Nucleus Class Header
 */
#ifndef SRC_INCLUDE_DEFORMEDNUCLEUS_H_
#define SRC_INCLUDE_DEFORMEDNUCLEUS_H_

#include "configuration.h"
#include "nucleus.h"
#include "threevector.h"

namespace Smash {

class DeformedNucleus : public Nucleus {
 public:
  DeformedNucleus();

  // Return the deformed Woods-Saxon probability for
  // the current parameters.
  double deformed_woods_saxon(double r, double cosx) const;

  // Deformed Woods-Saxon sampling routine.
  void deformed_distribute_nucleon(ThreeVector& vec) const;

  // Sets the positions of the nuclei inside nucleus A.
  virtual void arrange_nucleons();

  // Sets the parameters of the Woods-Saxon distribution
  // according to the current mass number.
  //
  // Ref. for deformation parameters is "Nuclear Ground-State Masses and Deformations" 
  // by P. Möller, J. R. Nix, W. D. Myers, and W. J. Swiatecki.
  // Corrections to deformation parameter beta2 in Uranium come from arxiv:nuclth/0506088 
  // by A. Kuhlman and U. Heinz.
  // For finite nucleon size corrections to the nuclear density and radius, see ref.
  // arxiv:0904.4080 [nucl-th] by T. Hirano and Y. Nara for copper and gold, and 
  // arxiv:1010.6222 [nucl-th] by T. Hirano, P. Huovinen, and Y. Nara for uranium.
  // Currently misplaced reference for Lead.
  virtual size_t determine_nucleus();

  // Set parameters for Woods-Saxon by hand using the configuration file.
  // \see Nucleus::manual_nucleus
  virtual void manual_nucleus(bool is_projectile, Configuration &modus_cfg);

  // Shifts the nucleus to correct impact parameter and z displacement.
  // \see Nucleus::shift
  virtual void shift(bool is_projectile, double initial_z_displacement,
                     double x_offset, float simulation_time);

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

 private:
  // Maximum nucleon radius.
  double r_max_ = 0.0;
  // Deformation parameters.
  double beta2_ = 0.0;
  double beta4_ = 0.0;
  // Nucleus orientation (initial profile in xz plane).
  double nucleus_polar_angle_ = 0.0;
  double nucleus_azimuthal_angle_ = 0.0;
};

}

#endif /* SRC_INCLUDE_DEFORMEDNUCLEUS_H_ */