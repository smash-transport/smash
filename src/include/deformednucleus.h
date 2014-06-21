/*
 *    Andy's Deformed Nucleus Class Header
 */
#ifndef SRC_INCLUDE_DEFORMEDNUCLEUS_H_
#define SRC_INCLUDE_DEFORMEDNUCLEUS_H_

#include "threevector.h"
#include "nucleus.h"

namespace Smash {

class DeformedNucleus : public Nucleus {
 public:
  DeformedNucleus();

  // Return the deformed Woods-Saxon probability for
  // the current parameters.
  double deformed_woods_saxon(const double& r, const double& cosx);

  // Deformed Woods-Saxon sampling routine.
  void deformed_distribute_nucleon(ThreeVector& vec);

  // Sets the positions of the nuclei inside nucleus A.
  virtual void arrange_nucleons();

  // Sets the parameters of the Woods-Saxon distribution
  // according to the current mass number.
  //
  // Ref. for deformation parameters is "Nuclear Ground-State Masses and Deformations" 
  // by P. Möller, J. R. Nix, W. D. Myers, and W. J. Swiatecki.  Corrections to
  // deformation parameter beta2 in Uranium come from arxiv:nuclth/0506088 by
  // A. Kuhlman and U. Heinz.
  // For finite nucleon size corrections to the nuclear density and radius, see ref.
  // arxiv:0904.4080 [nucl-th] by T. Hirano and Y. Nara for copper and gold, and 
  // arxiv:1010.6222 [nucl-th] by T. Hirano, P. Huovinen, and Y. Nara for uranium.
  // Currently unsure of reference for Lead.
  void determine_nucleus();

  // Spherical harmonics Y_2_0 and Y_4_0.
  double y_l_0(const int l, const double& cosx) const;

  void set_saturation_density(const double& density){
    saturation_density_ = density;
  }
  
  void set_initial_radius(const double& radius){
    initial_radius_ = radius;
  }
  
  void set_deformation_params(const double& b2, const double& b4) {
    beta2_ = b2;
    beta4_ = b4;
  }

  void set_manual_deformation(bool x) {
    manual_deformation_ = x;
  }

 private:
  double manual_deformation_ = false;
  double saturation_density_ = 1.0;
  double initial_radius_ = 1.0;
  double beta2_ = 0.0;
  double beta4_ = 0.0;

};

}

#endif /* SRC_INCLUDE_DEFORMEDNUCLEUS_H_ */