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
#include "angles.h"

namespace Smash {

/** DeformedNucleus: Child of nucleus for deformed nuclei.
 *
 * To use this modus, choose (ex: for deformed projectile)
 * \code
 * Modi:
 *      Nucleus:
 *            Projectile:
 *                 Deformed: true
 * \endcode
 * in the configuration file.
 *
 * Options added by NucleusModus go in the "Modi"→"Nucleus"->"a nucleus" section of the
 * configuration, where "a nucleus" is either projectile or target.
 * All options from the nucleus will still apply. The deformed nucleus adds new or updated
 * features which are outlined below.
 *
 * The following deformed nucleus directives are understood:
 * -------------
 */
// !!USER:Input
/**
 * \if user
 * \page input_modi_nucleus_ Input Section Modi:Nucleus
 * \endif
 *
 * \li `AUTOMATIC:` Sets all necessary parameters based on the atomic number
 * of the input nucleus (true=automatic, false=manual, see additional directives).
 * \li `BETA_2` The deformation coefficient for the spherical harmonic Y_2_0 in the
 * \beta decomposition of the nuclear radius in the deformed woods-saxon distribution.
 * \li `BETA_4` The deformation coefficient for the spherical harmonic Y_4_0.
 * \li `SATURATION_DENISTY` The normalization coefficient in the Woods-Saxon distribution,
 * needed here (and not in nucleus) due to the accept/reject sampling used. Default is
 * given as the infinite nuclear matter value .168f.
 * \li `THETA` The polar angle by which to rotate the nucleus.
 * \li `PHI` The azimuthal angle by which to rotate the nucleus.
 **/
 // !!/USER:Input

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

  /** Set parameters for Woods-Saxon by hand using the configuration file.
   * \see Nucleus::set_parameters_from_config
   */
  virtual void set_parameters_from_config(bool is_projectile, Configuration &config);

  /** Rotates the nucleus according to members nucleus_polar_angle_
   * and nucleus_azimuthal_angle_ and updates nucleon positions.
   */
  virtual void rotate();

  /// Spherical harmonics Y_2_0 and Y_4_0.
  double y_l_0(int l, double cosx) const;
  
  /// Set deformation coefficient for Y_2_0.
  inline void set_beta_2(double b2) {
    beta2_ = b2;
  }
  /// Set deformation coefficient for Y_4_0.
  inline void set_beta_4(double b4) {
    beta4_ = b4;
  }
  /// Set the nucleus polar angle.
  inline void set_polar_angle(double theta) {
    nuclear_orientation_.set_theta(theta);
  }
  /// Set the nucleus azimuthal angle.
  inline void set_azimuthal_angle(double phi) {
    nuclear_orientation_.set_phi(phi);
  }

 private:
  /// Deformation parameter 2.
  double beta2_ = 0.0;
  /// Deformation parameter 4.
  double beta4_ = 0.0;
  /// Nucleus orientation (initial profile in xz plane).
  Angles nuclear_orientation_;

};

}

#endif /* SRC_INCLUDE_DEFORMEDNUCLEUS_H_ */