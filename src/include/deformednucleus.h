/*
 *    Copyright (c) 2014-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_DEFORMEDNUCLEUS_H_
#define SRC_INCLUDE_DEFORMEDNUCLEUS_H_

#include <map>

#include "angles.h"
#include "configuration.h"
#include "forwarddeclarations.h"
#include "nucleus.h"
#include "threevector.h"

namespace smash {

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
 * Options added by NucleusModus go in the "Modi"→"Nucleus"->"a nucleus" section
 * of the
 * configuration, where "a nucleus" is either projectile or target.
 * All options from the nucleus will still apply. The deformed nucleus adds new
 * or updated
 * features which are outlined below.
 *
 * The following deformed nucleus directives are understood:
 * -------------
 */

class DeformedNucleus : public Nucleus {
 public:
  DeformedNucleus(const std::map<PdgCode, int> &particle_list, int nTest);
  DeformedNucleus(Configuration &config, int nTest);

  /** Return the deformed Woods-Saxon probability for the given position.
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
  ThreeVector distribute_nucleon() const override;

  /** Sets the deformation parameters of the Woods-Saxon distribution
   * according to the current mass number.
   *
   * The deformation parameters are taken from \iref{Moller:1993ed}.
   * Corrections to the deformation parameter beta2 in Uranium come from
   * \iref{Kuhlman:2005ts}. For finite nucleon size corrections to the nuclear
   * density and radius, see \iref{Hirano:2009ah} for copper and gold,
   * and \iref{Hirano:2010jg} for uranium.
   */
  void set_parameters_automatic() override;

  /** Set parameters for Woods-Saxon by hand using the configuration file.
   * \see Nucleus::set_parameters_from_config
   */
  void set_parameters_from_config(Configuration &config) override;

  /** Rotates the nucleus according to members nucleus_polar_angle_
   * and nucleus_azimuthal_angle_ and updates nucleon positions.
   */
  void rotate() override;

  /**
   * Does not allow to generate Fermi-momenta for a deformed nucleus.
   **/
  void generate_fermi_momenta() override;

  /// Spherical harmonics Y_2_0 and Y_4_0.
  double y_l_0(int l, double cosx) const;

  /// Set deformation coefficient for Y_2_0.
  inline void set_beta_2(double b2) { beta2_ = b2; }
  /// Set deformation coefficient for Y_4_0.
  inline void set_beta_4(double b4) { beta4_ = b4; }
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

}  // namespace smash

#endif  // SRC_INCLUDE_DEFORMEDNUCLEUS_H_
