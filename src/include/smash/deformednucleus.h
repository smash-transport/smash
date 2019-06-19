/*
 *    Copyright (c) 2014-2019
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

/**
 * DeformedNucleus: Child of nucleus for deformed nuclei.
 *
 * All options from the nucleus will still apply. The deformed nucleus adds
 * new or updated features which are outlined below.
 */

class DeformedNucleus : public Nucleus {
 public:
  /**
   * Constructor for DeformedNucles which takes a particle list and the number
   * of testparticles. This constructor is only used for testing purposes.
   * \param[in] particle_list Map with PDGCode and number of particles which
   * make up the nucleus \param[in] nTest number of testparticles
   */
  DeformedNucleus(const std::map<PdgCode, int> &particle_list, int nTest);
  /**
   * Constructor for DeformedNucleus, that needs the configuration parameters
   * from the inputfile and the number of testparticles
   * \param[in] config contains the parameters from the inputfile on the
   * numbers of particles with a certain PDG code
   * \param[in] nTest number of testparticles
   * \param[in] auto_deformation whether or not deformation parameters
   * should be set automatically
   */
  DeformedNucleus(Configuration &config, int nTest, bool auto_deformation);

  /**
   * Return the deformed Woods-Saxon probability for the given position.
   *
   * \param[in] r The radius at which to sample
   * \param[in] cosx The cosine of the polar angle at which to sample
   * \return The Woods-Saxon probability
   */
  double deformed_woods_saxon(double r, double cosx) const;

  /**
   * Deformed Woods-Saxon sampling routine.
   *
   * \return Spatial position from uniformly sampling
   * the deformed woods-saxon distribution
   */
  ThreeVector distribute_nucleon() override;

  /**
   * Sets the deformation parameters of the radius according to the current
   * mass number.
   *
   * The deformation parameters are taken from \iref{Moller:1993ed}.
   * Corrections to the deformation parameter beta2 in Uranium come from
   * \iref{Kuhlman:2005ts}. For finite nucleon size corrections to the nuclear
   * density and radius for copper and gold, see \iref{Hirano:2009ah},
   * and \iref{Hirano:2010jg} for uranium.
   */
  void set_deformation_parameters_automatic();

  /**
   * Set parameters for spherical deformation of the nucleus from the values
   * specified in the configuration file.
   *
   * \param config The configuration for the deformation of this nucleus
   *        (projectile or target).
   */
  void set_deformation_parameters_from_config(Configuration &config);

  /**
   * Rotates the nucleus according to members nucleus_polar_angle_
   * and nucleus_azimuthal_angle_ and updates nucleon positions.
   */
  void rotate() override;

  /**
   * Does not allow to generate Fermi-momenta for a deformed nucleus.
   * \throws domain_error if this function is ever called
   */
  void generate_fermi_momenta() override;

  /**
   * Spherical harmonics Y_2_0 and Y_4_0.
   * \param[in] l Angular momentum value (2 and 4 are supported)
   * \param[in] cosx Cosine of the polar angle
   * \return Value of the corresponding spherical harmonic
   * \throws domain_error if unsupported l is encountered
   */
  double y_l_0(int l, double cosx) const;

  /**
   * Set deformation coefficient for Y_2_0.
   * \param[in] b2 deformation coefficient for l=2
   */
  inline void set_beta_2(double b2) { beta2_ = b2; }
  /**
   * Set deformation coefficient for Y_4_0.
   * \param[in] b4 deformation coefficient for l=4
   */
  inline void set_beta_4(double b4) { beta4_ = b4; }
  /**
   * Set the nucleus polar angle.
   * \param[in] theta Polar angle of position inside nucleus
   */
  inline void set_polar_angle(double theta) {
    nuclear_orientation_.set_theta(theta);
  }
  /**
   * Set the nucleus azimuthal angle.
   * \param[in] phi Azimuthal angle of position inside nucleus
   */
  inline void set_azimuthal_angle(double phi) {
    nuclear_orientation_.set_phi(phi);
  }

 private:
  /// Deformation parameter for angular momentum l=2.
  double beta2_ = 0.0;
  /// Deformation parameter for angular momentum l=4.
  double beta4_ = 0.0;
  /**
   * Nucleus orientation (initial profile in xz plane) in terms of
   * a pair of angles (theta, phi)
   */
  Angles nuclear_orientation_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_DEFORMEDNUCLEUS_H_
