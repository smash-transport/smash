/*
 *    Copyright (c) 2014-2022,2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_SMASH_DEFORMEDNUCLEUS_H_
#define SRC_INCLUDE_SMASH_DEFORMEDNUCLEUS_H_

#include <map>

#include "angles.h"
#include "configuration.h"
#include "forwarddeclarations.h"
#include "nucleus.h"
#include "threevector.h"

namespace smash {

/**
 * Spherical harmonics Y_2_0, Y_2_2, Y_3_0 and Y_4_0.
 * \param[in] l Angular momentum value (2 and 4 are supported)
 * \param[in] m projection value (l = 2 and m = 2 are supported)
 * \param[in] cosx Cosine of the polar angle
 * \param[in] phi Azimuthal angle
 * \return Value of the corresponding spherical harmonic
 * \throws domain_error if unsupported l is encountered
 */
double y_l_m(int l, int m, double cosx, double phi);

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
  DeformedNucleus(
      const std::map<PdgCode, int> &particle_list, int nTest,
      SpinInteractionType spin_interaction_type = SpinInteractionType::Off);
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
   * \return the saturation density of the deformed_nucleus
   * \see saturation_density_
   */
  inline double get_saturation_density() const { return saturation_density_; }
  /**
   * Return the deformed Woods-Saxon probability density for the given position.
   * This corresponds to the nuclear density at the very same position.
   *
   * \param[in] r The radius at which to sample
   * \param[in] cosx The cosine of the polar angle at which to sample
   * \param[in] phi The azimuthal angle at which to sample
   * \return The Woods-Saxon density
   */
  double nucleon_density(double r, double cosx, double phi) const override;
  /**
   * Return the unnormalized deformed Woods-Saxon distribution for the given
   * position.
   *
   * \param[in] r The radius
   * \param[in] cosx The cosine of the polar angle
   * \param[in] phi The azimuthal angle
   * \return The unnormalized Woods-Saxon distribution
   */
  double nucleon_density_unnormalized(double r, double cosx,
                                      double phi) const override;
  /**
   * Return the integral over the azimuthal angle phi
   *
   * \param[in] r The radius
   * \param[in] cosx The cosine of the polar angle
   * \return The unnormalized Woods-Saxon distribution integrated over dphi
   */
  double integrant_nucleon_density_phi(double r, double cosx) const;
  /**
   * \return the normalized ground state density for the corresponding
   * Woods-Saxon parameter. This is done by integrating the Woods-Saxon
   * distribution and setting the normalization such that the integral of the
   * Woods-Saxon distribution yields the number of particles in the nucleus
   * \f$\int\rho(r)d^3r = N_{particles}\f$.
   *
   */
  double calculate_saturation_density() const override;
  /**
   * Set deformation coefficient for Y_2_0.
   * \param[in] b2 deformation coefficient for l=2
   */
  inline void set_beta_2(double b2) { beta2_ = b2; }
  /**
   * Set the triaxiality coefficient gamma for Y_2_0 and Y_2_2.
   * \param[in] ga triaxiality coefficient for l=2
   */
  inline void set_gamma(double ga) { gamma_ = ga; }
  /**
   * Set deformation coefficient for Y_3_0.
   * \param[in] b3 deformation coefficient for l=3
   */
  inline void set_beta_3(double b3) { beta3_ = b3; }
  /**
   * Set deformation coefficient for Y_4_0.
   * \param[in] b4 deformation coefficient for l=4
   */
  inline void set_beta_4(double b4) { beta4_ = b4; }
  /**
   * return the beta2 value.
   */
  inline double get_beta2() { return beta2_; }
  /**
   * return the beta3 value.
   */
  inline double get_beta3() { return beta3_; }
  /**
   * return the beta4 value.
   */
  inline double get_beta4() { return beta4_; }

 private:
  /// Deformation parameter for angular momentum l=2.
  double beta2_ = 0.0;
  /// Triaxiality parameter for angular momentum l=2.
  double gamma_ = 0.0;
  /// Deformation parameter for angular momentum l=3.
  double beta3_ = 0.0;
  /// Deformation parameter for angular momentum l=4.
  double beta4_ = 0.0;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_DEFORMEDNUCLEUS_H_
