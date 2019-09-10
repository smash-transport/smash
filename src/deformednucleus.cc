/*
 *    Copyright (c) 2014-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#include "smash/deformednucleus.h"

#include <cmath>
#include <stdexcept>

#include "smash/angles.h"
#include "smash/configuration.h"
#include "smash/constants.h"
#include "smash/fourvector.h"
#include "smash/particledata.h"
#include "smash/random.h"
#include "smash/threevector.h"

namespace smash {

/*!\Userguide
 * \page projectile_and_target Projectile and Target
 *
 * \li \key Beta_2 (double, optional):\n
 * The deformation coefficient for the spherical harmonic Y_2_0 in the
 * beta decomposition of the nuclear radius in the deformed woods-saxon
 * distribution. \n
 *
 * \li \key Beta_4 (double, optional):\n
 * The deformation coefficient for the spherical harmonic Y_4_0. \n
 *
 * \li \key Theta (double, optional): \n
 * The polar angle by which to rotate the nucleus. \n
 *
 * \li \key Phi (double, optional):\n
 * The azimuthal angle by which to rotate the nucleus.
 *
 * \li \key Random_Rotation (bool, optional, default = false):\n
 * Determines whether the created nucleus object should be randomly rotated in
 * space.
 *
 * \n
 */

// For readability and layout issues parts of the customnucleus userguide
// needs to be located here. The example configuration can however be found in
// customnucleus.cc

/*!\Userguide
 * \page projectile_and_target Projectile and Target
 *
 * - \key Custom: \n
 *    \li \key File_Directory (path, required if \key Custom exists): \n
 *
 * The directory where the external list with the nucleon configurations
 * is located. Make sure to use an absolute path.\n
 *
 * \li \key File_Name (string, required if \key Custom exists): \n
 * The file name of the external list with the nucleon configurations.
 *
 */

/*!\Userguide
* \n
* Example: Configuring a deformed nucleus
* --------------
* To configure a fixed target heavy-ion collision with deformed nuclei, whose
* spherical deformation is explicitly declared, it can be done according to the
* following example. For explanatory (and not physics) reasons,
* the projectile's Woods Saxon distribution is initialized automatically and
* it's spherical deformation manually, while the target nucleus is configured
* just the opposite.
*\verbatim
Modi:
    Collider:
        Projectile:
            Particles:    {2212: 29, 2112: 34}
            Deformed:
                # Manually set deformation parameters
                Automatic: False
                Beta_2: 0.1
                Beta_4: 0.3
                Theta: 0.8
                Phi: 0.02
        Target:
            Particles:    {2212: 29, 2112: 34}
            # manually set woods saxon parameters
            Saturation_Density: 0.1968
            Diffusiveness: 0.8
            Radius: 2.0
            Deformed:
                # Automatically set deformation parameters
                Automatic: True
                # Randomly rotate nucleus
                Random_Rotation: True
        E_kin: 1.2
        Calculation_Frame: "fixed target"
\endverbatim
*/

DeformedNucleus::DeformedNucleus(const std::map<PdgCode, int> &particle_list,
                                 int nTest)
    : Nucleus(particle_list, nTest) {}

DeformedNucleus::DeformedNucleus(Configuration &config, int nTest,
                                 bool auto_deformation)
    : Nucleus(config, nTest) {
  if (auto_deformation) {
    set_deformation_parameters_automatic();
  } else {
    set_deformation_parameters_from_config(config);
  }
}

ThreeVector DeformedNucleus::distribute_nucleon() {
  double a_radius;
  Angles a_direction;
  // Set a sensible maximum bound for radial sampling.
  double radius_max =
      Nucleus::get_nuclear_radius() / Nucleus::get_diffusiveness() +
      Nucleus::get_nuclear_radius() * Nucleus::get_diffusiveness();

  // Sample the distribution.
  do {
    a_direction.distribute_isotropically();
    // sample r**2 dr
    a_radius = radius_max * std::cbrt(random::canonical());
  } while (random::canonical() >
           nucleon_density(a_radius, a_direction.costheta()) /
               Nucleus::get_saturation_density());

  // Update (x, y, z) positions.
  return a_direction.threevec() * a_radius;
}

void DeformedNucleus::set_deformation_parameters_automatic() {
  // Set the deformation parameters
  // reference for U, Pb, Au, Cu: \iref{Moller:1993ed}
  // reference for Zr and Ru: \iref{Schenke:2019ruo}
  switch (Nucleus::number_of_particles()) {
    case 238:  // Uranium
      set_beta_2(0.28);
      set_beta_4(0.093);
      break;
    case 208:  // Lead
      set_beta_2(0.0);
      set_beta_4(0.0);
      break;
    case 197:  // Gold
      set_beta_2(-0.131);
      set_beta_4(-0.031);
      break;
    case 63:  // Copper
      set_beta_2(0.162);
      set_beta_4(-0.006);
      break;
    case 96: {
      size_t n_protons = Nucleus::number_of_protons();
      if (n_protons == 40) {  // Zirconium
        set_beta_2(0.0);
        set_beta_4(0.0);
        break;
      } else if (n_protons == 44) {  // Ruthenium
        set_beta_2(0.158);
        set_beta_4(0.0);
        break;
      } else {
        throw std::domain_error(
            "Number of protons for nuclei with mass number A = 96 does not "
            "match that of Zirconium or Ruthenium. The deformation parameters "
            "for additional isobars are currently not implemented."
            " Please specify at least \"Beta_2\" and \"Beta_4\" "
            "manually and set \"Automatic: False.\" ");
        break;
      }
      break;
    }
    default:
      throw std::domain_error(
          "Mass number not listed for automatically setting deformation "
          "parameters. Please specify at least \"Beta_2\" and \"Beta_4\" "
          "manually and set \"Automatic: False.\" ");
  }

  // Set a random nuclear rotation.
  nuclear_orientation_.distribute_isotropically();
}

void DeformedNucleus::set_deformation_parameters_from_config(
    Configuration &config) {
  // Deformation parameters.
  if (config.has_value({"Deformed", "Beta_2"})) {
    set_beta_2(static_cast<double>(config.take({"Deformed", "Beta_2"})));
  }
  if (config.has_value({"Deformed", "Beta_4"})) {
    set_beta_4(static_cast<double>(config.take({"Deformed", "Beta_4"})));
  }
  if (config.has_value({"Deformed", "Theta"})) {
    if (config.has_value({"Deformed", "Random_Rotation"}) &&
        config.take({"Deformed", "Random_Rotation"})) {
      throw std::domain_error(
          "Random rotation of nuclei is activated although"
          " theta is provided. Please specify only either of them. ");
    } else {
      set_polar_angle(static_cast<double>(config.take({"Deformed", "Theta"})));
    }
  }
  if (config.has_value({"Deformed", "Phi"})) {
    if (config.has_value({"Deformed", "Random_Rotation"}) &&
        config.take({"Deformed", "Random_Rotation"})) {
      throw std::domain_error(
          "Random rotation of nuclei is activated although"
          " phi is provided. Please specify only either of them. ");
    } else {
      set_azimuthal_angle(
          static_cast<double>(config.take({"Deformed", "Phi"})));
    }
  }
  if (config.take({"Deformed", "Random_Rotation"}, false)) {
    // Randomly generate euler angles for theta and phi. Psi needs not be
    // assigned, as the nucleus objects are symmetric with respect to psi.
    Nucleus::random_euler_angles();
    set_azimuthal_angle(euler_phi_);
    set_polar_angle(euler_theta_);
  }
}

void DeformedNucleus::rotate() {
  for (auto &particle : *this) {
    /* Rotate every vector by the nuclear azimuth phi and polar angle
     * theta (the Euler angles). This means applying the matrix for a
     * rotation of phi about z, followed by the matrix for a rotation
     * theta about the rotated x axis. The third angle psi is 0 by symmetry.*/
    ThreeVector three_pos = particle.position().threevec();
    three_pos.rotate(nuclear_orientation_.phi(), nuclear_orientation_.theta(),
                     0.);
    particle.set_3position(three_pos);
  }
}

/**
 * Spherical harmonics Y_2_0 and Y_4_0.
 * \param[in] l Angular momentum value (2 and 4 are supported)
 * \param[in] cosx Cosine of the polar angle
 * \return Value of the corresponding spherical harmonic
 * \throws domain_error if unsupported l is encountered
 */
static double y_l_0(int l, double cosx) {
  if (l == 2) {
    return (1. / 4) * std::sqrt(5 / M_PI) * (3. * (cosx * cosx) - 1);
  } else if (l == 4) {
    return (3. / 16) * std::sqrt(1 / M_PI) *
           (35. * (cosx * cosx) * (cosx * cosx) - 30. * (cosx * cosx) + 3);
  } else {
    throw std::domain_error(
        "Not a valid angular momentum quantum number in y_l_0.");
  }
}

double DeformedNucleus::nucleon_density(double r, double cosx) {
  return Nucleus::get_saturation_density() /
         (1 + std::exp((r - Nucleus::get_nuclear_radius() *
                                (1 + beta2_ * y_l_0(2, cosx) +
                                 beta4_ * y_l_0(4, cosx))) /
                       Nucleus::get_diffusiveness()));
}

}  // namespace smash
