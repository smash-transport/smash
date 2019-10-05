/*
 *    Copyright (c) 2014-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#include "smash/deformednucleus.h"

#include <cmath>
#include <map>
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
 * \li \key Orientation
 * \n
 * Determines the orientation of the nucleus by rotations
 * which are performed about the axes of a coordinate system
 * that is fixed with respect to the nucleus and whose axes
 * are parallel to those of the computational frame before the first rotation.
 * Note that the nucleus is first rotated by phi and then by theta.
 *    - \key Phi (double, optional, default = 0):\n
 * The angle by which to rotate the nucleus about the z-axis.
 *    - \key Theta (double, optional, default = pi/2): \n
 * The angle by which to rotate the nucleus about the rotated x-axis.
 *    - \key Random_Rotation (bool, optional, default = false):\n
 * Determines whether the created nucleus object should be randomly rotated in
 * space. \n
 * \n
 */

// For readability and layout issues parts of the customnucleus userguide
// needs to be located here. The example configuration can however be found in
// customnucleus.cc

/*!\Userguide
 * \page projectile_and_target Projectile and Target
 *
 *   - \key Custom: \n
 * \li \key File_Directory (path, required if \key Custom exists): \n
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
* its spherical deformation manually, while the target nucleus is configured
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
                Orientation:
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
                Orientation:
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

  if (config.has_value({"Deformed", "Orientation"})) {
    Configuration subconfig = config["Deformed"]["Orientation"];
    set_orientation_from_config(subconfig);
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
  bool listed = 0;
  std::map<int, std::string> A_map = {
      {238, "Uranium"}, {208, "Lead"}, {197, "Gold"}, {63, "Copper"}};
  std::map<std::string, std::string> Z_map = {
      {"Uranium", "92"}, {"Lead", "82"}, {"Gold", "79"}, {"Copper", "29"}};
  int A = Nucleus::number_of_particles();
  int Z = Nucleus::number_of_protons();
  switch (A) {
    case 238:  // Uranium
      if (Z == 92) {
        set_beta_2(0.28);
        set_beta_4(0.093);
      } else {
        listed = true;
      }
      break;
    case 208:  // Lead
      if (Z == 82) {
        set_beta_2(0.0);
        set_beta_4(0.0);
      } else {
        listed = true;
      }
      break;
    case 197:  // Gold
      if (Z == 79) {
        set_beta_2(-0.131);
        set_beta_4(-0.031);
      } else {
        listed = true;
      }
      break;
    case 63:  // Copper
      if (Z == 29) {
        set_beta_2(0.162);
        set_beta_4(-0.006);
      } else {
        listed = true;
      }
      break;
    case 96:
      if (Z == 40) {  // Zirconium
        set_beta_2(0.0);
        set_beta_4(0.0);
      } else if (Z == 44) {  // Ruthenium
        set_beta_2(0.158);
        set_beta_4(0.0);
      } else if (Z == 44) {  // Ruthenium
        set_beta_2(0.158);
        set_beta_4(0.0);
      } else {
        throw std::domain_error(
            "Number of protons for nuclei with mass number A = 96 does not "
            "match that of Zirconium or Ruthenium. The deformation parameters "
            "for additional isobars are currently not implemented."
            " Please specify at least \"Beta_2\" and \"Beta_4\" "
            "manually and set \"Automatic: False.\" ");
      }
      break;
    default:
      throw std::domain_error(
          "Mass number not listed for automatically setting deformation "
          "parameters. Please specify at least \"Beta_2\" and \"Beta_4\" "
          "manually and set \"Automatic: False.\" ");
  }
  if (listed) {
    throw std::domain_error("Mass number is listed under " + A_map[A] +
                            " but the proton "
                            "number of " +
                            std::to_string(Z) +
                            " does not match "
                            "its " +
                            Z_map[A_map[A]] +
                            " protons."
                            "Please specify at least \"Beta_2\" and \"Beta_4\" "
                            "manually and set \"Automatic: False.\" ");
  }
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
}

void DeformedNucleus::set_orientation_from_config(
    Configuration &orientation_config) {
  // Read in orientation if provided, otherwise, the defaults are
  // theta = pi/2, phi = 0, as declared in the angles class

  if (orientation_config.has_value({"Theta"})) {
    if (orientation_config.has_value({"Random_Rotation"}) &&
        orientation_config.take({"Random_Rotation"})) {
      throw std::domain_error(
          "Random rotation of nuclei is activated although"
          " theta is provided. Please specify only either of them. ");
    } else {
      set_polar_angle(static_cast<double>(orientation_config.take({"Theta"})));
    }
  }

  if (orientation_config.has_value({"Phi"})) {
    if (orientation_config.has_value({"Random_Rotation"}) &&
        orientation_config.take({"Random_Rotation"})) {
      throw std::domain_error(
          "Random rotation of nuclei is activated although"
          " phi is provided. Please specify only either of them. ");
    } else {
      set_azimuthal_angle(
          static_cast<double>(orientation_config.take({"Phi"})));
    }
  }

  if (orientation_config.take({"Random_Rotation"}, false)) {
    random_rotation_ = true;
  }
}

void DeformedNucleus::rotate() {
  if (random_rotation_) {
    // Randomly generate euler angles for theta and phi. Psi needs not be
    // assigned, as the nucleus objects are symmetric with respect to psi.
    Nucleus::random_euler_angles();
    set_azimuthal_angle(euler_phi_);
    set_polar_angle(euler_theta_);
  }
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

double y_l_0(int l, double cosx) {
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
