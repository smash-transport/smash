/*
 *    Copyright (c) 2014-2018
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
 * Additional Parameters to specify the shape of the modified Woods-Saxon
 * distribution, necessary to configure the deformation of the nucleus if
 * \key Deformed: True and \key Automatic: False.: \n
 *
 * \li \key Beta_2 (double, optional):\n
 * The deformation coefficient for the spherical harmonic Y_2_0 in the
 * beta decomposition of the nuclear radius in the deformed woods-saxon
 * distribution. \n
 *
 * \li \key Beta_4 (double, optional):\n
 * The deformation coefficient for the spherical harmonic Y_4_0. \n
 *
 * \li \key Saturation_Density (double, optional, default = .168)\n
 * The normalization coefficient in the Woods-Saxon distribution,
 * needed here (and not in nucleus) due to the accept/reject sampling used.
 * Default is given as the infinite nuclear matter value. \n
 *
 * \li \key Theta (double, optional): \n
 * The polar angle by which to rotate the nucleus. \n
 *
 * \li \key Phi (double, optional):\n
 * The azimuthal angle by which to rotate the nucleus.
 *
 * \n
 * Example: Configuring a deformed nucleus
 * --------------
 * To configure a fixed target heavy-ion collision with deformed nuclei, whose
 * deformation is explicitly declared, it can be done according to the following
 * example:
 *\verbatim
 Modi:
     Collider:
         Projectile:
             Particles:    {2212: 29, 2112 :34}
             Deformed: True
             Automatic: False
             Beta_2: 0.0
             Beta_4: 0.0
             Saturation_Density: 0.168
             Theta: 0.0
             Phi: 0.0
         Target:
             Particles:    {2212: 29, 2112 :34}
             Deformed: True
             Automatic: False
             Beta_2: 0.0
             Beta_4: 0.0
             Saturation_Density: 0.168
             Theta: 0.0
             Phi: 0.0
         E_kin: 1.2
         Calculation_Frame: "fixed target"
 \endverbatim
 */

DeformedNucleus::DeformedNucleus(const std::map<PdgCode, int> &particle_list,
                                 int nTest)
    : Nucleus(particle_list, nTest) {}

DeformedNucleus::DeformedNucleus(Configuration &config, int nTest)
    : Nucleus(config, nTest) {}

double DeformedNucleus::deformed_woods_saxon(double r, double cosx) const {
  return Nucleus::get_saturation_density() /
         (1 + std::exp(r - Nucleus::get_nuclear_radius() *
                               (1 + beta2_ * y_l_0(2, cosx) +
                                beta4_ * y_l_0(4, cosx)) /
                               Nucleus::get_diffusiveness()));
}

ThreeVector DeformedNucleus::distribute_nucleon() const {
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
           deformed_woods_saxon(a_radius, a_direction.costheta()) /
               Nucleus::get_saturation_density());

  // Update (x, y, z) positions.
  return a_direction.threevec() * a_radius;
}

void DeformedNucleus::set_parameters_automatic() {
  // Initialize the inherited attributes.
  Nucleus::set_parameters_automatic();
  /// \todo In issue #4743 the corrected values should be optional
  // Set the deformation parameters.
  switch (Nucleus::number_of_particles()) {
    case 238:  // Uranium
      // default: Moeller
      set_beta_2(0.215);
      set_beta_4(0.093);
      // correction: Kuhlman, Heinz
      // set_beta_2(0.28);
      // set_beta_4(0.093);
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
    default:
      throw std::domain_error(
          "Mass number not listed in"
          " DeformedNucleus::set_parameters_automatic.");
  }

  // Set a random nuclear rotation.
  nuclear_orientation_.distribute_isotropically();
}

void DeformedNucleus::set_parameters_from_config(Configuration &config) {
  // Inherited nucleus parameters.
  Nucleus::set_parameters_from_config(config);
  // Deformation parameters.
  if (config.has_value({"Beta_2"})) {
    set_beta_2(static_cast<double>(config.take({"Beta_2"})));
  }
  if (config.has_value({"Beta_4"})) {
    set_beta_4(static_cast<double>(config.take({"Beta_4"})));
  }
  // Saturation density (normalization for accept/reject sampling)
  if (config.has_value({"Saturation_Density"})) {
    Nucleus::set_saturation_density(
        static_cast<double>(config.take({"Saturation_Density"})));
  }
  if (config.has_value({"Theta"})) {
    set_polar_angle(static_cast<double>(config.take({"Theta"})));
  }
  if (config.has_value({"Phi"})) {
    set_azimuthal_angle(static_cast<double>(config.take({"Phi"})));
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

void DeformedNucleus::generate_fermi_momenta() {
  throw std::domain_error(
      "Fermi momenta currently not implemented"
      " for a deformed nucleus.");
}

double DeformedNucleus::y_l_0(int l, double cosx) const {
  if (l == 2) {
    return (1. / 4) * std::sqrt(5 / M_PI) * (3. * (cosx * cosx) - 1);
  } else if (l == 4) {
    return (3. / 16) * std::sqrt(1 / M_PI) *
           (35. * (cosx * cosx) * (cosx * cosx) - 30. * (cosx * cosx) + 3);
  } else {
    throw std::domain_error(
        "Not a valid angular momentum quantum number in"
        "DeformedNucleus::y_l_0.");
  }
}

}  // namespace smash
