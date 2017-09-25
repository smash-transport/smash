/*
 *    Copyright (c) 2014-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#include "include/deformednucleus.h"

#include <cmath>
#include <stdexcept>

#include "include/angles.h"
#include "include/configuration.h"
#include "include/constants.h"
#include "include/fourvector.h"
#include "include/particledata.h"
#include "include/random.h"
#include "include/threevector.h"

namespace Smash {

/*!\Userguide
 * \page input_deformed_nucleus_ Deformed Nucleus
 *
 * \li \key Beta_2 (double, optional, default = 0.0):\n
 * The deformation coefficient for the spherical harmonic Y_2_0 in the
 * beta decomposition of the nuclear radius in the deformed woods-saxon
 * distribution.
 * \li \key Beta_4 (double, optional, default = 0.0):\n
 * The deformation coefficient for the spherical harmonic Y_4_0.
 * \li \key Saturation_Density (double, optional, default = .168)\n
 * The normalization coefficient in the Woods-Saxon distribution,
 * needed here (and not in nucleus) due to the accept/reject sampling used.
 * Default is
 * given as the infinite nuclear matter value.
 * \li \key Theta (double, optional): \n
 * The polar angle by which to rotate the nucleus.
 * \li \key Phi (double, optional):\n
 * The azimuthal angle by which to rotate the nucleus.
 */

DeformedNucleus::DeformedNucleus(const std::map<PdgCode, int> &particle_list,
                                 int nTest)
    : Nucleus(particle_list, nTest) {}

DeformedNucleus::DeformedNucleus(Configuration &config, int nTest)
    : Nucleus(config, nTest) {}

double DeformedNucleus::deformed_woods_saxon(double r, double cosx) const {
  // Return the deformed woods-saxon calculation
  // at the given location for the current system.
  return Nucleus::get_saturation_density() /
         (1 +
          std::exp(r -
                   Nucleus::get_nuclear_radius() *
                       (1 + beta2_ * y_l_0(2, cosx) + beta4_ * y_l_0(4, cosx)) /
                       Nucleus::get_diffusiveness()));
}

ThreeVector DeformedNucleus::distribute_nucleon() const {
  double a_radius;
  Angles a_direction;
  // Set a sensible max bound for radial sampling.
  double radius_max =
      Nucleus::get_nuclear_radius() / Nucleus::get_diffusiveness() +
      Nucleus::get_nuclear_radius() * Nucleus::get_diffusiveness();

  // Sample the distribution.
  do {
    a_direction.distribute_isotropically();
    a_radius = radius_max * std::cbrt(Random::canonical());  // sample r**2 dr
  } while (Random::canonical() >
           deformed_woods_saxon(a_radius, a_direction.costheta()) /
               Nucleus::get_saturation_density());

  // Update (x, y, z).
  return a_direction.threevec() * a_radius;
}

void DeformedNucleus::set_parameters_automatic() {
  // Initialize the inherited attributes.
  Nucleus::set_parameters_automatic();

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

  // Polar angle
  if (config.has_value({"Theta"})) {
    set_polar_angle(static_cast<double>(config.take({"Theta"})));
  }
  // Azimuth
  if (config.has_value({"Phi"})) {
    set_azimuthal_angle(static_cast<double>(config.take({"Phi"})));
  }
}

void DeformedNucleus::rotate() {
  for (auto &particle : *this) {
    // Rotate every vector by the nuclear azimuth phi and polar angle
    // theta (the Euler angles). This means applying the matrix for a
    // rotation of phi about z, followed by the matrix for a rotation
    // theta about the rotated x axis. The third angle psi is 0 by symmetry.
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

}  // namespace Smash
