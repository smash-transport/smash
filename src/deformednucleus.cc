/*
 *    Andy's Deformed Nucleus Class File
 */
#include "include/angles.h"
#include "include/configuration.h"
#include "include/constants.h"
#include "include/deformednucleus.h"
#include "include/fourvector.h"
#include "include/particledata.h"
#include "include/random.h"
#include "include/threevector.h"

#include <cmath>
#include <stdexcept>

namespace Smash {

DeformedNucleus::DeformedNucleus() {}

double DeformedNucleus::deformed_woods_saxon(double r, double cosx) const {
  // Return the deformed woods-saxon calculation
  // at the given location for the current system.
  return Nucleus::get_saturation_density() / (1 + std::exp(r - Nucleus::get_nuclear_radius() *
         (1 + beta2_ * y_l_0(2, cosx) + beta4_ * y_l_0(4, cosx)) / Nucleus::get_diffusiveness()));
}

ThreeVector DeformedNucleus::deformed_distribute_nucleon() const {
  double a_radius;
  Angles a_direction;
  // Set a sensible max bound for radial sampling.
  double radius_max = a_radius / Nucleus::get_diffusiveness() + 
                      a_radius * Nucleus::get_diffusiveness();

  // Sample the distribution.
  do {
    a_direction.distribute_isotropically();
    a_radius = Random::uniform(0.0, radius_max);

  } while (Random::canonical() > deformed_woods_saxon(a_radius,
           a_direction.costheta()));

  // Update (x, y, z).
  return a_direction.threevec() * a_radius;
}

void DeformedNucleus::arrange_nucleons() {
  for (auto i = Nucleus::begin(); i != Nucleus::end(); i++) {
    // Sampling a deformed W.S., get the radial
    // position and solid angle for the nucleon.
    ThreeVector pos = deformed_distribute_nucleon();

    // Set the position of the nucleon.
    i->set_position(FourVector(0.0, pos));

    // Update the radial bound of the nucleus.
    double r_tmp = pos.abs();
    r_max_ = (r_tmp > r_max_) ? r_tmp : r_max_;
  }
  // Recenter and rotate
  align_center();
  rotate();
}

void DeformedNucleus::set_parameters_automatic() {
  // Initialize the inherited attributes.
  Nucleus::set_parameters_automatic();

  // Set the deformation parameters.
  switch (Nucleus::number_of_particles()) {
    case 238:  // Uranium
      // Moeller et. al. - Default.
      set_beta_2(0.215);
      set_beta_4(0.093);
      // Kuhlman, Heinz - Correction.
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
      throw std::domain_error("Mass number not listed in DeformedNucleus::determine_nucleus.");
  }

  // Set a random nuclear rotation.
  Angles nuclear_rotation;
  nuclear_rotation.distribute_isotropically();
  set_polar_angle(nuclear_rotation.theta());
  set_azimuthal_angle(nuclear_rotation.phi());
}

void DeformedNucleus::set_parameters_from_config(bool is_projectile, Configuration &config) {
  // Inherited nucleus parameters.
  Nucleus::set_parameters_from_config(is_projectile, config);
  const char * nucleus_type = is_projectile ? "Projectile" : "Target";

  // Deformation parameters.
  if (config.has_value({nucleus_type, "BETA_2"})) {
    set_beta_2(static_cast<double>(config.take({nucleus_type, "BETA_2"})));
  }
  if (config.has_value({nucleus_type, "BETA_4"})) {
    set_beta_4(static_cast<double>(config.take({nucleus_type, "BETA_4"})));
  }

  // Polar angle
  if (config.has_value({nucleus_type, "THETA"})) {
    set_polar_angle(static_cast<double>(config.take({nucleus_type, "THETA"})));
  }
  // Azimuth
  if (config.has_value({nucleus_type, "PHI"})) {
    set_azimuthal_angle(static_cast<double>(config.take({nucleus_type, "PHI"})));
  }
}

void DeformedNucleus::rotate() {
  for (auto i = begin(); i != end(); i++) {
    FourVector this_position = i->position();
    double old_x = this_position.x1();
    double old_y = this_position.x2();
    double old_z = this_position.x3();

    // Rotate every vector by the private members azimuthal_phi_ and
    // polar_theta_ (Euler angles). This means applying the matrix for
    // a rotation of azimuthal_phi_ about z, followed by the matrix for
    // a rotation polar_theta_ about the rotated x axis.
    double cos_phi = std::cos(azimuthal_phi_);
    double sin_phi = std::sin(azimuthal_phi_);
    double cos_theta = std::cos(polar_theta_);
    double sin_theta = std::sin(polar_theta_);

    this_position.set_x1(cos_phi * old_x + sin_phi * old_y);
    this_position.set_x2(-cos_theta * sin_phi * old_x + cos_theta
                         * cos_phi * old_y + sin_theta * old_z);
    this_position.set_x3(sin_theta * sin_phi * old_x - sin_theta
                         * cos_phi * old_y + cos_theta * old_z);

    i->set_position(this_position);
  }
}

void DeformedNucleus::shift(bool is_projectile, double initial_z_displacement,
                            double x_offset, float simulation_time) {
  // The amount to shift the z coordinates. If is_projectile, we shift
  // back by -r_max_, else we shift forward r_max_.
  double z_offset = is_projectile ? -r_max_ : r_max_;
  // In the current system, the nuclei may touch. We want them to be
  // a little apart, so we include a slightly bigger offset.
  z_offset += initial_z_displacement;

  // Move the nucleons to the new x and z positions, and set the time.
  for (auto i = begin(); i != end(); i++) {
    FourVector this_position = i->position();
    this_position.set_x3(this_position.x3() + z_offset);
    this_position.set_x1(this_position.x1() + x_offset);
    this_position.set_x0(simulation_time);
    i->set_position(this_position);
  }
}

double DeformedNucleus::y_l_0(int l, double cosx) const {
  if (l == 2) {
    return (1./4) * std::sqrt(5/M_PI) * (3. * (cosx * cosx) - 1);
  } else if (l == 4) {
    return (3./16) * std::sqrt(1/M_PI) * (35. * (cosx * cosx) * (cosx * cosx) - 30. * (cosx * cosx) + 3);
  } else {
    throw std::domain_error("Not a valid angular momentum quantum number in DeformedNucleus::y_l_0.");
  }
}

} // namespace Smash
