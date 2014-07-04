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

void DeformedNucleus::deformed_distribute_nucleon(ThreeVector& vec) const {
  // !!! Efficiency.
  // !!! radius_max?

  double a_radius;
  double radius_max = 2*Nucleus::get_nuclear_radius();
  Angles a_direction;

  // Sample the distribution.
  do {
    a_direction.distribute_isotropically();
    a_radius = Random::uniform(0.0, radius_max);

  } while (Random::canonical() > deformed_woods_saxon(a_radius, 
           a_direction.costheta()));

  // Update (x, y, z).
  vec = ThreeVector(a_radius * a_direction.x(), a_radius * 
                    a_direction.y(), a_radius * a_direction.z());

}

void DeformedNucleus::arrange_nucleons() {
  for (auto i = Nucleus::begin(); i != Nucleus::end(); i++) { 
    // Using a deformed W.S., get the radial 
    // position and solid angle for the nucleon.
    ThreeVector pos;
    deformed_distribute_nucleon(pos);

    // Set the position of the current nucleon.
    i->set_position(FourVector(0.0, pos));

    // Update the radial bound of the nucleus.
    double r_tmp = pos.abs();
    r_max_ = (r_tmp > r_max_) ? r_tmp : r_max_;
  }
}

size_t DeformedNucleus::determine_nucleus() {
  // Establish the current system.
  size_t mass_number = Nucleus::determine_nucleus();
  // Set deformation parameters.
  // Uranium
  if (mass_number == 238) {
    // Moeller et. al. - Default.
    DeformedNucleus::set_beta_2(0.215);
    DeformedNucleus::set_beta_4(0.093);
    // Kuhlman, Heinz - Correction.
    // DeformedNucleus::set_beta_2(0.28);
    // DeformedNucleus::set_beta_4(0.093);
  }  
  // Lead
    else if (mass_number == 208) {
      DeformedNucleus::set_beta_2(0.0);
      DeformedNucleus::set_beta_4(0.0);
  }
  // Gold 
    else if (mass_number == 197) {
      DeformedNucleus::set_beta_2(-0.131);
      DeformedNucleus::set_beta_4(-0.031);
  }
  // Copper
    else if (mass_number == 63) {
      DeformedNucleus::set_beta_2(0.162);
      DeformedNucleus::set_beta_4(-0.006);
  } else {
      throw std::domain_error("Mass number not listed in DeformedNucleus::determine_nucleus.");
  }
  return mass_number;
}

void DeformedNucleus::manual_nucleus(bool is_projectile, Configuration &config) {
  // Regular nucleus parameters.
  Nucleus::manual_nucleus(is_projectile, config);

  const char * nucleus_type = is_projectile ? "Projectile" : "Target";
  // Deformation parameters.
  if (config.has_value({nucleus_type, "BETA_2"})) {
    set_beta_2(static_cast<double>(config.take({nucleus_type, "BETA_2"})));
  }
  if (config.has_value({nucleus_type, "BETA_4"})) {
    set_beta_4(static_cast<double>(config.take({nucleus_type, "BETA_4"})));
  }
}

void DeformedNucleus::shift(bool is_projectile, double initial_z_displacement,
                            double x_offset, float simulation_time) {
  // Determine the nucleus rotation.
  Angles nuclear_rotation;
  nuclear_rotation.distribute_isotropically();
  nucleus_polar_angle_ = nuclear_rotation.theta();
  nucleus_azimuthal_angle_ = nuclear_rotation.phi();

  // The amount to shift the z coordinates. If is_projectile, we shift
  // back by -r_max_, else we shift forward r_max_.
  double z_offset = is_projectile ? -r_max_ : r_max_;
  // In the current system, the nuclei may touch. We want them to be
  // a little apart, so we include a slightly bigger offset.
  z_offset += initial_z_displacement;

  // Move the nucleons to their new x, y, and z positions, and set the time.
  for (auto i = begin(); i != end(); i++) {
    FourVector this_position = i->position();

    // Get current position.
    double old_radius = std::sqrt(this_position.x1() * this_position.x1()
                           + this_position.x2() * this_position.x2()
                           + this_position.x3() * this_position.x3());
    double old_polar_angle = std::acos(this_position.x3()/old_radius);
    double old_azimuthal_angle = std::atan(this_position.x2()/this_position.x1());

    // Rotate
    this_position.set_x1(old_radius * std::sin(old_polar_angle + nucleus_polar_angle_)
                         * std::cos(old_azimuthal_angle + nucleus_azimuthal_angle_));
    this_position.set_x2(old_radius * std::sin(old_polar_angle + nucleus_polar_angle_) 
                          * std::sin(old_azimuthal_angle + nucleus_azimuthal_angle_));
    this_position.set_x3(old_radius * std::cos(old_polar_angle + nucleus_polar_angle_));

    // Translate
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
    return (3./16) * std::sqrt(1/M_PI) * (35. * (cosx * cosx) * (cosx * cosx)- 30. * (cosx * cosx) + 3);
  } else {
    throw std::domain_error("Not a valid angular momentum quantum number in DeformedNucleus::y_l_0.");
  }
}

} // namespace Smash