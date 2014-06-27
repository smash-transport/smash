/*
 *    Andy's Deformed Nucleus Class File
 */
#include "include/angles.h"
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

  double radius;
  double radius_max = 20;
  Angles dir;

  // Sample the distribution.
  do {
    dir.distribute_isotropically();
    radius = Random::uniform(0.0, radius_max);

  } while (Random::canonical() > deformed_woods_saxon(radius, 
           dir.costheta()));

  // Update (x, y, z).
  vec = ThreeVector(radius * dir.x(), radius * dir.y(), radius * dir.z());

}

void DeformedNucleus::arrange_nucleons() {
  for (auto i = Nucleus::begin(); i != Nucleus::end(); i++) { 
    // Using a deformed W.S., get the radial 
    // position and solid angle for the nucleon.
    ThreeVector pos;
    deformed_distribute_nucleon(pos);

    // Set the position of the current nucleon.
    i->set_position(FourVector(0.0, pos));

    // Update max/min bounds of nucleus
    // !!! Need 3D coords.
  }
}

size_t DeformedNucleus::determine_nucleus() {
  // Establish the current system.
  size_t mass_number = Nucleus::determine_nucleus();
  // Set deformation parameters.
  // Uranium
  if (mass_number == 238){
    // Moeller et. al. - Default.
    DeformedNucleus::set_deformation_params(0.215, 0.093);
    // Kuhlman, Heinz - Correction.
    // DeformedNucleus::set_deformation_params(0.28, 0.093);
  }  
  // Lead !!! Reference?
    else if (mass_number == 208) {
      DeformedNucleus::set_deformation_params(0.0, 0.0);
  }
  // Gold 
    else if (mass_number == 197) {
      DeformedNucleus::set_deformation_params(-0.131, -0.031);
  }
  // Copper
    else if (mass_number == 63){
      DeformedNucleus::set_deformation_params(0.162, -0.006);
  } else {
      throw std::domain_error("Mass number not listed in DeformedNucleus::determine_nucleus.");
  }
  return mass_number;
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

}