/*
 *    Andy's Deformed Nucleus Class File
 */
#include "include/deformednucleus.h"
#include "include/threevector.h"
#include "include/fourvector.h"
#include "include/constants.h"
#include "include/random.h"
#include "include/particledata.h"

#include <cmath>
#include <stdexcept>

namespace Smash {

DeformedNucleus::DeformedNucleus() {}

double DeformedNucleus::deformed_woods_saxon(const double& r, const double& cosx) {
  // Return the deformed woods-saxon calculation
  // at the given location for the current system.
  return saturation_density_ / (1 + std::exp(r - initial_radius_ * (1 + beta2_ * 
         y_l_0(2, cosx) + beta4_ * y_l_0(4, cosx)) / Nucleus::get_diffusiveness()));
}

void DeformedNucleus::deformed_distribute_nucleon(ThreeVector& vec) {
  // !!! Efficiency.
  // !!! Angles object?
  
  // Sample the distribution.
  double radius;
  double radius_max = 20; // !!! Initial value?
  double polar_angle;

  do {
    radius = Random::uniform(0.0, radius_max);
    polar_angle = Random::uniform(0.0, twopi / 2);

  } while (Random::canonical() > deformed_woods_saxon(radius, 
           std::cos(polar_angle)));

  double azimuthal_angle = Random::uniform(0.0, twopi);
  double sinx = std::sin(polar_angle);

  // Update (x, y, z).
  vec = ThreeVector(radius * sinx * std::cos(azimuthal_angle), radius *
        sinx * std::sin(azimuthal_angle), radius * std::cos(polar_angle));

}

void DeformedNucleus::arrange_nucleons() {
  for (auto i = Nucleus::begin(); i != Nucleus::end(); i++) { 
    // Using a deformed W.S., get both the radial 
    // position and solid angle for the nucleon.
    ThreeVector pos;
    deformed_distribute_nucleon(pos);

    // Set the position of the current nucleon.
    i->set_position(FourVector(0.0, pos));

    // Update max/min bounds of nucleus
    // !!! Need 3D coords.
  }
}

void DeformedNucleus::determine_nucleus() {

  // Establish the current system.
  size_t mass_number = Nucleus::number_of_particles();

  // Uranium
  if (mass_number == 238){
    // Hirano, Huovinen, Nara - Correction.
    Nucleus::set_diffusiveness(0.44);
    DeformedNucleus::set_initial_radius(6.86);
    // Previous values. !!! Unsure.
    // Nucleus::set_diffusiveness(0.54);
    // DeformedNucleus::set_initial_radius(6.86);

    DeformedNucleus::set_saturation_density(0.166);

    // Kuhlman, Heinz - Correction.
    DeformedNucleus::set_deformation_params(0.28, 0.093);
    // Moeller et. al. - Previous.
    // DeformedNucleus::set_deformation_params(0.215, 0.093);
  }  
  // Lead !!! Reference?
    else if (mass_number == 208) {
    Nucleus::set_diffusiveness(0.44);
    DeformedNucleus::set_initial_radius(6.67);
    DeformedNucleus::set_saturation_density(0.161);
    DeformedNucleus::set_deformation_params(0.0, 0.0);

  }
  // Gold 
    else if (mass_number == 197) {
    // Hirano, Nara - Corrections.
    Nucleus::set_diffusiveness(0.44);
    DeformedNucleus::set_initial_radius(6.42);
    // Previous values.
    // Nucleus::set_diffusiveness(0.535);
    // DeformedNucleus::set_initial_radius(6.38);

    DeformedNucleus::set_saturation_density(0.1695);
    DeformedNucleus::set_deformation_params(-0.131, -0.031);
  }
  // Copper
    else if (mass_number == 63){
    // Hirano, Nara - Corrections.
    Nucleus::set_diffusiveness(0.50);
    DeformedNucleus::set_initial_radius(4.28);    
    // Previous values.
    // Nucleus::set_diffusiveness(0.5977);
    // DeformedNucleus::set_initial_radius(4.20641);

    DeformedNucleus::set_saturation_density(0.1686);  
    DeformedNucleus::set_deformation_params(0.162, -0.006);
  } else {
    throw std::domain_error("Mass number not listed in DeformedNucleus::determine_nucleus.");
  }
}

double DeformedNucleus::y_l_0(const int l, const double& cosx) const {
  if (l == 2) {
    return (1./4) * std::sqrt(10/twopi) * (3. * cosx * cosx - 1);
  } else if (l == 4) {
    return (3./16) * std::sqrt(2/twopi) * (35. * std::pow(cosx, 4) - 30. * cosx * cosx + 3);
  } else {
    throw std::domain_error("Not a valid angular momentum quantum number in DeformedNucleus::y_l_0.");
  }
}

}