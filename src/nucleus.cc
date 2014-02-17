/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "include/nucleus.h"
#include "include/angles.h"

float Nucleus::get_mass() const {
  float mass = 0.f;
  for (std::map<int, ParticleData>::const_iterator i = cbegin();
       i != cend(); i++) {
    mass += particle_type(i->second.pdgcode()).mass();
  }
  return mass;
}

float inline Nucleus::nuclear_radius() const {
  return proton_radius_*pow(size(), 1./3.);
}

/* returns a length r with a relative probability of $$\frac{dN}{dr} =
 * \frac{r^2}{\exp\left(\frac{r-R}{d}\right) + 1}$$ where $d$ is a class
 * parameter and $R$ is the function parameter ''radius''. */
float Nucleus::distribution_nucleons() const {
  float radius_scaled = nuclear_radius()/softness_;
  float prob_range1 = 1.0;
  float prob_range2 = 3. / radius_scaled;
  float prob_range3 = 2. * prob_range2 / radius_scaled;
  float prob_range4 = 1. * prob_range3 / radius_scaled;
  float all_ranges = prob_range1 + prob_range2 + prob_range3 + prob_range4;
  float t;
  do {
    float which_range = drand48() * all_ranges - prob_range1;
    if (which_range < 0.0) {
      t = radius_scaled * (pow(drand48(), 1./3.) - 1.);
    } else {
      t = -log(drand48());
      if (which_range >= prob_range2) {
        t += -log(drand48());
        if (which_range >= prob_range2 + prob_range3) {
          t += -log(drand48());
        }
      }
    }
  } while (drand48() > 1./(1. + exp(-fabs(t)) ) );
  float position_scaled = t + radius_scaled;
  float position = position_scaled * softness_;
  return position;
}

void Nucleus::arrange_nucleons() {
  for (auto i = begin(); i != end(); i++) {
    float r = distribution_nucleons();
    Angles dir;
    dir.distribute_isotropically();
    i->second.set_position(0.0, r*dir.x(), r*dir.y(), r*dir.z());
  }
}

void Nucleus::boost(const double& beta_squared) {
  double one_over_gamma = sqrt(1.0 - beta_squared);
  double gamma = 1.0/one_over_gamma;
  double gammabeta = sqrt(beta_squared)*gamma;
  FourVector u_mu(gamma, 0.0, 0.0, gammabeta);
  for (auto i = begin(); i != end(); i++) {
    // a real Lorentz Transformation would leave the particles at
    // different times here, which we would then have to propagate back
    // to equal times. Since we know the result, we can simply multiply
    // the z-value with 1/gamma.
    FourVector this_position = i->second.position();
    this_position.set_x3(this_position.x3() * one_over_gamma);
    i->second.set_position(this_position);
    // for momenta, though, we CAN do normal Lorentz Boosts, since we
    // *do* want to transform the zero-component (i.e., the energy).
    FourVector this_momentum = i->second.momentum();
    this_momentum.LorentzBoost(u_mu);
    i->second.set_momentum(this_momentum);
  }
  // we also need to update max_z_:
  max_z_ *= one_over_gamma;
  return;
}
