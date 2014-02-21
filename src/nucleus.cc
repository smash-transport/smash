/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include <limits>

#include "include/nucleus.h"
#include "include/angles.h"

float Nucleus::mass() const {
  float total_mass = 0.f;
  for (std::map<int, ParticleData>::const_iterator i = cbegin();
       i != cend(); i++) {
    total_mass += particle_type(i->second.pdgcode()).mass();
  }
  return total_mass;
}

/// Nuclear radius is calculated with the proton radius times the third
/// root of the number of nucleons.
float inline Nucleus::nuclear_radius() const {
  return proton_radius_*pow(size(), 1./3.);
}

float Nucleus::distribution_nucleons() const {
  // softness_ zero or negative? Use hard sphere.
  if (softness_ < std::numeric_limits<float>::min()) {
    return nuclear_radius()*(pow(drand48(), 1./3.));
  }
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
    // get radial position of current nucleon:
    float r = distribution_nucleons();
    // get solid angle for current nucleon:
    Angles dir;
    dir.distribute_isotropically();
    double z = r*dir.z();
    double x = r*dir.x();
    // set position of current nucleon:
    i->second.set_position(0.0, x, r*dir.y(), z);
    // update maximal and minimal z values
    z_max_ = (z > z_max_) ? z : z_max_;
    z_min_ = (z < z_min_) ? z : z_min_;
    // update maximal and minimal x values
    x_max_ = (x > x_max_) ? x : x_max_;
    x_min_ = (x < x_min_) ? x : x_min_;
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
  // we also need to update z_max_:
  z_max_ *= one_over_gamma;
  return;
}

void Nucleus::fill_from_list(const std::map<int, int>& particle_list) {
  for (auto n = particle_list.cbegin(); n != particle_list.cend(); n++) {
    for (int i = 1; i < n->second; i++) {
      create(n->first);
    }
  }
}

void Nucleus::set_softness(const float& soft) {
  softness_ = soft;
}

void Nucleus::auto_set_masses(const Particles *particles) {
  for (auto p = begin(); p != end(); p++) {
    p->second.set_momentum(
      particles->particle_type(p->second.pdgcode()).mass(), 0.0, 0.0, 0.0);
  }
}

void Nucleus::shift(const bool is_projectile,
                    const double& initial_z_displacement,
                    const double& x_offset,
                    const double& simulation_time) {
  // amount to shift z value. If is_projectile, we shift to -z_max_,
  // else we shift to -z_min_ (z_min_ itself should be negative).
  double z_offset = is_projectile ? -z_max_ : -z_min_;
  // now, the nuclei would touch. We want them to be a little apart, so
  // we need a bigger offset.
  z_offset += initial_z_displacement;
  for (auto i = begin(); i != end(); i++) {
    FourVector this_position = i->second.position();
    this_position.set_x3(this_position.x3() + z_offset);
    this_position.set_x1(this_position.x1() + x_offset);
    this_position.set_x0(simulation_time);
    i->second.set_position(this_position);
  }
}

void Nucleus::copy_particles(Particles* particles) {
  for (auto p = begin(); p != end(); p++) {
    particles->add_data(p->second);
  }
}
