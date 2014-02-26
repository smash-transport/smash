/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include <limits>

#include "include/nucleus.h"
#include "include/angles.h"

Nucleus::Nucleus() {}

float Nucleus::mass() const {
  float total_mass = 0.f;
  for (auto i = cbegin(); i != cend(); i++) {
    total_mass += sqrt(i->momentum().Dot());
  }
  return total_mass;
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
    i->set_position(0.0, x, r*dir.y(), z);
    // update maximal and minimal z values
    z_max_ = (z > z_max_) ? z : z_max_;
    z_min_ = (z < z_min_) ? z : z_min_;
    // update maximal and minimal x values
    x_max_ = (x > x_max_) ? x : x_max_;
    x_min_ = (x < x_min_) ? x : x_min_;
  }
}

void Nucleus::boost(const double& beta_squared) {
  double sign = beta_squared >= 0 ? 1 : -1;
  double beta_squared_abs = std::abs(beta_squared);
  double one_over_gamma = sqrt(1.0 - beta_squared_abs);
  //double gamma = 1.0/one_over_gamma;
  //double gammabeta = sign*sqrt(beta_squared_abs)*gamma;
  FourVector u_mu(1, 0.0, 0.0, sign*sqrt(beta_squared_abs));
  for (auto i = begin(); i != end(); i++) {
    // a real Lorentz Transformation would leave the particles at
    // different times here, which we would then have to propagate back
    // to equal times. Since we know the result, we can simply multiply
    // the z-value with 1/gamma.
    FourVector this_position = i->position();
    this_position.set_x3(this_position.x3() * one_over_gamma);
    i->set_position(this_position);
    // for momenta, though, we CAN do normal Lorentz Boosts, since we
    // *do* want to transform the zero-component (i.e., the energy).
    FourVector this_momentum = i->momentum();
    this_momentum = this_momentum.LorentzBoost(u_mu);
    i->set_momentum(this_momentum);
  }
  printf("Velocity vector: %12e %12e %12e %12e\n",
         u_mu.x0(), u_mu.x1(), u_mu.x2(), u_mu.x3()
        );
  // we also need to update z_max_:
  z_max_ *= one_over_gamma;
  return;
}

void Nucleus::fill_from_list(const std::map<int, int>& particle_list) {
  for (auto n = particle_list.cbegin(); n != particle_list.cend(); n++) {
    for (int i = 0; i < n->second; i++) {
      // append particle to list
      particles.push_back({});
      // set this particle's PDG code.
      particles.back().set_pdgcode(n->first);
    }
  }
}

void Nucleus::set_softness(const float& soft) {
  softness_ = soft;
}

void Nucleus::auto_set_masses(const Particles *external_particles) {
  for (auto p = begin(); p != end(); p++) {
    p->set_momentum(
      external_particles->particle_type(p->pdgcode()).mass(), 0.0, 0.0, 0.0);
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
    FourVector this_position = i->position();
    this_position.set_x3(this_position.x3() + z_offset);
    this_position.set_x1(this_position.x1() + x_offset);
    this_position.set_x0(simulation_time);
    i->set_position(this_position);
  }
}

void Nucleus::copy_particles(Particles* external_particles) {
  for (auto p = begin(); p != end(); p++) {
    external_particles->add_data(*p);
  }
}
