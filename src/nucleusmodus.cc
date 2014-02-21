/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#include <cmath>
#include <cstdlib>
#include <cstring>
// #include <list>

#include "include/nucleusmodus.h"
#include "include/angles.h"
#include "include/configuration.h"
#include "include/experimentparameters.h"
#include "include/outputroutines.h"

NucleusModus::NucleusModus(Configuration &config)
    : sqrt_s_NN_(config.take({"Nucleus", "SQRTSNN"})) {
  std::vector<int> sqrts_n = config.take({"Nucleus", "SQRTS_N"});
  pdg_sNN_1_ = sqrts_n[0];
  pdg_sNN_2_ = sqrts_n[1];
  // fill nuclei with particles
  std::map<int, int> pro = config.take({"Nucleus", "Projectile", "PARTICLES"});
  projectile_.fill_from_list(pro);
  std::map<int, int> tar = config.take({"Nucleus", "Target", "PARTICLES"});
  target_.fill_from_list(tar);
  // set softness of the nuclei if given (else take the default value)
  if (config.has_value({"Nucleus", "Projectile", "SOFTNESS"})) {
    projectile_.set_softness(static_cast<float>(
                config.take({"Nucleus", "Projectile", "SOFTNESS"})));
  }
  if (config.has_value({"Nucleus", "Target", "SOFTNESS"})) {
    target_.set_softness(static_cast<float>(
                config.take({"Nucleus", "Target", "SOFTNESS"})));
  }

  // Impact paramter setting: Either "VALUE", "RANGE" or "MAX".
  if (config.has_value({"Nucleus", "Impact", "VALUE"})) {
    impact_ = config.take({"Nucleus", "Impact", "VALUE"});
  } else {
    bool sampling_quadratically = true;
    float min = 0.0;
    float max = 0.0;
    // if not value, we need to take a look at the sampling:
    if (config.has_value({"Nucleus", "Impact", "SAMPLE"})) {
      std::string sampling_method = config.take({"Nucleus",
                                                 "IMPACT", "SAMPLE"});
      if (sampling_method.compare(0, 6, "linear") == 0) {
        sampling_quadratically = false;
      }
    }
    if (config.has_value({"Nucleus", "Impact", "RANGE"})) {
       std::vector<float> range = config.take({"Nucleus", "Impact", "RANGE"});
       min = range.at(0);
       max = range.at(1);
    }
    if (config.has_value({"Nucleus", "Impact", "MAX"})) {
       min = 0.0;
       max = config.take({"Nucleus", "Impact", "MAX"});
    }
    sample_impact(sampling_quadratically, min, max);
  }
}

/* initial_conditions - sets particle data for @particles */
void NucleusModus::initial_conditions(Particles *particles,
                                       const ExperimentParameters &) {
  // WHAT'S MISSING:
  //
  // Nuclei can be non-spherical. If they are, then they may be randomly
  // aligned. Excentricities and angles should be configurable.
  projectile_.auto_set_masses(particles);
  target_.auto_set_masses(particles);
  float mass_projec = projectile_.mass();
  float mass_target = target_.mass();
  float mass_1 = particles->particle_type(pdg_sNN_1_).mass();
  float mass_2 = particles->particle_type(pdg_sNN_2_).mass();
  if (sqrt_s_NN_ < (mass_1 + mass_2)*(mass_1 + mass_2)) {
    throw "Error in input: sqrt(s_NN) is smaller than masses.";
  }
  float total_mandelstam_s = (sqrt_s_NN_ - mass_1*mass_1 - mass_2*mass_2)
                             * mass_projec*mass_target
                             / (mass_1*mass_2)
                           + mass_projec*mass_projec + mass_target*mass_target;
  float velocity_squared = (total_mandelstam_s - (mass_projec+mass_target))
                         / (total_mandelstam_s - (mass_projec-mass_target));
  // populate the nuclei with appropriately distributed nucleons
  projectile_.arrange_nucleons();
  target_.arrange_nucleons();
  // boost the nuclei to the appropriate velocity (target is in opposite
  // direction!
  projectile_.boost(velocity_squared);
  target_.boost(-velocity_squared);
  // shift the nuclei along the z-axis so that they are 2*1 fm apart
  // touch and along the x-axis to get the right impact parameter.
  // Projectile hits at positive x.
  // Also, it sets the time of the particles to
  // -initial_z_displacement/sqrt(velocity_squared).
  double simulation_time = -initial_z_displacement/sqrt(velocity_squared);
  projectile_.shift(true, -initial_z_displacement, +impact_/2.0
                                                 , simulation_time);
  target_.shift(false,    +initial_z_displacement, -impact_/2.0
                                                 , simulation_time);
  // now, put the particles in the nuclei into particles.
  projectile_.copy_particles(particles);
  target_.copy_particles(particles);
}

void NucleusModus::sample_impact(const bool s, const float min,
                                               const float max) {
  // chi is the random number
  double chi = drand48();
  if (s) {
    // quadratic sampling: Note that for min > max, this still yields
    // the correct distribution (only that chi = 0 then is the upper
    // end, not the lower).
    impact_ = sqrt(min*min + chi * (max*max - min*min));
  } else {
    // linear sampling. Still, min > max works fine.
    impact_ = min + chi * (max - min);
  }
}
