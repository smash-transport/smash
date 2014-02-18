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
  // PHYSICS:
  //    First, calculate the masses of both nuclei MA and MB.
  //    Second: gamma*beta = sqrt( s-(MA+MB)**2 )/sqrt(4*MA*MB)
  //    Third: Get momenta from p_i = +-gamma*beta*mass_i
  //
  //    That is if a sqrt(s)_total is given and we want calculations to
  //    happen in frame of equal velocity.
  //
  //    if sqrt(s)_(ij) is given (i and j being the indices of two
  //    particles), then
  //    s_tot = (s_ij - mi**2 - mj**2) * MA*MB/(mi*mj) + MA**2 + MB**2.
  //
  //    One could also think of defining sqrt(s)_NN as the collision of
  //    two particles whose mass is mA = MA/NA, where NA is the number
  //    of particles in nucleus A, and correspondingly for mB.
  //
  //    Momenta/positions
  //    First, initialize nuclei at (0,0,0) at rest
  //    Second, boost nuclei
  //    Third, shift them so that they barely touch each other
  //    Fourth, set the time when they /will/ touch to 0.0.
  projectile_.auto_set_masses(particles);
  target_.auto_set_masses(particles);
  float mass_projec = projectile_.mass();
  float mass_target = target_.mass();
  float mass_1 = particles->particle_type(pdg_sNN_1_).mass();
  float mass_2 = particles->particle_type(pdg_sNN_2_).mass();
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
  // shift the nuclei along the z-axis so that they barely touch and
  // along the x-axis to get the right impact parameter:
  // projectile_.shift(projectile_.z_max_, impact_/2.0);
  // target_.shift(target_.z_min_, impact_/2.0);
}

void NucleusModus::sample_impact(const bool s, float min, float max) {
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
