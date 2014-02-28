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

NucleusModus::NucleusModus(Configuration modus_config) {
  Configuration modus_cfg = modus_config["Nucleus"];
  sqrt_s_NN_ = modus_cfg.take({"SQRTSNN"});
  std::vector<int> sqrts_n = modus_cfg.take({"SQRTS_N"});
  pdg_sNN_1_ = sqrts_n[0];
  pdg_sNN_2_ = sqrts_n[1];
  // fill nuclei with particles
  std::map<int, int> pro = modus_cfg.take({"Projectile", "PARTICLES"});
  projectile_.fill_from_list(pro);
  std::map<int, int> tar = modus_cfg.take({"Target", "PARTICLES"});
  target_.fill_from_list(tar);
  // set softness of the nuclei if given (else take the default value)
  if (modus_cfg.has_value({"Projectile", "SOFTNESS"})) {
    projectile_.set_softness(static_cast<float>(
                modus_cfg.take({"Projectile", "SOFTNESS"})));
  }
  if (modus_cfg.has_value({"Target", "SOFTNESS"})) {
    target_.set_softness(static_cast<float>(
                modus_cfg.take({"Target", "SOFTNESS"})));
  }

  // Impact paramter setting: Either "VALUE", "RANGE" or "MAX".
  if (modus_cfg.has_value({"Impact", "VALUE"})) {
    impact_ = modus_cfg.take({"Impact", "VALUE"});
  } else {
    bool sampling_quadratically = true;
    float min = 0.0;
    float max = 0.0;
    // if not value, we need to take a look at the sampling:
    if (modus_cfg.has_value({"Impact", "SAMPLE"})) {
      std::string sampling_method = modus_cfg.take({"IMPACT", "SAMPLE"});
      if (sampling_method.compare(0, 6, "uniform") == 0) {
        sampling_quadratically = false;
      }
    }
    if (modus_cfg.has_value({"Impact", "RANGE"})) {
       std::vector<float> range = modus_cfg.take({"Impact", "RANGE"});
       min = range.at(0);
       max = range.at(1);
    }
    if (modus_cfg.has_value({"Impact", "MAX"})) {
       min = 0.0;
       max = modus_cfg.take({"Impact", "MAX"});
    }
    sample_impact(sampling_quadratically, min, max);
  }
  if (modus_cfg.has_value({"INITIAL_DISTANCE"})) {
    initial_z_displacement_ = modus_cfg.take({"INITIAL_DISTANCE"});
    // the displacement is half the distance (both nuclei are shifted
    // initial_z_displacement_ away from origin)
    initial_z_displacement_ /= 2.0;
  }
}

void NucleusModus::print_startup() {
  printf("Nucleus initialized:\n");
  printf("sqrt_s_NN = %g GeV (pairs of %d and %d)\n", sqrt_s_NN_,
                                                      pdg_sNN_1_, pdg_sNN_2_);
  printf("Impact parameter: %g fm\n", impact_);
  printf("Initial distance betw nuclei: %g fm\n", 2.0*initial_z_displacement_);
  printf("Projectile has been initialized with %lu particles\n",
                                                           projectile_.size());
  printf("Target     has been initialized with %lu particles\n",
                                                               target_.size());
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
  printf("Masses of Nuclei: %g GeV %g GeV\n", projectile_.mass(), target_.mass());
  printf("Radii of Nuclei: %g fm %g fm\n", projectile_.nuclear_radius(), target_.nuclear_radius());
  float mass_1, mass_2;
  // set the masses used in sqrt_sNN. mass1 corresponds to the
  // projectile.
  if (pdg_sNN_1_ != 0) {
    // If PDG Code is given, use mass of this particle type.
    mass_1 = particles->particle_type(pdg_sNN_1_).mass();
  } else if(projectile_.size() > 0) {
    // else, use average mass of a particle in that nucleus
    mass_1 = projectile_.mass()/projectile_.size();
  } else {
    throw "Projectile nucleus empty!\n";
  }
  // same logic for mass2 and target as for projectile directly above.
  if (pdg_sNN_2_ != 0) {
    mass_2 = particles->particle_type(pdg_sNN_2_).mass();
  } else if(target_.size() > 0) {
    mass_2 = target_.mass()/target_.size();
  } else {
    throw "Target nucleus empty!\n";
  }
  double s_NN = sqrt_s_NN_*sqrt_s_NN_;
  if (s_NN < (mass_1 + mass_2)*(mass_1 + mass_2)) {
    throw "Error in input: sqrt(s_NN) is smaller than masses.";
  }
  float total_mandelstam_s = (s_NN - mass_1*mass_1 - mass_2*mass_2)
                             * mass_projec*mass_target
                             / (mass_1*mass_2)
                           + mass_projec*mass_projec + mass_target*mass_target;
  float velocity_squared = (total_mandelstam_s - (mass_projec + mass_target)
                                               * (mass_projec + mass_target))
                         / (total_mandelstam_s - (mass_projec - mass_target)
                                               * (mass_projec - mass_target));
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
  // -initial_z_displacement_/sqrt(velocity_squared).
  double simulation_time = -initial_z_displacement_/sqrt(velocity_squared);
  projectile_.shift(true, -initial_z_displacement_, +impact_/2.0
                                                  , simulation_time);
  target_.shift(false,    +initial_z_displacement_, -impact_/2.0
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
