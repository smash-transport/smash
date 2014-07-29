/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include "boost/tuple/tuple.hpp"
 
#include "include/nucleusmodus.h"
#include "include/angles.h"
#include "include/configuration.h"
//#include "include/constants.h"
#include "include/experimentparameters.h"
#include "include/numerics.h"
#include "include/outputroutines.h"
#include "include/particles.h"
#include "include/pdgcode.h"
#include "include/random.h"

namespace Smash {

NucleusModus::NucleusModus(Configuration modus_config,
                           const ExperimentParameters &params) {
  Configuration modus_cfg = modus_config["Nucleus"];

  // TODO: Allow user to select different energy inputs (energy/mom. of beam, etc.)
  // Get energy input.
  sqrt_s_NN_ = modus_cfg.take({"SQRTSNN"});
  std::vector<PdgCode> sqrts_n = modus_cfg.take({"SQRTS_N"});
  pdg_sNN_1_ = sqrts_n[0];
  pdg_sNN_2_ = sqrts_n[1];

  // Get the reference frame for the calculation.
  // TODO: Improve frame user interface.
  if (modus_cfg.has_value({"CALCULATION_FRAME"})) {
    frame_ = modus_cfg.take({"CALCULATION_FRAME"});
  }

  // Decide which type of nucleus (deformed or not).
  if (modus_cfg.has_value({"Projectile", "DEFORMED"}) &&
      modus_cfg.take({"Projectile", "DEFORMED"})) {
    projectile_ = new DeformedNucleus();
  } else {
    projectile_ = new Nucleus();
  }
  if (modus_cfg.has_value({"Target", "DEFORMED"}) &&
      modus_cfg.take({"Target", "DEFORMED"})) {
    target_ = new DeformedNucleus();
  } else {
    target_ = new Nucleus();
  }  

  // Fill nuclei with particles.
  std::map<PdgCode, int> pro = modus_cfg.take({"Projectile", "PARTICLES"});
  projectile_->fill_from_list(pro, params.testparticles);
  std::map<PdgCode, int> tar = modus_cfg.take({"Target", "PARTICLES"});
  target_->fill_from_list(tar, params.testparticles);

  // Ask to construct nuclei based on atomic number; otherwise, look
  // for the user defined values or take the default parameters.
  if (modus_cfg.has_value({"Projectile", "AUTOMATIC"}) && 
      modus_cfg.take({"Projectile", "AUTOMATIC"})) {
    projectile_->set_parameters_automatic();
  } else {
    projectile_->set_parameters_from_config(true, modus_cfg);
  }
  if (modus_cfg.has_value({"Target", "AUTOMATIC"}) &&
      modus_cfg.take({"Target", "AUTOMATIC"})) {
    target_->set_parameters_automatic();
  } else {
    target_->set_parameters_from_config(false, modus_cfg);
  }

  // Impact parameter setting: Either "VALUE", "RANGE", or "MAX".
  if (modus_cfg.has_value({"Impact", "VALUE"})) {
    impact_ = modus_cfg.take({"Impact", "VALUE"});
  } else {
    bool sampling_quadratically = true;
    float min = 0.0;
    float max = 0.0;
    // If impact is not supplied by value, inspect sampling parameters:
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
    // Get an impact parameter.
    sample_impact(sampling_quadratically, min, max);
  }

  // Determine the initial separation between nuclei.
  if (modus_cfg.has_value({"INITIAL_DISTANCE"})) {
    initial_z_displacement_ = modus_cfg.take({"INITIAL_DISTANCE"});
    // the displacement is half the distance (both nuclei are shifted
    // initial_z_displacement_ away from origin)
    initial_z_displacement_ /= 2.0;
  }
}

void NucleusModus::print_startup() {
  printf("Nucleus initialized:\n");
  printf("sqrt_s_NN = %g GeV (pairs of %s and %s)\n", sqrt_s_NN_,
         pdg_sNN_1_.string().c_str(), pdg_sNN_2_.string().c_str());
  printf("Impact parameter: %g fm\n", impact_);
  printf("Initial distance betwn nuclei: %g fm\n", 2.0*initial_z_displacement_);
  printf("Projectile initialized with %zu particles (%zu test particles)\n",
          projectile_->number_of_particles(), projectile_->size());
  printf("Target initialized with %zu particles (%zu test particles)\n",
          target_->number_of_particles(), target_->size());
}

float NucleusModus::initial_conditions(Particles *particles,
                                      const ExperimentParameters&) {
  float mass_projec = projectile_->mass();
  float mass_target = target_->mass();
  printf("Masses of Nuclei: %g GeV %g GeV\n", projectile_->mass(),
                                              target_->mass());
  printf("Radii of Nuclei: %g fm %g fm\n", projectile_->get_nuclear_radius(),
                                           target_->get_nuclear_radius());

  // Populate the nuclei with appropriately distributed nucleons.
  // If deformed, this includes rotating the nucleus.
  projectile_->arrange_nucleons();
  target_->arrange_nucleons();

  // set the masses used in sqrt_sNN. mass1 corresponds to the
  // projectile.
  float mass_1, mass_2;
  if (pdg_sNN_1_ != 0) {
    // If PDG Code is given, use mass of this particle type.
    mass_1 = ParticleType::find(pdg_sNN_1_).mass();
  } else if (projectile_->size() > 0) {
    // else, use average mass of a particle in that nucleus
    mass_1 = projectile_->mass()/projectile_->size();
  } else {
    throw NucleusEmpty("Projectile nucleus is empty!");
  }
  if (pdg_sNN_2_ != 0) {
    mass_2 = ParticleType::find(pdg_sNN_2_).mass();
  } else if (target_->size() > 0) {
    mass_2 = target_->mass()/target_->size();
  } else {
    throw NucleusEmpty("Target nucleus is empty!");
  }

  // Check if energy is valid.
  double s_NN = sqrt_s_NN_ * sqrt_s_NN_;
  if (sqrt_s_NN_ < mass_1 + mass_2) {
    throw ModusDefault::InvalidEnergy(
                      "Error in input: sqrt(s_NN) is smaller than masses:\n"
                      + std::to_string(sqrt_s_NN_) + " GeV < "
                      + std::to_string(mass_1) + " GeV + "
                      + std::to_string(mass_2) + " GeV.");
  }
  // Calculate the lorentz invariant mandelstam variable for the total system.
  float total_mandelstam_s = (s_NN - mass_1 * mass_1 - mass_2 * mass_2)
                             * mass_projec * mass_target / (mass_1 * mass_2)
                             + mass_projec * mass_projec + mass_target * mass_target;
  // Use the total mandelstam variable to get the frame-denendent velocity for
  // each nucleus. Position 1 is projectile, position 2 is target.
  boost::tuple<double, double> velocities = get_velocities(total_mandelstam_s, 
                                                           mass_projec, mass_target);                           
  // If velocities are too close to 1 for our calculations, throw an exception.
  if (almost_equal(std::abs(1.0 - velocities.get<0>()), 0.0)
      || almost_equal(std::abs(1.0 - velocities.get<1>()), 0.0)) {
    throw std::domain_error("Found velocity equal to 1 in nucleusmodus::initial_conditions.\n"
                            "Consider using the center of velocity reference frame.");
  }

  // Shift the nuclei into starting positions.
  // Keep the pair separated in z by some small distance,
  // and shift in x by the impact parameter. (Projectile is 
  // chosen to hit at positive x.)
  // For regular nuclei, the shift is along the z-axis so that
  // the nuclei are 2*1 fm apart.
  // For deformed nuclei, movement is also along z but due to
  // geometry, initial separation may include extra space.
  // After shifting, set the time component of the particles to
  // -initial_z_displacement_/sqrt(velocity_squared).
  float avg_velocity = sqrt(velocities.get<0>() * velocities.get<0>() 
                            + velocities.get<1>() * velocities.get<1>());
  float simulation_time = -initial_z_displacement_ / avg_velocity;
  projectile_->shift(true, -initial_z_displacement_, +impact_/2.0,
                     simulation_time);
  target_->shift(false, initial_z_displacement_, -impact_/2.0,
                 simulation_time);

  // Boost the nuclei to the appropriate velocity.
  projectile_->boost(velocities.get<0>());
  target_->boost(velocities.get<1>());

  // Put the particles in the nuclei into particles.
  projectile_->copy_particles(particles);
  target_->copy_particles(particles);
  return simulation_time;
}

void NucleusModus::sample_impact(bool s, float min, float max) {
  if (s) {
    // quadratic sampling: Note that for min > max, this still yields
    // the correct distribution (only that canonical() = 0 then is the
    // upper end, not the lower).
    impact_ = sqrt(min*min + Random::canonical() * (max*max - min*min));
  } else {
    // linear sampling. Still, min > max works fine.
    impact_ = Random::uniform(min, max);
  }
}

boost::tuple<double, double> NucleusModus::get_velocities(float s, float m1, float m2) {
  double v1 = 0.0;
  double v2 = 0.0;
  // Frame dependent calculations of velocities. Assume v1 >= 0, v2 <= 0.
  switch (frame_) {
    case 1:  // Center of velocity.
        v1 = sqrt((s - (m1 + m2) * (m1 + m2))
                  / (s - (m1 - m2) * (m1 - m2)));
        v2 = - v1;
      break;
    case 2:  // Center of mass.
      {
        double A = (s -(m1 - m2) * (m1 - m2)) * (s -(m1 + m2) * (m1 + m2));
        double B = - 8 * (m1 * m1) * m1 * (m2 * m2) * m2 - ((m1 * m1) + (m2 * m2))
                   * (s - (m1 * m1) - (m2 * m2)) * (s - (m1 * m1) - (m2 * m2));
        double C = (m1 * m1) * (m2 * m2) * A;
        // Compute positive center of mass momentum.
        double abs_p = sqrt((-B - sqrt(B * B - 4 * A * C)) / (2 * A));
        v1 = abs_p / m1;
        v2 = -abs_p / m2;
      }
      break;
    case 3:  // Target at rest.
      v1 = sqrt(1 - 4 * (m1 * m1) * (m2 * m2) / ((s - (m1 * m1) - (m2 * m2)) 
                * (s - (m1 * m1) - (m2 * m2))));
      break;
    default:
      throw std::domain_error("Invalid reference frame in NucleusModus::get_velocities.");
  }
  return boost::make_tuple(v1, v2);
}

}  // namespace Smash
