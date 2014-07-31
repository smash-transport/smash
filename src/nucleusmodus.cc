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
#include <tuple>
#include <utility>
 
#include "include/nucleusmodus.h"
#include "include/angles.h"
#include "include/configuration.h"
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

  // Get the reference frame for the collision calculation.
  if (modus_cfg.has_value({"CALCULATION_FRAME"})) {
    frame_ = modus_cfg.take({"CALCULATION_FRAME"});
  }

  // Decide which type of nucleus: deformed or not (default).
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
  if (projectile_->size() < 1) {
    throw NucleusEmpty("Input Error: Projectile nucleus is empty.");
  }
  std::map<PdgCode, int> tar = modus_cfg.take({"Target", "PARTICLES"});
  target_->fill_from_list(tar, params.testparticles);
  if (target_->size() < 1) {
    throw NucleusEmpty("Input Error: Target nucleus is empty.");
  }

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

  // Get the total nucleus-nucleus collision energy. Since there is 
  // no meaningful choice for a default energy, we require the user to
  // give one (and only one) energy input from the available options.
  int energy_input = 0;
  float mass_projec = projectile_->mass();
  float mass_target = target_->mass();
  // Option 1: Center of mass energy.
  if (modus_cfg.has_value({"SQRTSNN"})) {
      float sqrt_s_NN = modus_cfg.take({"SQRTSNN"});
      // Note that \f$\sqrt{s_{NN}}\f$ is different for neutron-neutron and
      // proton-proton collisions (because of the different masses).
      // Therefore,representative particles are needed to specify which
      // two particles' collisions have this \f$\sqrt{s_{NN}}\f$. The vector
      // specifies a pair of PDG codes for the two particle species we want to use.
      // The default is otherwise the average nucleon mass for each nucleus.
      PdgCode id_1 = 0, id_2 = 0;
      if (modus_cfg.has_value({"SQRTS_REPS"})) {
        std::vector<PdgCode> sqrts_reps = modus_cfg.take({"SQRTS_REPS"});
        id_1 = sqrts_reps[0];
        id_2 = sqrts_reps[1];
      }
      float mass_1, mass_2;
      if (id_1 != 0) {
        // If PDG Code is given, use mass of this particle type.
        mass_1 = ParticleType::find(id_1).mass();
      } else {
        // else, use average mass of a particle in that nucleus
        mass_1 = projectile_->mass()/projectile_->size();
      }
      if (id_2 != 0) {
        mass_2 = ParticleType::find(id_2).mass();
      } else {
        mass_2 = target_->mass()/target_->size();
      } 
      // Check that input satisfies the lower bound (everything at rest).
      if (sqrt_s_NN < mass_1 + mass_2) {
        throw ModusDefault::InvalidEnergy(
                "Input Error: sqrt(s_NN) is smaller than masses:\n"
                + std::to_string(sqrt_s_NN) + " GeV < "
                + std::to_string(mass_1) + " GeV + "
                + std::to_string(mass_2) + " GeV.");
      }
      // Set the total nucleus-nucleus collision energy.
      total_s_= (sqrt_s_NN * sqrt_s_NN - mass_1 * mass_1 - mass_2 * mass_2)
                * mass_projec * mass_target / (mass_1 * mass_2)
                + mass_projec * mass_projec + mass_target * mass_target;
      energy_input++;
  }
  // Option 2: Energy of the projectile nucleus (target at rest).
  if (modus_cfg.has_value({"E_LAB"})) {
      int e_lab = modus_cfg.take({"E_LAB"});
      // Check that energy is nonnegative.
      if (e_lab < 0) {
        throw ModusDefault::InvalidEnergy("Input Error: E_LAB must be nonnegative.");
      }
      // Set the total nucleus-nucleus collision energy.
      total_s_ = (mass_projec * mass_projec) + (mass_target * mass_target)
                 + 2 * e_lab * mass_target;
      energy_input++;
  }
  // Option 3: Momentum of the projectile nucleus (target at rest).
  if (modus_cfg.has_value({"P_LAB"})) {
      int p_lab = modus_cfg.take({"P_LAB"});
      // Check upper bound (projectile mass).
      if (p_lab * p_lab > mass_projec * mass_projec) {
        throw ModusDefault::InvalidEnergy(
                "Input Error: P_LAB squared is greater than projectile mass squared: \n"
                + std::to_string(p_lab * p_lab) + " GeV > "
                + std::to_string(mass_projec * mass_projec) + " GeV");
      }
      // Set the total nucleus-nucleus collision energy.
      total_s_ = (mass_projec * mass_projec) + (mass_target * mass_target)
                  + 2 * (mass_projec * mass_projec) * mass_target
                  / sqrt((mass_projec * mass_projec) - (p_lab * p_lab));
      energy_input++;
  }
  if (energy_input != 1){
    throw std::domain_error("Input Error: Redundant or nonexistant collision energy.");
  } 

  // Impact parameter setting: Either "VALUE", "RANGE", or "MAX".
  // Unspecified means 0 impact parameter.
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
    // Sample impact parameter distribution.
    sample_impact(sampling_quadratically, min, max);
  }

  // Look for user-defined initial separation between nuclei.
  if (modus_cfg.has_value({"INITIAL_DISTANCE"})) {
    initial_z_displacement_ = modus_cfg.take({"INITIAL_DISTANCE"});
    // the displacement is half the distance (both nuclei are shifted
    // initial_z_displacement_ away from origin)
    initial_z_displacement_ /= 2.0;
  }
}

void NucleusModus::print_startup() {
  printf("Nucleus initialized:\n");
  printf("S (nucleus-nucleus) = %g GeV \n", total_s_);
  printf("Impact parameter = %g fm\n", impact_);
  printf("Initial distance between nuclei: %g fm\n", 2.0*initial_z_displacement_);
  printf("Projectile initialized with %zu particles (%zu test particles)\n",
          projectile_->number_of_particles(), projectile_->size());
  printf("Target initialized with %zu particles (%zu test particles)\n",
          target_->number_of_particles(), target_->size());
  printf("Masses (projectile, target): %g GeV %g GeV\n",
          projectile_->mass(), target_->mass());
  printf("Radii (projectile, target): %g fm %g fm\n",
          projectile_->get_nuclear_radius(), target_->get_nuclear_radius());
}

float NucleusModus::initial_conditions(Particles *particles,
                                      const ExperimentParameters&) {
  // Populate the nuclei with appropriately distributed nucleons.
  // If deformed, this includes rotating the nucleus.
  projectile_->arrange_nucleons();
  target_->arrange_nucleons();

  // Use the total mandelstam variable to get the frame-dependent velocity for
  // each nucleus. Position 1 is projectile, position 2 is target.
  double v1, v2;
  std::tie(v1, v2) = get_velocities(total_s_, projectile_->mass(), target_->mass());

  // If velocities are too close to 1 for our calculations, throw an exception.
  if (almost_equal(std::abs(1.0 - v1), 0.0)
      || almost_equal(std::abs(1.0 - v2), 0.0)) {
    throw std::domain_error("Found velocity equal to 1 in nucleusmodus::initial_conditions.\n"
                            "Consider using the center of velocity reference frame.");
  }

  // Shift the nuclei into starting positions. Keep the pair separated
  // in z by some small distance and shift in x by the impact parameter. 
  // (Projectile is chosen to hit at positive x.)
  // After shifting, set the time component of the particles to
  // -initial_z_displacement_/average_velocity.
  float avg_velocity = sqrt(v1 * v1 
                            + v2 * v2);
  float simulation_time = -initial_z_displacement_ / avg_velocity;
  projectile_->shift(true, -initial_z_displacement_, +impact_/2.0,
                     simulation_time);
  target_->shift(false, initial_z_displacement_, -impact_/2.0,
                 simulation_time);

  // Boost the nuclei to the appropriate velocity.
  projectile_->boost(v1);
  target_->boost(v2);

  // Put the particles in the nuclei into code particles.
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

std::pair<double, double> NucleusModus::get_velocities(float s, float m1, float m2) {
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
  return std::make_pair(v1, v2);
}

}  // namespace Smash
