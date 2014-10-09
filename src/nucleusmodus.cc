/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "include/nucleusmodus.h"
#include "include/angles.h"
#include "include/configuration.h"
#include "include/experimentparameters.h"
#include "include/logging.h"
#include "include/numerics.h"
#include "include/particles.h"
#include "include/pdgcode.h"
#include "include/random.h"

namespace Smash {
    
/*!\Userguide
 * \page input_modi_nucleus_ Nucleus
 *
 * Possible Incident Energies, only one can be given:
 * \key SQRTSNN (float, optional, no default): \n
 * Defines the energy of the collision as center-of-mass
 * energy in the collision of one participant each from both nuclei.
 * Optional: Since not all participants have the same mass, and hence
 * \f$\sqrt{s_{\rm NN}}\f$ is different for \f$NN\f$ = proton+proton and
 * \f$NN\f$=neutron+neutron, you can specify which \f$NN\f$-pair you
 * want this to refer to with `SQRTS_REPS`. This expects a vector of two
 * PDG Codes, e.g. `SQRTS_N: [2212, 2212]` for proton-proton.
 * Default is the average nucleon mass per nucleus.
 *
 * \key E_LAB (float, optional, no default): \n
 * Defines the energy of the collision by the initial energy in GeV? of
 * the projectile nucleus.  This assumes the target nucleus is at rest.
 *
 * \key P_LAB (float, optional, no default): \n
 * Defines the energy of the collision by the initial momentum
 * of the projectile nucleus in GeV?.  This assumes the target nucleus is at rest.
 *
 * \key CALCULATION FRAME (int, required, default = 1): \n
 * The frame in which the collision is calculated.\n
 * 1 - center of velocity frame\n
 * 2 - center of mass frame\n
 * 3 - fixed target frame
 *
 * \key Projectile: \n
 * Section for projectile nucleus. The projectile will
 * start at \f$z < 0\f$ and fly in positive \f$z\f$-direction, at \f$x
 * \ge 0\f$.
 *
 * \key Target: \n
 * Section for target nucleus. The target will start at \f$z
 * > 0\f$ and fly in negative \f$z\f$-direction, at \f$x \le 0\f$.
 *
 *
 * \key Projectile: and \key Target: \n
 * \li \key PARTICLES (int:int, int:int, required):\n
 * A map in which the keys are PDG codes and the
 * values are number of particles with that PDG code that should be in
 * the current nucleus. E.g.\ `PARTICLES: {2212: 82, 2112: 126}` for a
 * lead-208 nucleus (82 protons and 126 neutrons = 208 nucleons), and
 * `PARTICLES: {2212: 1, 2112: 1, 3122: 1}` for Hyper-Triton (one
 * proton, one neutron and one Lambda).
 * \li \key DEFORMED (bool, required): \n
 * true - deformed nucleus is initialized
 * false - spherical nucleus is initialized
 * \li \key AUTOMATIC (bool, required): \n
 * Whether or not to use default values based on the
 * current nucleus atomic number (true/false).
 *
 * Additional Woods-Saxon Parameters: \n
 * There are also many other
 * parameters for specifying the shape of the Woods-Saxon distribution,
 * and other nucleus specific properties. See NUCLEUS.CC and
 * DEFORMEDNUCLEUS.CC for more on these choices.
 *
 * \key Impact: \n
 * A section for the impact parameter (= distance of the two
 * straight lines that the center of masses of the nuclei travel on).
 *
 * \li \key VALUE (float, optional, optional, default = 0.f): fixed value for 
 * the impact parameter. No other \key Impact: directive is looked at.
 * \li \key SAMPLE (string, optional, default = quadratic sampling): \n
 * if \key uniform, use uniform sampling of the impact parameter 
 * (\f$dP(b) = db\f$). If else, use areal input sampling
 * (the probability of an input parameter range is proportional to the
 * area corresponding to that range, \f$dP(b) = b\cdot db\f$).
 * \li \key RANGE (float, float, optional, default = 0.0f):\n
 * A vector of minimal and maximal impact parameters
 * between which b should be chosen. (The order of these is not
 * important.)
 * \li \key MAX (float, optional, default = 0.0f): 
 * Like `RANGE: [0.0, MAX]`. Note that if both \key RANGE and
 * \key MAX are specified, \key MAX takes precedence.
 *
 * Note that there are no safeguards to prevent you from specifying
 * negative impact parameters. The value chosen here is simply the
 * x-component of \f$\vec b\f$. The result will be that the projectile
 * and target will have switched position in x.
 *
 * \key INITIAL_DISTANCE (float, optional, default = no displacement): \n
 * The initial distance of the two nuclei. That
 * means \f$z_{\rm min}^{\rm target} - z_{\rm max}^{\rm projectile}\f$.
 */

    
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
    projectile_ = std::unique_ptr<DeformedNucleus>(new DeformedNucleus());
  } else {
    projectile_ = std::unique_ptr<Nucleus>(new Nucleus());

  }
  if (modus_cfg.has_value({"Target", "DEFORMED"}) &&
      modus_cfg.take({"Target", "DEFORMED"})) {
    target_ = std::unique_ptr<DeformedNucleus>(new DeformedNucleus());
  } else {
    target_ = std::unique_ptr<Nucleus>(new Nucleus());
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
    projectile_->set_parameters_from_config("Projectile", modus_cfg);
  }
  if (modus_cfg.has_value({"Target", "AUTOMATIC"}) &&
      modus_cfg.take({"Target", "AUTOMATIC"})) {
    target_->set_parameters_automatic();
  } else {
    target_->set_parameters_from_config("Target", modus_cfg);
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
      float e_lab = modus_cfg.take({"E_LAB"});
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
      float p_lab = modus_cfg.take({"P_LAB"});
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
      std::string sampling_method = modus_cfg.take({"Impact", "SAMPLE"});
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

std::ostream &operator<<(std::ostream &out, const NucleusModus &m) {
  return out << "-- Nucleus Modus:\n"
                "S (nucleus-nucleus) = " << format(m.total_s_, "GeV")
             << "\nImpact parameter = " << format(m.impact_, "fm")
             << "\nInitial distance between nuclei: "
             << format(2 * m.initial_z_displacement_, "fm")
             << "\nProjectile:\n" << *m.projectile_
             << "\nTarget:\n" << *m.target_;
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
