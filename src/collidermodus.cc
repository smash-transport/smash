/*
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "include/collidermodus.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "include/angles.h"
#include "include/configuration.h"
#include "include/cxx14compat.h"
#include "include/experimentparameters.h"
#include "include/interpolation.h"
#include "include/kinematics.h"
#include "include/logging.h"
#include "include/numerics.h"
#include "include/particles.h"
#include "include/pdgcode.h"
#include "include/random.h"

namespace Smash {

/*!\Userguide
 * \page input_modi_collider_ Collider
 *
 * Possible Incident Energies, only one can be given:
 *
 * \key Sqrtsnn (float, optional, no default): \n
 * Defines the energy of the collision as center-of-mass
 * energy in the collision of one participant each from both nuclei.
 * Optional: Since not all participants have the same mass, and hence
 * \f$\sqrt{s_{\rm NN}}\f$ is different for \f$NN\f$ = proton+proton and
 * \f$NN\f$=neutron+neutron, you can specify which \f$NN\f$-pair you
 * want this to refer to with `Sqrts_Reps`. This expects a vector of two
 * PDG Codes, e.g. `Sqrts_n: [2212, 2212]` for proton-proton.
 * Default is the average nucleon mass in the given nucleus.
 *
 * \key E_Kin (float, optional, no default): \n
 * Defines the energy of the collision by the kinetic energy per nucleon of
 * the projectile nucleus (in AGeV). This assumes the target nucleus is at rest.
 *
 * \key P_Lab (float, optional, no default): \n
 * Defines the energy of the collision by the initial momentum per nucleon
 * of the projectile nucleus (in AGeV). This assumes the target nucleus is at
 * rest.
 *
 * \key Calculation_Frame (int, required, default = 1): \n
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
 * \li \key Particles (int:int, int:int, required):\n
 * A map in which the keys are PDG codes and the
 * values are number of particles with that PDG code that should be in
 * the current nucleus. E.g.\ `Particles: {2212: 82, 2112: 126}` for a
 * lead-208 nucleus (82 protons and 126 neutrons = 208 nucleons), and
 * `Particles: {2212: 1, 2112: 1, 3122: 1}` for Hyper-Triton (one
 * proton, one neutron and one Lambda).
 * \li \key Deformed (bool, optional, default = false): \n
 * true - deformed nucleus is initialized
 * false - spherical nucleus is initialized
 * \li \key Automatic (bool, optional, default = true): \n
 * true - sets all necessary parameters based on the atomic number
 * of the input nucleus \n
 * false - manual values according to deformed nucleus (see below)
 *
 * Additional Woods-Saxon Parameters: \n
 * There are also many other
 * parameters for specifying the shape of the Woods-Saxon distribution,
 * and other nucleus specific properties. See nucleus.cc and
 * deformednucleus.cc for more on these choices.
 *
 * \li \subpage input_deformed_nucleus_
 *
 * \key Impact: \n
 * A section for the impact parameter (= distance (in fm) of the two
 * straight lines that the center of masses of the nuclei travel on).
 *
 * \li \key Value (float, optional, optional, default = 0.f): fixed value for
 * the impact parameter. No other \key Impact: directive is looked at.
 * \li \key Sample (string, optional, default = \key quadratic): \n
 * if \key uniform, use uniform sampling of the impact parameter
 * (\f$dP(b) = db\f$). If \key quadratic use areal (aka quadratic) input
 * sampling (the probability of an input parameter range is proportional to the
 * area corresponding to that range, \f$dP(b) = b\cdot db\f$). If \key custom,
 * use \key Values and \key Yields to interpolate the impact parameter
 * distribution and use rejection sampling.
 * \li \key Values (floats, optional):
 * Values of the impact parameter, corresponding to \key Yields. Must be same
 * length as \key Yields. Required for \key Sample = "custom".
 * \li \key Yields (floats, optional):
 * Values of the particle yields, corresponding to \key Values. Must be same
 * lenght as \key Values. Required for \key Sample = "custom".
 *
 * \li \key Range (float, float, optional, default = 0.0f):\n
 * A vector of minimal and maximal impact parameters
 * between which b should be chosen. (The order of these is not
 * important.)
 * \li \key Max (float, optional, default = 0.0f):
 * Like `Range: [0.0, Max]`. Note that if both \key Range and
 * \key Max are specified, \key Max takes precedence.
 *
 * Note that there are no safeguards to prevent you from specifying
 * negative impact parameters. The value chosen here is simply the
 * x-component of \f$\vec b\f$. The result will be that the projectile
 * and target will have switched position in x.
 *
 * \key Initial_Distance (float, optional, default = 2.0): \n
 * The initial distance of the two nuclei (in fm). That
 * means \f$z_{\rm min}^{\rm target} - z_{\rm max}^{\rm projectile}\f$.\n
 *
 * Note that this distance is applied before the Lorentz boost
 * to chosen calculation frame, and thus the actual distance may be different.
 *
 * \key Fermi_Motion (bool, optional, default = false): \n
 * Defines if Fermi motion is included. Note that Fermi motion
 * is senseless physicswise if potentials are off: without potentials
 * nucleons will just fly apart.
 */

ColliderModus::ColliderModus(Configuration modus_config,
                             const ExperimentParameters &params) {
  Configuration modus_cfg = modus_config["Collider"];

  // Get the reference frame for the collision calculation.
  if (modus_cfg.has_value({"Calculation_Frame"})) {
    frame_ = modus_cfg.take({"Calculation_Frame"});
  }

  // Decide which type of nucleus: deformed or not (default).
  if (modus_cfg.has_value({"Projectile", "Deformed"}) &&
      modus_cfg.take({"Projectile", "Deformed"})) {
    projectile_ = std::unique_ptr<DeformedNucleus>(new DeformedNucleus());
  } else {
    projectile_ = std::unique_ptr<Nucleus>(new Nucleus());
  }
  if (modus_cfg.has_value({"Target", "Deformed"}) &&
      modus_cfg.take({"Target", "Deformed"})) {
    target_ = std::unique_ptr<DeformedNucleus>(new DeformedNucleus());
  } else {
    target_ = std::unique_ptr<Nucleus>(new Nucleus());
  }

  // Fill nuclei with particles.
  std::map<PdgCode, int> pro = modus_cfg.take({"Projectile", "Particles"});
  projectile_->fill_from_list(pro, params.testparticles);
  if (projectile_->size() < 1) {
    throw ColliderEmpty("Input Error: Projectile nucleus is empty.");
  }
  std::map<PdgCode, int> tar = modus_cfg.take({"Target", "Particles"});
  target_->fill_from_list(tar, params.testparticles);
  if (target_->size() < 1) {
    throw ColliderEmpty("Input Error: Target nucleus is empty.");
  }

  // Ask to construct nuclei based on atomic number; otherwise, look
  // for the user defined values or take the default parameters.
  if (modus_cfg.has_value({"Projectile", "Automatic"}) &&
      modus_cfg.take({"Projectile", "Automatic"})) {
    projectile_->set_parameters_automatic();
  } else {
    projectile_->set_parameters_from_config("Projectile", modus_cfg);
  }
  if (modus_cfg.has_value({"Target", "Automatic"}) &&
      modus_cfg.take({"Target", "Automatic"})) {
    target_->set_parameters_automatic();
  } else {
    target_->set_parameters_from_config("Target", modus_cfg);
  }

  // Consider an option to include Fermi-motion
  fermi_motion_ = modus_cfg.take({"Fermi_Motion"}, false);

  // Get the total nucleus-nucleus collision energy. Since there is
  // no meaningful choice for a default energy, we require the user to
  // give one (and only one) energy input from the available options.
  int energy_input = 0;
  float mass_projec = projectile_->mass();
  float mass_target = target_->mass();
  // Option 1: Center of mass energy.
  if (modus_cfg.has_value({"Sqrtsnn"})) {
    float sqrt_s_NN = modus_cfg.take({"Sqrtsnn"});
    /* Note that \f$\sqrt{s_{NN}}\f$ is different for neutron-neutron and
     * proton-proton collisions (because of the different masses). Therefore,
     * representative particles are needed to specify which two particles'
     * collisions have this \f$\sqrt{s_{NN}}\f$. The vector specifies a pair
     * of PDG codes for the two particle species we want to use. The default
     * is otherwise the average nucleon mass for each nucleus.  */
    PdgCode id_a = 0, id_b = 0;
    if (modus_cfg.has_value({"Sqrts_Reps"})) {
      std::vector<PdgCode> sqrts_reps = modus_cfg.take({"Sqrts_Reps"});
      id_a = sqrts_reps[0];
      id_b = sqrts_reps[1];
    }
    float mass_a, mass_b;
    if (id_a != 0) {
      // If PDG Code is given, use mass of this particle type.
      mass_a = ParticleType::find(id_a).mass();
    } else {
      // else, use average mass of a particle in that nucleus
      mass_a = projectile_->mass() / projectile_->size();
    }
    if (id_b != 0) {
      mass_b = ParticleType::find(id_b).mass();
    } else {
      mass_b = target_->mass() / target_->size();
    }
    // Check that input satisfies the lower bound (everything at rest).
    if (sqrt_s_NN < mass_a + mass_b) {
      throw ModusDefault::InvalidEnergy(
          "Input Error: sqrt(s_NN) is smaller than masses:\n" +
          std::to_string(sqrt_s_NN) + " GeV < " + std::to_string(mass_a) +
          " GeV + " + std::to_string(mass_b) + " GeV.");
    }
    // Set the total nucleus-nucleus collision energy.
    total_s_ = (sqrt_s_NN * sqrt_s_NN - mass_a * mass_a - mass_b * mass_b) *
                   mass_projec * mass_target / (mass_a * mass_b) +
               mass_projec * mass_projec + mass_target * mass_target;
    energy_input++;
  }
  /* Option 2: Kinetic energy per nucleon of the projectile nucleus
   * (target at rest).  */
  if (modus_cfg.has_value({"E_Kin"})) {
    float e_kin = modus_cfg.take({"E_Kin"});
    // Check that energy is nonnegative.
    if (e_kin < 0) {
      throw ModusDefault::InvalidEnergy(
          "Input Error: "
          "E_Kin must be nonnegative.");
    }
    // Set the total nucleus-nucleus collision energy.
    total_s_ = s_from_Ekin(e_kin * projectile_->number_of_particles(),
                           mass_projec, mass_target);
    energy_input++;
  }
  // Option 3: Momentum of the projectile nucleus (target at rest).
  if (modus_cfg.has_value({"P_Lab"})) {
    float p_lab = modus_cfg.take({"P_Lab"});
    // Check that p_lab is nonnegative.
    if (p_lab < 0) {
      throw ModusDefault::InvalidEnergy(
          "Input Error: "
          "P_Lab must be nonnegative.");
    }
    // Set the total nucleus-nucleus collision energy.
    total_s_ = s_from_plab(p_lab * projectile_->number_of_particles(),
                           mass_projec, mass_target);
    energy_input++;
  }
  if (energy_input == 0) {
    throw std::domain_error(
        "Input Error: Non-existent collision energy. "
        "Please provide one of Sqrtsnn/E_Kin/P_Lab.");
  }
  if (energy_input > 1) {
    throw std::domain_error(
        "Input Error: Redundant collision energy. "
        "Please provide only one of Sqrtsnn/E_Kin/P_Lab.");
  }

  // Impact parameter setting: Either "Value", "Range", "Max" or "Sample".
  // Unspecified means 0 impact parameter.
  if (modus_cfg.has_value({"Impact", "Value"})) {
    impact_ = modus_cfg.take({"Impact", "Value"});
    imp_min_ = impact_;
    imp_max_ = impact_;
  } else {
    // If impact is not supplied by value, inspect sampling parameters:
    if (modus_cfg.has_value({"Impact", "Sample"})) {
      sampling_ = modus_cfg.take({"Impact", "Sample"});
      if (sampling_ == Sampling::CUSTOM) {
        if (!(modus_cfg.has_value({"Impact", "Values"}) ||
              modus_cfg.has_value({"Impact", "Yields"}))) {
          throw std::domain_error(
              "Input Error: Need impact parameter spectrum for custom "
              "sampling. "
              "Please provide Values and Yields.");
        }
        const std::vector<float> impacts = modus_cfg.take({"Impact", "Values"});
        const std::vector<float> yields = modus_cfg.take({"Impact", "Yields"});
        if (impacts.size() != yields.size()) {
          throw std::domain_error(
              "Input Error: Need as many impact parameter values as yields. "
              "Please make sure that Values and Yields have the same length.");
        }
        impact_interpolation_ = make_unique<InterpolateData<float>>(
            InterpolateData<float>(impacts, yields));

        auto imp_minmax = std::minmax_element(impacts.begin(), impacts.end());
        imp_min_ = *imp_minmax.first;
        imp_max_ = *imp_minmax.second;
        yield_max_ = *std::max_element(yields.begin(), yields.end());
      }
    }
    if (modus_cfg.has_value({"Impact", "Range"})) {
      std::array<float, 2> range = modus_cfg.take({"Impact", "Range"});
      imp_min_ = range[0];
      imp_max_ = range[1];
    }
    if (modus_cfg.has_value({"Impact", "Max"})) {
      imp_min_ = 0.0;
      imp_max_ = modus_cfg.take({"Impact", "Max"});
    }
  }

  // Look for user-defined initial separation between nuclei.
  if (modus_cfg.has_value({"Initial_Distance"})) {
    initial_z_displacement_ = modus_cfg.take({"Initial_Distance"});
    // the displacement is half the distance (both nuclei are shifted
    // initial_z_displacement_ away from origin)
    initial_z_displacement_ /= 2.0;
  }
}

std::ostream &operator<<(std::ostream &out, const ColliderModus &m) {
  return out << "-- Collider Modus:\n"
                "sqrt(S) (nucleus-nucleus) = "
             << format(std::sqrt(m.total_s_), "GeV")
             << "\nInitial distance between nuclei: "
             << format(2 * m.initial_z_displacement_, "fm") << "\nProjectile:\n"
             << *m.projectile_ << "\nTarget:\n" << *m.target_;
}

float ColliderModus::initial_conditions(Particles *particles,
                                        const ExperimentParameters &) {
  const auto &log = logger<LogArea::Collider>();
  // Sample impact parameter distribution.
  sample_impact();

  log.info() << "Impact parameter = " << format(impact_, "fm");
  // Populate the nuclei with appropriately distributed nucleons.
  // If deformed, this includes rotating the nucleus.
  projectile_->arrange_nucleons();
  target_->arrange_nucleons();

  // Use the total mandelstam variable to get the frame-dependent velocity for
  // each nucleus. Position a is projectile, position b is target.
  double v_a, v_b;
  std::tie(v_a, v_b) =
      get_velocities(total_s_, projectile_->mass(), target_->mass());

  // If velocities are too close to 1 for our calculations, throw an exception.
  if (almost_equal(std::abs(1.0 - v_a), 0.0) ||
      almost_equal(std::abs(1.0 - v_b), 0.0)) {
    throw std::domain_error(
        "Found velocity equal to 1 in "
        "nucleusmodus::initial_conditions.\nConsider using"
        "the center of velocity reference frame.");
  }

  // Shift the nuclei into starting positions. Keep the pair separated
  // in z by some small distance and shift in x by the impact parameter.
  // (Projectile is chosen to hit at positive x.)
  // After shifting, set the time component of the particles to
  // -initial_z_displacement_/average_velocity.
  float avg_velocity = std::sqrt(v_a * v_a + v_b * v_b);
  float simulation_time = -initial_z_displacement_ / avg_velocity;
  projectile_->shift(true, -initial_z_displacement_, +impact_ / 2.0,
                     simulation_time);
  target_->shift(false, initial_z_displacement_, -impact_ / 2.0,
                 simulation_time);

  // Generate Fermi momenta if necessary
  if (fermi_motion_) {
    log.info() << "Fermi motion is ON";
    projectile_->generate_fermi_momenta();
    target_->generate_fermi_momenta();
  }

  // Boost the nuclei to the appropriate velocity.
  projectile_->boost(v_a);
  target_->boost(v_b);

  // Put the particles in the nuclei into code particles.
  projectile_->copy_particles(particles);
  target_->copy_particles(particles);
  return simulation_time;
}

void ColliderModus::sample_impact() {
  if (sampling_ == Sampling::QUADRATIC) {
    // quadratic sampling: Note that for bmin > bmax, this still yields
    // the correct distribution (only that canonical() = 0 then is the
    // upper end, not the lower).
    impact_ = std::sqrt(imp_min_ * imp_min_ +
                        Random::canonical() *
                            (imp_max_ * imp_max_ - imp_min_ * imp_min_));
  } else if (sampling_ == Sampling::CUSTOM) {
    // rejection sampling based on given distribution
    assert(impact_interpolation_ != nullptr);
    float probability_random = 1;
    float probability = 0;
    float b;
    while (probability_random > probability) {
      b = Random::uniform(imp_min_, imp_max_);
      probability = (*impact_interpolation_)(b) / yield_max_;
      assert(probability < 1.0f);
      probability_random = Random::uniform(0.f, 1.f);
    }
    impact_ = b;
  } else {
    // linear sampling. Still, min > max works fine.
    impact_ = Random::uniform(imp_min_, imp_max_);
  }
}

std::pair<double, double> ColliderModus::get_velocities(float s, float m_a,
                                                        float m_b) {
  double v_a = 0.0;
  double v_b = 0.0;
  // Frame dependent calculations of velocities. Assume v_a >= 0, v_b <= 0.
  switch (frame_) {
    case CalculationFrame::CenterOfVelocity:
      v_a = center_of_velocity_v(s, m_a, m_b);
      v_b = -v_a;
      break;
    case CalculationFrame::CenterOfMass:
      {
        // Compute center of mass momentum.
        double abs_p = pCM_from_s(s, m_a, m_b);
        v_a = abs_p / m_a;
        v_b = -abs_p / m_b;
      }
      break;
    case CalculationFrame::FixedTarget:
      v_a = fixed_target_projectile_v(s, m_a, m_b);
      break;
    default:
      throw std::domain_error(
          "Invalid reference frame in "
          "ColliderModus::get_velocities.");
  }
  return std::make_pair(v_a, v_b);
}

}  // namespace Smash
