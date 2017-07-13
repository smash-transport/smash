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
#include "include/fourvector.h"
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
 * \key Sqrtsnn (double, optional, no default): \n
 * Defines the energy of the collision as center-of-mass
 * energy in the collision of one participant each from both nuclei
 * (using the average participant mass in the given nucleus).
 *
 * \key E_Kin (double, optional, no default): \n
 * Defines the energy of the collision by the kinetic energy per nucleon of
 * the projectile nucleus (in AGeV). This assumes the target nucleus is at rest.
 *
 * \key P_Lab (double, optional, no default): \n
 * Defines the energy of the collision by the initial momentum per nucleon
 * of the projectile nucleus (in AGeV). This assumes the target nucleus is at
 * rest.
 *
 * \key Calculation_Frame (string, required, default = "center of velocity"): \n
 * The frame in which the collision is calculated.\n
 * "center of velocity", "center of mass" or "fixed target"
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
 * \li \key Value (double, optional, optional, default = 0.f): fixed value for
 * the impact parameter. No other \key Impact: directive is looked at.
 * \li \key Sample (string, optional, default = \key quadratic): \n
 * if \key uniform, use uniform sampling of the impact parameter
 * (\f$dP(b) = db\f$). If \key quadratic use areal (aka quadratic) input
 * sampling (the probability of an input parameter range is proportional to the
 * area corresponding to that range, \f$dP(b) = b\cdot db\f$). If \key custom,
 * use \key Values and \key Yields to interpolate the impact parameter
 * distribution and use rejection sampling.
 * \li \key Values (doubles, optional):
 * Values of the impact parameter, corresponding to \key Yields. Must be same
 * length as \key Yields. Required for \key Sample = "custom".
 * \li \key Yields (doubles, optional):
 * Values of the particle yields, corresponding to \key Values. Must be same
 * length as \key Values. Required for \key Sample = "custom".
 *
 * \li \key Range (double, double, optional, default = 0.0f):\n
 * A vector of minimal and maximal impact parameters
 * between which b should be chosen. (The order of these is not
 * important.)
 * \li \key Max (double, optional, default = 0.0f):
 * Like `Range: [0.0, Max]`. Note that if both \key Range and
 * \key Max are specified, \key Max takes precedence.
 *
 * Note that there are no safeguards to prevent you from specifying
 * negative impact parameters. The value chosen here is simply the
 * x-component of \f$\vec b\f$. The result will be that the projectile
 * and target will have switched position in x.
 *
 * \key Initial_Distance (double, optional, default = 2.0): \n
 * The initial distance of the two nuclei (in fm). That
 * means \f$z_{\rm min}^{\rm target} - z_{\rm max}^{\rm projectile}\f$.\n
 *
 * Note that this distance is applied before the Lorentz boost
 * to chosen calculation frame, and thus the actual distance may be different.
 *
 * \key Fermi_Motion (string, optional, default = "off"): \n
 * Defines if Fermi motion is included.\n
 * "off", "on" or "frozen"\n
 * Use "frozen" if you want to use Fermi motion
 * without potentials. Use "on" if you want to use Fermi motion
 * in combination with potentials.\n
 *
 */

ColliderModus::ColliderModus(Configuration modus_config,
                             const ExperimentParameters &params) {
  Configuration modus_cfg = modus_config["Collider"];

  // Get the reference frame for the collision calculation.
  if (modus_cfg.has_value({"Calculation_Frame"})) {
    frame_ = modus_cfg.take({"Calculation_Frame"});
  }

  /// Determine whether to avoid the first collsions within the same nucleus
  if (modus_cfg.has_value({"Collisions_Within_Nucleus"})) {
    cll_in_nucleus_ = modus_cfg.take({"Collisions_Within_Nucleus"});
  }

  // Set up the projectile nucleus
  Configuration proj_cfg = modus_cfg["Projectile"];
  if (proj_cfg.has_value({"Deformed"}) && proj_cfg.take({"Deformed"})) {
    projectile_ = make_unique<DeformedNucleus>(proj_cfg, params.testparticles);
  } else {
    projectile_ = make_unique<Nucleus>(proj_cfg, params.testparticles);
  }
  if (projectile_->size() < 1) {
    throw ColliderEmpty("Input Error: Projectile nucleus is empty.");
  }

  // Set up the target nucleus
  Configuration targ_cfg = modus_cfg["Target"];
  if (targ_cfg.has_value({"Deformed"}) && targ_cfg.take({"Deformed"})) {
    target_ = make_unique<DeformedNucleus>(targ_cfg, params.testparticles);
  } else {
    target_ = make_unique<Nucleus>(targ_cfg, params.testparticles);
  }
  if (target_->size() < 1) {
    throw ColliderEmpty("Input Error: Target nucleus is empty.");
  }

  // Get the Fermi-Motion input (off, on, frozen)
  if (modus_cfg.has_value({"Fermi_Motion"})) {
    // We only read the value, because it is still required by the experiment
    // class to make sure we don't use frozen Fermi momenta with potentials.
    fermi_motion_ = modus_cfg.read({"Fermi_Motion"});
  }

  // Get the total nucleus-nucleus collision energy. Since there is
  // no meaningful choice for a default energy, we require the user to
  // give one (and only one) energy input from the available options.
  int energy_input = 0;
  const double mass_projec = projectile_->mass();
  const double mass_target = target_->mass();
  // average mass of a particle in that nucleus
  const double mass_a = projectile_->mass() /
                        projectile_->number_of_particles();
  const double mass_b = target_->mass() / target_->number_of_particles();
  // Option 1: Center of mass energy.
  if (modus_cfg.has_value({"Sqrtsnn"})) {
    sqrt_s_NN_ = modus_cfg.take({"Sqrtsnn"});
    // Check that input satisfies the lower bound (everything at rest).
    if (sqrt_s_NN_ <= mass_a + mass_b) {
      throw ModusDefault::InvalidEnergy(
          "Input Error: sqrt(s_NN) is not larger than masses:\n" +
          std::to_string(sqrt_s_NN_) + " GeV <= " + std::to_string(mass_a) +
          " GeV + " + std::to_string(mass_b) + " GeV.");
    }
    // Set the total nucleus-nucleus collision energy.
    total_s_ = (sqrt_s_NN_ * sqrt_s_NN_ - mass_a * mass_a - mass_b * mass_b) *
                   mass_projec * mass_target / (mass_a * mass_b) +
               mass_projec * mass_projec + mass_target * mass_target;
    energy_input++;
  }
  /* Option 2: Kinetic energy per nucleon of the projectile nucleus
   * (target at rest).  */
  if (modus_cfg.has_value({"E_Kin"})) {
    const double e_kin = modus_cfg.take({"E_Kin"});
    // Check that energy is nonnegative.
    if (e_kin < 0) {
      throw ModusDefault::InvalidEnergy(
          "Input Error: "
          "E_Kin must be nonnegative.");
    }
    // Set the total nucleus-nucleus collision energy.
    total_s_ = s_from_Ekin(e_kin * projectile_->number_of_particles(),
                           mass_projec, mass_target);
    sqrt_s_NN_ = std::sqrt(s_from_Ekin(e_kin, mass_a, mass_b));
    energy_input++;
  }
  // Option 3: Momentum of the projectile nucleus (target at rest).
  if (modus_cfg.has_value({"P_Lab"})) {
    const double p_lab = modus_cfg.take({"P_Lab"});
    // Check that p_lab is nonnegative.
    if (p_lab < 0) {
      throw ModusDefault::InvalidEnergy(
          "Input Error: "
          "P_Lab must be nonnegative.");
    }
    // Set the total nucleus-nucleus collision energy.
    total_s_ = s_from_plab(p_lab * projectile_->number_of_particles(),
                           mass_projec, mass_target);
    sqrt_s_NN_ = std::sqrt(s_from_plab(p_lab, mass_a, mass_b));
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
      if (sampling_ == Sampling::Custom) {
        if (!(modus_cfg.has_value({"Impact", "Values"}) ||
              modus_cfg.has_value({"Impact", "Yields"}))) {
          throw std::domain_error(
              "Input Error: Need impact parameter spectrum for custom "
              "sampling. "
              "Please provide Values and Yields.");
        }
        const std::vector<double> impacts = modus_cfg.take({"Impact", "Values"});
        const std::vector<double> yields = modus_cfg.take({"Impact", "Yields"});
        if (impacts.size() != yields.size()) {
          throw std::domain_error(
              "Input Error: Need as many impact parameter values as yields. "
              "Please make sure that Values and Yields have the same length.");
        }
        impact_interpolation_ = make_unique<InterpolateDataLinear<double>>(
            InterpolateDataLinear<double>(impacts, yields));

        const auto imp_minmax = std::minmax_element(impacts.begin(),
                                                    impacts.end());
        imp_min_ = *imp_minmax.first;
        imp_max_ = *imp_minmax.second;
        yield_max_ = *std::max_element(yields.begin(), yields.end());
      }
    }
    if (modus_cfg.has_value({"Impact", "Range"})) {
      const std::array<double, 2> range = modus_cfg.take({"Impact", "Range"});
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
             << "sqrt(S) (nucleus-nucleus) = "
             << format(std::sqrt(m.total_s_), "GeV\n")
             << "sqrt(S) (nucleon-nucleon) = "
             << format(m.sqrt_s_NN_, "GeV\n")
             << "Initial distance between nuclei: "
             << format(2 * m.initial_z_displacement_, "fm") << "\nProjectile:\n"
             << *m.projectile_ << "\nTarget:\n" << *m.target_;
}

double ColliderModus::initial_conditions(Particles *particles,
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

  // If velocities are larger or equal to 1, throw an exception.
  if (v_a >= 1.0 || v_b >= 1.0) {
    throw std::domain_error(
        "Found velocity equal to or larger than 1 in "
        "ColliderModus::initial_conditions.\nConsider using "
        "the center of velocity reference frame.");
  }

  // Calculate the beam velocity of the projectile and the target, which will be
  // used to calculate the beam momenta in experiment.cc
  if (fermi_motion_ == FermiMotion::Frozen) {
      velocity_projectile_ = v_a;
      velocity_target_ = v_b;
  }

  // Generate Fermi momenta if necessary
  if (fermi_motion_ == FermiMotion::On ||
      fermi_motion_ == FermiMotion::Frozen) {
    // Frozen: Fermi momenta will be ignored during the propagation to
    // avoid that the nuclei will fly apart.
    projectile_->generate_fermi_momenta();
    target_->generate_fermi_momenta();
    if (fermi_motion_ == FermiMotion::On) {
      log.info() << "Fermi motion is ON.";
    } else {
      log.info() << "FROZEN Fermi motion is on.";
    }
  } else if (fermi_motion_ == FermiMotion::Off) {
    // No Fermi-momenta are generated in this case
    log.info() << "Fermi motion is OFF.";
  } else {
    throw std::domain_error("Invalid Fermi_Motion input.");
  }

  // Boost the nuclei to the appropriate velocity.
  projectile_->boost(v_a);
  target_->boost(v_b);

  // Shift the nuclei into starting positions. Contracted spheres with
  // nuclear radii should touch exactly at t=0. Modus starts at negative
  // time corresponding to additinal initial displacement.
  const double d_a = std::max(0.0f, projectile_->get_diffusiveness());
  const double d_b = std::max(0.0f, target_->get_diffusiveness());
  const double r_a = projectile_->get_nuclear_radius();
  const double r_b = target_->get_nuclear_radius();
  const double dz = initial_z_displacement_;

  const double simulation_time = -dz / std::abs(v_a);
  const double proj_z = -dz -
                        std::sqrt(1.0 - v_a*v_a) * (r_a + d_a);
  const double targ_z = +dz * std::abs(v_b/v_a) +
                        std::sqrt(1.0 - v_b*v_b) * (r_b + d_b);
  projectile_->shift(proj_z, +impact_ / 2.0, simulation_time);
  target_->    shift(targ_z, -impact_ / 2.0, simulation_time);

  // Put the particles in the nuclei into code particles.
  projectile_->copy_particles(particles);
  target_->copy_particles(particles);
  return simulation_time;
}

void ColliderModus::sample_impact() {
  switch (sampling_) {
    case Sampling::Quadratic: {
      // quadratic sampling: Note that for bmin > bmax, this still yields
      // the correct distribution (only that canonical() = 0 then is the
      // upper end, not the lower).
      impact_ = std::sqrt(imp_min_ * imp_min_ +
                          Random::canonical() *
                              (imp_max_ * imp_max_ - imp_min_ * imp_min_));
    }
    break;
    case Sampling::Custom: {
      // rejection sampling based on given distribution
      assert(impact_interpolation_ != nullptr);
      double probability_random = 1;
      double probability = 0;
      double b;
      while (probability_random > probability) {
        b = Random::uniform(imp_min_, imp_max_);
        probability = (*impact_interpolation_)(b) / yield_max_;
        assert(probability < 1.0f);
        probability_random = Random::uniform(0.f, 1.f);
      }
      impact_ = b;
    }
    break;
    case Sampling::Uniform: {
      // linear sampling. Still, min > max works fine.
      impact_ = Random::uniform(imp_min_, imp_max_);
    }
  }
}

std::pair<double, double> ColliderModus::get_velocities(double s, double m_a,
                                                        double m_b) {
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
        double pCM = pCM_from_s(s, m_a, m_b);
        v_a = pCM / std::sqrt(m_a*m_a + pCM*pCM);
        v_b = -pCM / std::sqrt(m_b*m_b + pCM*pCM);
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
