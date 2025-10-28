/*
 *    Copyright (c) 2014-2022,2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#include "smash/deformednucleus.h"

#include <cmath>
#include <map>
#include <stdexcept>

#include "smash/configuration.h"
#include "smash/constants.h"
#include "smash/fourvector.h"
#include "smash/input_keys.h"
#include "smash/random.h"
#include "smash/threevector.h"

namespace smash {

DeformedNucleus::DeformedNucleus(const std::map<PdgCode, int> &particle_list,
                                 int nTest,
                                 SpinInteractionType spin_interaction_type)
    : Nucleus(particle_list, nTest, spin_interaction_type) {}

DeformedNucleus::DeformedNucleus(Configuration &config, int nTest,
                                 bool auto_deformation)
    : Nucleus(config, nTest) {
  if (auto_deformation) {
    set_deformation_parameters_automatic();
    set_saturation_density(calculate_saturation_density());
  } else {
    set_deformation_parameters_from_config(config);
  }
  /* If the config does not contain (anymore) a target or a projectile
  sub-section, this code should not be executed, because e.g. the
  is_about_projectile function would fail. */
  if (has_projectile_or_target(config)) {
    const auto &orientation_section = [&config]() {
      return is_about_projectile(config) ? InputSections::m_c_p_orientation
                                         : InputSections::m_c_t_orientation;
    }();
    if (config.has_section(orientation_section)) {
      Configuration sub_conf =
          config.extract_complete_sub_configuration(orientation_section);
      set_orientation_from_config(sub_conf);
    }
  }
}

ThreeVector DeformedNucleus::distribute_nucleon() {
  double a_radius;
  Angles a_direction;
  // Set a sensible maximum bound for radial sampling.
  double radius_max =
      Nucleus::get_nuclear_radius() / Nucleus::get_diffusiveness() +
      Nucleus::get_nuclear_radius() * Nucleus::get_diffusiveness();

  // Sample the distribution.
  do {
    a_direction.distribute_isotropically();
    // sample r**2 dr
    a_radius = radius_max * std::cbrt(random::canonical());
  } while (random::canonical() > nucleon_density(a_radius,
                                                 a_direction.costheta(),
                                                 a_direction.phi()) /
                                     Nucleus::get_saturation_density());

  // Update (x, y, z) positions.
  return a_direction.threevec() * a_radius;
}

void DeformedNucleus::set_deformation_parameters_automatic() {
  // Set the deformation parameters
  // reference for U, Pb, Au, Cu: \iref{Moller:1993ed}
  // reference for Zr and Ru: \iref{Schenke:2019ruo}
  // reference for Xe: \iref{Moller:2015fba}
  bool listed = 0;
  const std::map<int, std::string> A_map = {{238, "Uranium"},
                                            {208, "Lead"},
                                            {197, "Gold"},
                                            {63, "Copper"},
                                            {129, "Xenon"}};
  const std::map<std::string, std::string> Z_map = {{"Uranium", "92"},
                                                    {"Lead", "82"},
                                                    {"Gold", "79"},
                                                    {"Copper", "29"},
                                                    {"Xenon", "54"}};
  int A = Nucleus::number_of_particles();
  int Z = Nucleus::number_of_protons();
  switch (A) {
    case 238:  // Uranium
      if (Z == 92) {
        set_beta_2(0.28);
        set_beta_4(0.093);
      } else {
        listed = true;
      }
      break;
    case 208:  // Lead
      if (Z == 82) {
        set_beta_2(0.0);
        set_beta_4(0.0);
      } else {
        listed = true;
      }
      break;
    case 197:  // Gold
      if (Z == 79) {
        set_beta_2(-0.131);
        set_beta_4(-0.031);
      } else {
        listed = true;
      }
      break;
    case 129:  // Xenon
      if (Z == 54) {
        set_beta_2(0.162);
        set_beta_4(-0.003);
      } else {
        listed = true;
      }
      break;
    case 63:  // Copper
      if (Z == 29) {
        set_beta_2(0.162);
        set_beta_4(-0.006);
      } else {
        listed = true;
      }
      break;
    case 96:
      if (Z == 40) {  // Zirconium
        set_beta_2(0.0);
        set_beta_4(0.0);
      } else if (Z == 44) {  // Ruthenium
        set_beta_2(0.158);
        set_beta_4(0.0);
      } else {
        throw std::domain_error(
            "Number of protons for nuclei with mass number A = 96 does not "
            "match that of Zirconium or Ruthenium. The deformation parameters "
            "for additional isobars are currently not implemented."
            " Please specify at least \"Beta_2\" and \"Beta_4\" "
            "manually and set \"Automatic: False.\" ");
      }
      break;
    default:
      throw std::domain_error(
          "Mass number not listed for automatically setting deformation "
          "parameters. Please specify at least \"Beta_2\" and \"Beta_4\" "
          "manually and set \"Automatic: False.\" ");
  }
  if (listed) {
    throw std::domain_error("Mass number is listed under " + A_map.at(A) +
                            " but the proton "
                            "number of " +
                            std::to_string(Z) +
                            " does not match "
                            "its " +
                            Z_map.at(A_map.at(A)) +
                            " protons."
                            "Please specify at least \"Beta_2\" and \"Beta_4\" "
                            "manually and set \"Automatic: False.\" ");
  }
}

void DeformedNucleus::set_deformation_parameters_from_config(
    Configuration &config) {
  if (has_projectile_or_target(config)) {
    const bool is_projectile = is_about_projectile(config);
    const auto &[beta2_key, beta3_key, beta4_key,
                 gamma_key] = [&is_projectile]() {
      return is_projectile
                 ? std::make_tuple(
                       InputKeys::modi_collider_projectile_deformed_beta2,
                       InputKeys::modi_collider_projectile_deformed_beta3,
                       InputKeys::modi_collider_projectile_deformed_beta4,
                       InputKeys::modi_collider_projectile_deformed_gamma)
                 : std::make_tuple(
                       InputKeys::modi_collider_target_deformed_beta2,
                       InputKeys::modi_collider_target_deformed_beta3,
                       InputKeys::modi_collider_target_deformed_beta4,
                       InputKeys::modi_collider_target_deformed_gamma);
    }();
    // Deformation parameters
    if (config.has_value(beta2_key)) {
      beta2_ = config.take(beta2_key);
    }
    if (config.has_value(gamma_key)) {
      gamma_ = config.take(gamma_key);
    }
    if (config.has_value(beta3_key)) {
      beta3_ = config.take(beta3_key);
    }
    if (config.has_value(beta4_key)) {
      beta4_ = config.take(beta4_key);
    }
  }
}

double y_l_m(int l, int m, double cosx, double phi) {
  if (l == 2 && m == 0) {
    return (1. / 4) * std::sqrt(5 / M_PI) * (3. * (cosx * cosx) - 1);
  } else if (l == 2 && m == 2) {
    double sinx2 = 1. - cosx * cosx;
    return (1. / 4) * std::sqrt(15 / (2. * M_PI)) * sinx2 * std::cos(2. * phi);
  } else if (l == 3 && m == 0) {
    return (1. / 4) * std::sqrt(7 / M_PI) *
           (5. * cosx * (cosx * cosx) - 3. * cosx);
  } else if (l == 4 && m == 0) {
    return (3. / 16) * std::sqrt(1 / M_PI) *
           (35. * (cosx * cosx) * (cosx * cosx) - 30. * (cosx * cosx) + 3);
  } else {
    throw std::domain_error(
        "Not a valid angular momentum quantum number in y_l_m.");
  }
}

double DeformedNucleus::nucleon_density(double r, double cosx,
                                        double phi) const {
  return Nucleus::get_saturation_density() /
         (1 + std::exp((r - Nucleus::get_nuclear_radius() *
                                (1 +
                                 beta2_ * (std::cos(gamma_) *
                                               y_l_m(2, 0, cosx, phi) +
                                           std::sqrt(2) * std::sin(gamma_) *
                                               y_l_m(2, 2, cosx, phi)) +
                                 beta3_ * y_l_m(3, 0, cosx, phi) +
                                 beta4_ * y_l_m(4, 0, cosx, phi))) /
                       Nucleus::get_diffusiveness()));
}

double DeformedNucleus::nucleon_density_unnormalized(double r, double cosx,
                                                     double phi) const {
  return 1.0 /
         (1 + std::exp((r - Nucleus::get_nuclear_radius() *
                                (1 +
                                 beta2_ * (std::cos(gamma_) *
                                               y_l_m(2, 0, cosx, phi) +
                                           std::sqrt(2) * std::sin(gamma_) *
                                               y_l_m(2, 2, cosx, phi)) +
                                 beta3_ * y_l_m(3, 0, cosx, phi) +
                                 beta4_ * y_l_m(4, 0, cosx, phi))) /
                       Nucleus::get_diffusiveness()));
}

double DeformedNucleus::integrant_nucleon_density_phi(double r,
                                                      double cosx) const {
  Integrator integrate;
  // Perform the phi integration. This is needed if the triaxiality coefficient
  // gamma is included, which includes a dependency around the phi axis.
  // Unfortunately the Integrator class does not support 3d integration which is
  // why this intermediate integral is needed. It has been checked that the
  // integral factorizes.
  const auto result = integrate(0.0, 2.0 * M_PI, [&](double phi) {
    return nucleon_density_unnormalized(r, cosx, phi);
  });
  return result.value();
}

double DeformedNucleus::calculate_saturation_density() const {
  Integrator2d integrate;
  // Transform integral from (0, oo) to (0, 1) via r = (1 - t) / t.
  // To prevent overflow, the integration is only performed to t = 0.01 which
  // corresponds to r = 99fm. Additionally the precision settings in the
  // Integrator2d scheme are equally important. However both these point affect
  // the result only after the seventh digit which should not be relevant here.
  if (gamma_ == 0.0) {
    const auto result = integrate(0.01, 1, -1, 1, [&](double t, double cosx) {
      const double r = (1 - t) / t;
      return twopi * std::pow(r, 2.0) *
             nucleon_density_unnormalized(r, cosx, 0.0) / std::pow(t, 2.0);
    });
    const auto rho0 = number_of_particles() / result.value();
    return rho0;
  } else {
    const auto result = integrate(0.01, 1, -1, 1, [&](double t, double cosx) {
      const double r = (1 - t) / t;
      return std::pow(r, 2.0) * integrant_nucleon_density_phi(r, cosx) /
             std::pow(t, 2.0);
    });
    const auto rho0 = number_of_particles() / result.value();
    return rho0;
  }
}

}  // namespace smash
