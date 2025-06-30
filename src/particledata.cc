/*
 *
 *    Copyright (c) 2014-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/particledata.h"

#include <iomanip>
#include <iostream>
#include <optional>
#include <vector>

#include "smash/constants.h"
#include "smash/iomanipulators.h"
#include "smash/logging.h"
#include "smash/numerics.h"

namespace smash {

double ParticleData::effective_mass() const {
  const double m_pole = pole_mass();
  if (m_pole < really_small) {
    // prevent numerical problems with massless or very light particles
    return m_pole;
  } else {
    return momentum().abs();
  }
}

void ParticleData::set_history(int ncoll, uint32_t pid, ProcessType pt,
                               double time_last_coll,
                               const ParticleList &plist) {
  if ((pt != ProcessType::Wall) && (pt != ProcessType::FluidizationNoRemoval)) {
    history_.collisions_per_particle = ncoll;
    history_.time_last_collision = time_last_coll;
  }
  history_.id_process = pid;
  history_.process_type = pt;
  switch (pt) {
    case ProcessType::Decay:
    case ProcessType::Wall:
      // only store one parent
      history_.p1 = plist[0].pdgcode();
      history_.p2 = 0x0;
      break;
    case ProcessType::Elastic:
    case ProcessType::FailedString:
    case ProcessType::Fluidization:
    case ProcessType::FluidizationNoRemoval:
      // Parent particles are not updated by the elastic scatterings,
      // failed string processes, or fluidizations
      break;
    case ProcessType::TwoToOne:
    case ProcessType::TwoToTwo:
    case ProcessType::TwoToThree:
    case ProcessType::TwoToFour:
    case ProcessType::TwoToFive:
    case ProcessType::StringSoftSingleDiffractiveAX:
    case ProcessType::StringSoftSingleDiffractiveXB:
    case ProcessType::StringSoftDoubleDiffractive:
    case ProcessType::StringSoftAnnihilation:
    case ProcessType::StringSoftNonDiffractive:
    case ProcessType::StringHard:
    case ProcessType::Bremsstrahlung:
      // store two parent particles
      history_.p1 = plist[0].pdgcode();
      history_.p2 = plist[1].pdgcode();
      break;
    case ProcessType::Thermalization:
    case ProcessType::MultiParticleThreeMesonsToOne:
    case ProcessType::MultiParticleThreeToTwo:
    case ProcessType::MultiParticleFourToTwo:
    case ProcessType::MultiParticleFiveToTwo:
    case ProcessType::Freeforall:
    case ProcessType::None:
      // nullify parents
      history_.p1 = 0x0;
      history_.p2 = 0x0;
      break;
  }
}

void ParticleData::set_unpolarized_spin_vector() {
  // Check whether the velocity of a particle is set and not nan
  if (std::isnan(velocity().x1()) || std::isnan(velocity().x2()) ||
      std::isnan(velocity().x3())) {
    std::cout << "Velocity not set!" << std::endl;
  }

  if (!std::isnan(spin_vector_.x0()) || !std::isnan(spin_vector_.x1()) ||
      !std::isnan(spin_vector_.x2()) || !std::isnan(spin_vector_.x3())) {
    std::cout << "Spin vector already set for particle index: " << index_ << std::endl;
  }

  // Set the mean and standard deviation for the normal distribution
  // of the spin vector components. This random sampling is solely for
  // initializing particle spin vectors and is not physical. Hence, the
  // arbitrary standard deviation is set to a small value for testing.
  double mean = 0.0;
  double sigma = 0.75;

  const int total_spin = spin();
  // For spin 0 all 4 components of the Pauli-Lubanski vector are zero
  if (total_spin == 0) {
    spin_vector_ = FourVector(0., 0., 0., 0.);
    /**
     * For finite spin states, we need to choose the spin vector components such
     * that the average polarization is zero. This is achieved by using \f$S^0 =
     * 0$ in the particle rest frame. The spatial components of the spin vector
     * are sampled from a normal distribution with mean 0. After sampling, we
     * boost the spin vector to lab frame.
     */
  } else {
    FourVector spin_vector(0., random::normal(mean, sigma),
                           random::normal(mean, sigma),
                           random::normal(mean, sigma));
    // Boost the spin vector to the lab frame
    spin_vector_ = spin_vector.lorentz_boost(velocity());
  }
}

double ParticleData::xsec_scaling_factor(double delta_time) const {
  // if formation times are NaNs simply return a NaN to avoid floating-point
  // FE_INVALID errors (that would be raised by unordered comparisons with NaNs)
  if (std::isnan(formation_time_) || std::isnan(begin_formation_time_)) {
    return smash_NaN<double>;
  }

  double time_of_interest = position_.x0() + delta_time;
  // cross section scaling factor at the time_of_interest
  double scaling_factor;

  if (formation_power_ <= 0.) {
    // use a step function to form particles
    if (time_of_interest < formation_time_) {
      // particles will not be fully formed at time of interest
      scaling_factor = initial_xsec_scaling_factor_;
    } else {
      // particles are fully formed at time of interest
      scaling_factor = 1.;
    }
  } else {
    // use smooth function to scale cross section (unless particles are already
    // fully formed at desired time or will start to form later)
    if (formation_time_ <= time_of_interest) {
      // particles are fully formed when colliding
      scaling_factor = 1.;
    } else if (begin_formation_time_ >= time_of_interest) {
      // particles will start formimg later
      scaling_factor = initial_xsec_scaling_factor_;
    } else {
      // particles are in the process of formation at the given time
      scaling_factor =
          initial_xsec_scaling_factor_ +
          (1. - initial_xsec_scaling_factor_) *
              std::pow((time_of_interest - begin_formation_time_) /
                           (formation_time_ - begin_formation_time_),
                       formation_power_);
    }
  }
  return scaling_factor;
}

std::ostream &operator<<(std::ostream &out, const ParticleData &p) {
  out.fill(' ');
  return out << p.type().name() << " (" << std::setw(5) << p.type().pdgcode()
             << ")" << std::right << "{id:" << field<6> << p.id()
             << ", process:" << field<4> << p.id_process()
             << ", pos [fm]:" << p.position() << ", mom [GeV]:" << p.momentum()
             << ", formation time [fm]:" << p.formation_time()
             << ", cross section scaling factor:" << p.xsec_scaling_factor()
             << "}";
}

std::ostream &operator<<(std::ostream &out, const ParticleList &particle_list) {
  auto column = out.tellp();
  out << '[';
  for (const auto &p : particle_list) {
    if (out.tellp() - column >= 201) {
      out << '\n';
      column = out.tellp();
      out << ' ';
    }
    out << std::setw(5) << std::setprecision(3) << p.momentum().abs3()
        << p.type().name();
  }
  return out << ']';
}

std::ostream &operator<<(std::ostream &out,
                         const PrintParticleListDetailed &particle_list) {
  bool first = true;
  out << '[';
  for (const auto &p : particle_list.list) {
    if (first) {
      first = false;
    } else {
      out << "\n ";
    }
    out << p;
  }
  return out << ']';
}

double ParticleData::formation_power_ = 0.0;

ParticleData create_valid_smash_particle_matching_provided_quantities(
    PdgCode pdgcode, double mass, const FourVector &four_position,
    const FourVector &four_momentum, int log_area, bool &mass_warning,
    bool &on_shell_warning) {
  // Check input position and momentum for nan values
  if (is_any_nan(four_position) || is_any_nan(four_momentum)) {
    logg[log_area].fatal() << "Input particle has at least one nan value in "
                              "position and/or momentum four vector.";
    throw std::invalid_argument(
        "Invalid input (nan) for particle position or momentum.");
  }

  // Some preliminary tool to avoid duplication later
  static const auto emph = einhard::Yellow_t_::ANSI();
  static const auto restore_default = einhard::NoColor_t_::ANSI();
  auto prepare_needed_warnings = [&mass_warning, &on_shell_warning, &mass,
                                  &four_momentum](const ParticleData &p) {
    std::array<std::optional<std::string>, 2> warnings{};
    if (mass_warning) {
      warnings[0] = "Provided mass of stable particle " + p.type().name() +
                    " = " + std::to_string(mass) +
                    " [GeV] is inconsistent with value = " +
                    std::to_string(p.pole_mass()) + " [GeV] from " +
                    "particles file.\nForcing E = sqrt(p^2 + m^2)" +
                    ", where m is the mass contained in the particles file." +
                    "\nFurther warnings about discrepancies between the " +
                    "input mass and the mass contained in the particles file" +
                    " will be suppressed.\n" + emph + "Please make sure" +
                    " that changing input particle properties is an " +
                    "acceptable behavior." + restore_default;
    }
    if (on_shell_warning) {
      std::stringstream ss{};
      ss << four_momentum;
      warnings[1] =
          "Provided 4-momentum " + ss.str() + " [GeV] and mass " +
          std::to_string(mass) + " [GeV] do not satisfy E^2 - p^2 = m^2.\n" +
          "This may originate from the lack of numerical" +
          " precision in the input. Setting E to sqrt(p^2 + " +
          "m^2).\nFurther warnings about E != sqrt(p^2 + m^2) will" +
          " be suppressed.\n" + emph + "Please make sure that setting " +
          "particles back on the mass shell is an acceptable behavior." +
          restore_default;
    }
    return warnings;
  };
  auto warn_if_needed = [&log_area](bool &flag,
                                    const std::optional<std::string> &message) {
    if (flag) {
      logg[log_area].warn(message.value());
      flag = false;
    }
  };
  auto is_particle_stable_and_with_invalid_mass =
      [&mass](const ParticleData &p) {
        return p.type().is_stable() &&
               std::abs(mass - p.pole_mass()) > really_small;
      };
  auto is_particle_off_its_mass_shell = [&mass](const ParticleData &p) {
    return std::abs(p.momentum().sqr() - mass * mass) > really_small;
  };

  // Actual implementation
  ParticleData smash_particle{ParticleType::find(pdgcode)};
  const auto warnings = prepare_needed_warnings(smash_particle);
  if (is_particle_stable_and_with_invalid_mass(smash_particle)) {
    warn_if_needed(mass_warning, warnings[0]);
    smash_particle.set_4momentum(smash_particle.pole_mass(),
                                 four_momentum.threevec());
    //smash_particle.set_unpolarized_spin_vector();
  } else {
    smash_particle.set_4momentum(four_momentum);
    //smash_particle.set_unpolarized_spin_vector();
    if (is_particle_off_its_mass_shell(smash_particle)) {
      warn_if_needed(on_shell_warning, warnings[1]);
      smash_particle.set_4momentum(mass, four_momentum.threevec());
      //smash_particle.set_unpolarized_spin_vector();
    }
  }

  // Set spatial coordinates, they will later be backpropagated if needed
  smash_particle.set_4position(four_position);
  smash_particle.set_formation_time(four_position.x0());
  smash_particle.set_cross_section_scaling_factor(1.0);

  return smash_particle;
}

bool are_particles_identical_at_given_time(const ParticleData &p1,
                                           const ParticleData &p2,
                                           double time) {
  if (p1.pdgcode() != p2.pdgcode()) {
    return false;
  } else {
    if (p1.momentum() != p2.momentum()) {
      return false;
    }
    auto get_propagated_position = [&time](const ParticleData &p) {
      const double t = p.position().x0();
      const FourVector u(1.0, p.velocity());
      return p.position() + u * (time - t);
    };
    return get_propagated_position(p1) == get_propagated_position(p2);
  }
}

}  // namespace smash
