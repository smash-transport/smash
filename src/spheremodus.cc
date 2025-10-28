/*
 *
 *    Copyright (c) 2012-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/spheremodus.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <map>
#include <utility>
#include <vector>

#include "smash/angles.h"
#include "smash/chemicalpotential.h"
#include "smash/configuration.h"
#include "smash/constants.h"
#include "smash/cxx17compat.h"
#include "smash/experimentparameters.h"
#include "smash/fourvector.h"
#include "smash/hadgas_eos.h"
#include "smash/logging.h"
#include "smash/particles.h"
#include "smash/quantumsampling.h"
#include "smash/random.h"
#include "smash/threevector.h"

namespace smash {
static constexpr int LSphere = LogArea::Sphere::id;

SphereModus::SphereModus(Configuration modus_config,
                         const ExperimentParameters &)
    : radius_(modus_config.take(InputKeys::modi_sphere_radius)),
      sphere_temperature_(
          modus_config.take(InputKeys::modi_sphere_temperature)),
      start_time_(modus_config.take(InputKeys::modi_sphere_startTime)),
      use_thermal_(
          modus_config.take(InputKeys::modi_sphere_useThermalMultiplicities)),
      mub_(modus_config.take(InputKeys::modi_sphere_baryonChemicalPotential)),
      mus_(modus_config.take(InputKeys::modi_sphere_strangeChemicalPotential)),
      muq_(modus_config.take(InputKeys::modi_sphere_chargeChemicalPotential)),
      hf_multiplier_(
          modus_config.take(InputKeys::modi_sphere_heavyFlavorMultiplier)),
      account_for_resonance_widths_(
          modus_config.take(InputKeys::modi_sphere_accountResonanceWidths)),
      init_multipl_(use_thermal_
                        ? std::map<PdgCode, int>()
                        : modus_config.take(
                              InputKeys::modi_sphere_initialMultiplicities)),
      init_distr_(modus_config.take(InputKeys::modi_sphere_initialCondition)),
      radial_velocity_(
          modus_config.take(InputKeys::modi_sphere_addRadialVelocity)),
      radial_velocity_exponent_(
          modus_config.take(InputKeys::modi_sphere_addRadialVelocityExponent)),
      /* Note that it is crucial not to take other keys from the Jet section
       * before Jet_PDG, since we want here the take to throw in case the user
       * had a Jet section without the mandatory Jet_PDG key. If all other keys
       * are taken first, the section is removed from the config because empty,
       * and has_section(InputSections::m_s_jet) method would return false.
       */
      jet_pdg_(modus_config.has_section(InputSections::m_s_jet)
                   ? make_optional<PdgCode>(
                         modus_config.take(InputKeys::modi_sphere_jet_jetPdg))
                   : std::nullopt),
      jet_mom_(modus_config.take(InputKeys::modi_sphere_jet_jetMomentum)),
      jet_pos_(modus_config.take(InputKeys::modi_sphere_jet_jetPosition)),
      jet_back_(modus_config.take(InputKeys::modi_sphere_jet_backToBack)),
      jet_back_separation_(
          jet_back_ ? modus_config.take(
                          InputKeys::modi_sphere_jet_backToBackSeparation)
                    : 0) {
  if (!jet_back_ &&
      modus_config.has_value(InputKeys::modi_sphere_jet_backToBackSeparation)) {
    throw std::invalid_argument(
        "In order to specify 'Back_To_Back_Separation', 'Back_To_Back' must be "
        "true.");
  }
}

/* console output on startup of sphere specific parameters */
std::ostream &operator<<(std::ostream &out, const SphereModus &m) {
  out << "-- Sphere Modus:\nRadius of the sphere: " << m.radius_ << " fm\n";
  if (m.use_thermal_) {
    out << "Thermal multiplicities (T = " << m.sphere_temperature_
        << " GeV, muB = " << m.mub_ << " GeV, muS = " << m.mus_
        << " GeV, muQ = " << m.muq_ << " GeV)\n";
  } else {
    for (const auto &p : m.init_multipl_) {
      ParticleTypePtr ptype = &ParticleType::find(p.first);
      out << ptype->name() << " initial multiplicity " << p.second << '\n';
    }
  }
  switch (m.init_distr_) {
    case SphereInitialCondition::ThermalMomentaBoltzmann:
      out << "Boltzmann momentum distribution with T = "
          << m.sphere_temperature_ << " GeV.\n";
      break;
    case SphereInitialCondition::ThermalMomentaQuantum:
      out << "Fermi/Bose momentum distribution with T = "
          << m.sphere_temperature_ << " GeV.\n";
      break;
    case SphereInitialCondition::IC_ES:
      out << "Sphere Initial Condition is IC_ES";
      break;
    case SphereInitialCondition::IC_1M:
      out << "Sphere Initial Condition is IC_1M";
      break;
    case SphereInitialCondition::IC_2M:
      out << "Sphere Initial Condition is IC_2M";
      break;
    case SphereInitialCondition::IC_Massive:
      out << "Sphere Initial Condition is IC_Massive";
      break;
  }
  if (m.jet_pdg_) {
    ParticleTypePtr ptype = &ParticleType::find(m.jet_pdg_.value());
    const auto pos = m.jet_pos_;
    if (m.jet_back_) {
      ParticleTypePtr anti =
          ptype->has_antiparticle() ? ptype->get_antiparticle() : ptype;
      out << "Adding a dijet " << ptype->name() << anti->name()
          << " centered at (" << pos.x1() << ", " << pos.x2() << ", "
          << pos.x3() << ") separated by " << m.jet_back_separation_
          << " fm,\neach with " << m.jet_mom_ << " GeV of initial momentum.\n";
    } else {
      out << "Adding a " << ptype->name() << " as a jet at (" << pos.x1()
          << ", " << pos.x2() << ", " << pos.x3() << ") fm with " << m.jet_mom_
          << " GeV of initial momentum.\n";
    }
  }
  return out;
}

/* initial_conditions - sets particle data for @particles */
double SphereModus::initial_conditions(Particles *particles,
                                       const ExperimentParameters &parameters) {
  FourVector momentum_total(0, 0, 0, 0);
  const double T = this->sphere_temperature_;
  const double V = 4.0 / 3.0 * M_PI * radius_ * radius_ * radius_;
  /* Create NUMBER OF PARTICLES according to configuration */
  if (use_thermal_) {
    if (average_multipl_.empty()) {
      for (const ParticleType &ptype : ParticleType::list_all()) {
        const bool is_eos_particle = HadronGasEos::is_eos_particle(ptype);
        const bool use_heavy_flavor = ptype.pdgcode().is_heavy_flavor() &&
                                      (hf_multiplier_ > really_small);
        if (is_eos_particle || use_heavy_flavor) {
          const double n = HadronGasEos::partial_density(
              ptype, T, mub_, mus_, muq_, account_for_resonance_widths_);
          average_multipl_[ptype.pdgcode()] = n * V * parameters.testparticles;
          if (ptype.pdgcode().is_heavy_flavor()) {
            average_multipl_[ptype.pdgcode()] *= hf_multiplier_;
          }
        }
      }
    }
    double nb_init = 0.0, ns_init = 0.0, nq_init = 0.0;
    for (const auto &mult : average_multipl_) {
      const int thermal_mult_int = random::poisson(mult.second);
      particles->create(thermal_mult_int, mult.first);
      nb_init += mult.second * mult.first.baryon_number();
      ns_init += mult.second * mult.first.strangeness();
      nq_init += mult.second * mult.first.charge();
      logg[LSphere].debug(mult.first, " initial multiplicity ",
                          thermal_mult_int);
    }
    logg[LSphere].info("Initial hadron gas baryon density ", nb_init);
    logg[LSphere].info("Initial hadron gas strange density ", ns_init);
    logg[LSphere].info("Initial hadron gas charge density ", nq_init);
  } else {
    for (const auto &p : init_multipl_) {
      particles->create(p.second * parameters.testparticles, p.first);
      logg[LSphere].debug("Particle ", p.first, " initial multiplicity ",
                          p.second);
    }
  }
  std::unique_ptr<QuantumSampling> quantum_sampling;
  if (this->init_distr_ == SphereInitialCondition::ThermalMomentaQuantum) {
    quantum_sampling = std::make_unique<QuantumSampling>(init_multipl_, V, T);
  }
  /* loop over particle data to fill in momentum and position information */
  for (ParticleData &data : *particles) {
    Angles phitheta;
    /* thermal momentum according Maxwell-Boltzmann distribution */
    double momentum_radial = 0.0, mass = data.pole_mass();
    /* assign momentum_radial according to requested distribution */
    switch (init_distr_) {
      case (SphereInitialCondition::IC_ES):
        momentum_radial = sample_momenta_IC_ES(T);
        break;
      case (SphereInitialCondition::IC_1M):
        momentum_radial = sample_momenta_1M_IC(T, mass);
        break;
      case (SphereInitialCondition::IC_2M):
        momentum_radial = sample_momenta_2M_IC(T, mass);
        break;
      case (SphereInitialCondition::IC_Massive):
        momentum_radial = sample_momenta_non_eq_mass(T, mass);
        break;
      case (SphereInitialCondition::ThermalMomentaBoltzmann):
      default:
        mass = (!account_for_resonance_widths_)
                   ? data.type().mass()
                   : HadronGasEos::sample_mass_thermal(data.type(), 1.0 / T);
        momentum_radial = sample_momenta_from_thermal(T, mass);
        break;
      case (SphereInitialCondition::ThermalMomentaQuantum):
        /*
         * **********************************************************************
         * Sampling the thermal momentum according Bose/Fermi/Boltzmann
         * distribution.
         * We take the pole mass as the mass.
         * **********************************************************************
         */
        mass = data.type().mass();
        momentum_radial = quantum_sampling->sample(data.pdgcode());
        break;
    }
    phitheta.distribute_isotropically();
    logg[LSphere].debug(data.type().name(), "(id ", data.id(),
                        ") radial momentum ", momentum_radial, ", direction",
                        phitheta);
    data.set_4momentum(mass, phitheta.threevec() * momentum_radial);
    momentum_total += data.momentum();
    /* uniform sampling in a sphere with radius r */
    double position_radial;
    position_radial = std::cbrt(random::canonical()) * radius_;
    Angles pos_phitheta;
    pos_phitheta.distribute_isotropically();
    data.set_4position(
        FourVector(start_time_, pos_phitheta.threevec() * position_radial));
    data.set_formation_time(start_time_);
  }

  /* Boost in radial direction with an underlying velocity field of the form
   * u_r = u_0 * (r / R)^n
   */
  if (radial_velocity_ > 0.0) {
    if (radial_velocity_ > 1.0) {
      throw std::invalid_argument(
          "Additional velocity cannot be greater than 1!");
    }
    if (radial_velocity_exponent_ < 0.0) {
      throw std::invalid_argument("Flow velocity exponent cannot be negative!");
    }
    for (ParticleData &data : *particles) {
      double particle_radius = std::sqrt(data.position().sqr3());
      auto e_r = data.position().threevec() / particle_radius;
      auto radial_velocity =
          -1.0 * radial_velocity_ * e_r *
          std::pow(particle_radius / radius_, radial_velocity_exponent_);
      data.set_4momentum(data.momentum().lorentz_boost(radial_velocity));
      momentum_total += data.momentum();
    }
  }

  /* Make total 3-momentum 0 */
  for (ParticleData &data : *particles) {
    data.set_4momentum(data.momentum().abs(),
                       data.momentum().threevec() -
                           momentum_total.threevec() / particles->size());
  }

  /* Add a single highly energetic particle in the center of the sphere (jet) */
  if (jet_pdg_) {
    auto &pdg = jet_pdg_.value();
    auto &jet_particle = particles->create(pdg);
    auto displacement = ThreeVector(jet_back_separation_ / 2., 0., 0.);
    jet_particle.set_formation_time(start_time_);
    jet_particle.set_4position(
        FourVector(start_time_, jet_pos_ + displacement));
    jet_particle.set_4momentum(ParticleType::find(pdg).mass(),
                               ThreeVector(jet_mom_, 0., 0.));
    if (jet_back_) {
      auto &anti_pdg = pdg.has_antiparticle() ? pdg.get_antiparticle() : pdg;
      auto &jet_antiparticle = particles->create(anti_pdg);
      jet_antiparticle.set_formation_time(start_time_);
      jet_antiparticle.set_4position(
          FourVector(start_time_, jet_pos_ - displacement));
      jet_antiparticle.set_4momentum(ParticleType::find(anti_pdg).mass(),
                                     ThreeVector(-jet_mom_, 0., 0.));
    }
  }

  /* Recalculate total momentum */
  momentum_total = FourVector(0, 0, 0, 0);
  for (ParticleData &data : *particles) {
    momentum_total += data.momentum();
    /* IC: debug checks */
    logg[LSphere].debug() << data;
  }
  /* allows to check energy conservation */
  logg[LSphere].debug() << "Sphere initial total 4-momentum [GeV]: "
                        << momentum_total;
  return start_time_;
}
}  // namespace smash
