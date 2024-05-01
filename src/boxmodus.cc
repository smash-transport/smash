/*
 *    Copyright (c) 2012-2020,2022-2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#include "smash/boxmodus.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <map>
#include <utility>
#include <vector>

#include "smash/algorithms.h"
#include "smash/angles.h"
#include "smash/constants.h"
#include "smash/cxx17compat.h"
#include "smash/experimentparameters.h"
#include "smash/input_keys.h"
#include "smash/logging.h"
#include "smash/quantumsampling.h"
#include "smash/random.h"
#include "smash/threevector.h"
#include "smash/wallcrossingaction.h"

namespace smash {
static constexpr int LBox = LogArea::Box::id;

/* console output on startup of box specific parameters */
std::ostream &operator<<(std::ostream &out, const BoxModus &m) {
  out << "-- Box Modus:\nSize of the box: (" << m.length_ << " fm)³\n";
  if (m.use_thermal_) {
    out << "Thermal multiplicities "
        << "(T = " << m.temperature_ << " GeV, muB = " << m.mub_
        << " GeV, muS = " << m.mus_ << " GeV, muQ = " << m.muq_ << " GeV)\n";
  } else {
    for (const auto &p : m.init_multipl_) {
      ParticleTypePtr ptype = &ParticleType::find(p.first);
      out << ptype->name() << " initial multiplicity " << p.second << '\n';
    }
  }
  switch (m.initial_condition_) {
    case BoxInitialCondition::PeakedMomenta:
      out << "All initial momenta = 3T = " << 3 * m.temperature_ << " GeV\n";
      break;
    case BoxInitialCondition::ThermalMomentaBoltzmann:
      out << "Boltzmann momentum distribution with T = " << m.temperature_
          << " GeV.\n";
      break;
    case BoxInitialCondition::ThermalMomentaQuantum:
      out << "Fermi/Bose momentum distribution with T = " << m.temperature_
          << " GeV.\n";
      break;
  }
  if (m.jet_pdg_) {
    ParticleTypePtr ptype = &ParticleType::find(m.jet_pdg_.value());
    out << "Adding a " << ptype->name() << " as a jet in the middle "
        << "of the box with " << m.jet_mom_ << " GeV initial momentum.\n";
  }
  return out;
}

BoxModus::BoxModus(Configuration modus_config,
                   const ExperimentParameters &parameters)
    : initial_condition_(
          modus_config.take(InputKeys::modi_box_initialCondition)),
      length_(modus_config.take(InputKeys::modi_box_length)),
      equilibration_time_(
          modus_config.take(InputKeys::modi_box_equilibrationTime)),
      temperature_(modus_config.take(InputKeys::modi_box_temperature)),
      start_time_(modus_config.take(InputKeys::modi_box_startTime)),
      use_thermal_(
          modus_config.take(InputKeys::modi_box_useThermalMultiplicities)),
      mub_(modus_config.take(InputKeys::modi_box_baryonChemicalPotential)),
      mus_(modus_config.take(InputKeys::modi_box_strangeChemicalPotential)),
      muq_(modus_config.take(InputKeys::modi_box_chargeChemicalPotential)),
      account_for_resonance_widths_(
          modus_config.take(InputKeys::modi_box_accountResonanceWidths)),
      init_multipl_(
          use_thermal_
              ? std::map<PdgCode, int>()
              : modus_config.take(InputKeys::modi_box_initialMultiplicities)),
      /* Note that it is crucial not to take other keys from the Jet section
       * before Jet_PDG, since we want here the take to throw in case the user
       * had a Jet section without the mandatory Jet_PDG key. If all other keys
       * are taken first, the section is removed from the config because empty,
       * and has_section(InputSections::m_b_jet) method would return false.
       */
      jet_pdg_(modus_config.has_section(InputSections::m_b_jet)
                   ? make_optional<PdgCode>(
                         modus_config.take(InputKeys::modi_box_jet_jetPdg))
                   : std::nullopt),

      jet_mom_(modus_config.take(InputKeys::modi_box_jet_jetMomentum)) {
  if (parameters.res_lifetime_factor < 0.) {
    throw std::invalid_argument(
        "Resonance lifetime modifier cannot be negative!");
  }
  // Check consistency, just in case
  if (std::abs(length_ - parameters.box_length) > really_small) {
    throw std::runtime_error("Box length inconsistency");
  }
}

double BoxModus::initial_conditions(Particles *particles,
                                    const ExperimentParameters &parameters) {
  double momentum_radial = 0.0, mass = 0.0;
  Angles phitheta;
  FourVector momentum_total(0, 0, 0, 0);
  auto uniform_length = random::make_uniform_distribution(0.0, this->length_);
  const double T = this->temperature_;
  const double V = length_ * length_ * length_;
  /* Create NUMBER OF PARTICLES according to configuration, or thermal case */
  if (use_thermal_) {
    if (average_multipl_.empty()) {
      for (const ParticleType &ptype : ParticleType::list_all()) {
        if (HadronGasEos::is_eos_particle(ptype)) {
          const double lifetime_factor =
              ptype.is_stable() ? 1. : parameters.res_lifetime_factor;
          const double n = lifetime_factor * HadronGasEos::partial_density(
                                                 ptype, T, mub_, mus_, muq_,
                                                 account_for_resonance_widths_);
          average_multipl_[ptype.pdgcode()] = n * V * parameters.testparticles;
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
      logg[LBox].debug(mult.first, " initial multiplicity ", thermal_mult_int);
    }
    logg[LBox].info("Initial hadron gas baryon density ", nb_init);
    logg[LBox].info("Initial hadron gas strange density ", ns_init);
    logg[LBox].info("Initial hadron gas charge density ", nq_init);
  } else {
    for (const auto &p : init_multipl_) {
      particles->create(p.second * parameters.testparticles, p.first);
      logg[LBox].debug("Particle ", p.first, " initial multiplicity ",
                       p.second);
    }
  }
  std::unique_ptr<QuantumSampling> quantum_sampling;
  if (this->initial_condition_ == BoxInitialCondition::ThermalMomentaQuantum) {
    quantum_sampling = std::make_unique<QuantumSampling>(init_multipl_, V, T);
  }
  for (ParticleData &data : *particles) {
    /* Set MOMENTUM SPACE distribution */
    if (this->initial_condition_ == BoxInitialCondition::PeakedMomenta) {
      /* initial thermal momentum is the average 3T */
      momentum_radial = 3.0 * T;
      mass = data.pole_mass();
    } else {
      if (this->initial_condition_ ==
          BoxInitialCondition::ThermalMomentaBoltzmann) {
        /* thermal momentum according Maxwell-Boltzmann distribution */
        mass = (!account_for_resonance_widths_)
                   ? data.type().mass()
                   : HadronGasEos::sample_mass_thermal(data.type(), 1.0 / T);
        momentum_radial = sample_momenta_from_thermal(T, mass);
      } else if (this->initial_condition_ ==
                 BoxInitialCondition::ThermalMomentaQuantum) {
        /*
         * Sampling the thermal momentum according Bose/Fermi/Boltzmann
         * distribution.
         * We take the pole mass as the mass.
         */
        mass = data.type().mass();
        momentum_radial = quantum_sampling->sample(data.pdgcode());
      }
    }
    phitheta.distribute_isotropically();
    logg[LBox].debug(data.type().name(), "(id ", data.id(),
                     ") radial momentum ", momentum_radial, ", direction",
                     phitheta);
    data.set_4momentum(mass, phitheta.threevec() * momentum_radial);
    momentum_total += data.momentum();

    /* Set COORDINATE SPACE distribution */
    ThreeVector pos{uniform_length(), uniform_length(), uniform_length()};
    data.set_4position(FourVector(start_time_, pos));
    /// Initialize formation time
    data.set_formation_time(start_time_);
  }

  /* Make total 3-momentum of the box 0 and initialize an unpolarized spin
   * vector */
  for (ParticleData &data : *particles) {
    data.set_4momentum(data.momentum().abs(),
                       data.momentum().threevec() -
                           momentum_total.threevec() / particles->size());
    // Initialize spin vector
    data.set_unpolarized_spin_vector();
  }

  /* Add a single highly energetic particle in the center of the box (jet) */
  if (jet_pdg_) {
    auto &jet_particle = particles->create(jet_pdg_.value());
    jet_particle.set_formation_time(start_time_);
    jet_particle.set_4position(FourVector(start_time_, 0., 0., 0.));
    jet_particle.set_4momentum(ParticleType::find(jet_pdg_.value()).mass(),
                               ThreeVector(jet_mom_, 0., 0.));
    jet_particle.set_unpolarized_spin_vector();
  }

  /* Recalculate total momentum */
  momentum_total = FourVector(0, 0, 0, 0);
  for (ParticleData &data : *particles) {
    momentum_total += data.momentum();
    /* IC: debug checks */
    logg[LBox].debug() << data;
  }
  /* allows to check energy conservation */
  logg[LBox].debug() << "Initial total 4-momentum [GeV]: " << momentum_total;
  return start_time_;
}

int BoxModus::impose_boundary_conditions(Particles *particles,
                                         const OutputsList &output_list) {
  int wraps = 0;

  for (ParticleData &data : *particles) {
    FourVector position = data.position();
    bool wall_hit = enforce_periodic_boundaries(position.begin() + 1,
                                                position.end(), length_);
    if (wall_hit) {
      const ParticleData incoming_particle(data);
      data.set_4position(position);
      ++wraps;
      ActionPtr action =
          std::make_unique<WallcrossingAction>(incoming_particle, data);
      for (const auto &output : output_list) {
        if (!output->is_dilepton_output() && !output->is_photon_output()) {
          output->at_interaction(*action, 0.);
        }
      }
    }
  }
  logg[LBox].debug("Moved ", wraps, " particles back into the box.");
  return wraps;
}

}  // namespace smash
