/*
 *    Copyright (c) 2012-2020,2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <map>
#include <utility>
#include <vector>

#include "smash/algorithms.h"
#include "smash/angles.h"
#include "smash/boxmodus.h"
#include "smash/constants.h"
#include "smash/experimentparameters.h"
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
  if (m.insert_jet_) {
    ParticleTypePtr ptype = &ParticleType::find(m.jet_pdg_);
    out << "Adding a " << ptype->name() << " as a jet in the middle "
        << "of the box with " << m.jet_mom_ << " GeV initial momentum.\n";
  }
  return out;
}

/*!\Userguide
 * \page input_modi_box_ Box
 *
 * \key Initial_Condition (string, required, no default): \n
 * Controls initial momentum distribution of particles.
 * \li \key "peaked momenta" - All particles have momentum \f$p = 3 \cdot T\f$,
 * where T is the temperature. Directions of momenta are uniformly distributed.
 * \li \key "thermal momenta" - Momenta are sampled from a Maxwell-Boltzmann
 * distribution.
 *
 * \key Length (double, required): \n
 * Length of the cube's edge, in fm.
 *
 * \key Temperature (double, required): \n
 * Temperature in the box, in GeV.
 *
 * \key Start_Time (double, required): \n
 * Starting time of the simulation.
 * All particles in the box are initialized with \f$x^0\f$ = Start_Time.
 *
 * \key Equilibration_Time (double, optional): \n
 * Time after which the output of the box is written out. The first time however
 * will be printed. This is useful if one wants to simulate boxes for very long
 * times and knows at which time the box reaches its thermal and chemical
 * equilibrium.
 *
 * \key Init_Multiplicities (int, required): \n
 * Map of PDG number and quantity of this PDG number.
 * Controls how many particles of each sort will be initialized.
 *
 * \key Use_Thermal_Multiplicities (bool, optional, default = false): \n
 * If this option is set to true then Init_Multiplicities are ignored and the
 * box is initialized with all particle species of the particle table that
 * belong to the hadron gas equation of state, see
 * HadronGasEos::is_eos_particle(). The multiplicities are sampled from
 * Poisson distributions \f$ Poi(n_i V) \f$, where \f$ n_i \f$ are the
 * grand-canonical thermal densities of the corresponding species and \f$ V \f$
 * is the box volume. This option simulates the grand-canonical ensemble, where
 * the number of particles is not fixed from event to event.
 *
 * \key Baryon_Chemical_Potential (double, optional, default = 0.0): \n
 * Baryon chemical potential \f$ \mu_B \f$ used in case if
 * Use_Thermal_Multiplicities is true to compute thermal densities \f$ n_i \f$.
 *
 * \key Strange_Chemical_Potential (double, optional, default = 0.0): \n
 * Strangeness chemical potential \f$ \mu_S \f$ used in case if
 * Use_Thermal_Multiplicities is true to compute thermal densities \f$ n_i \f$.
 *
 * \key Charge_Chemical_Potential (double, optional, default = 0.0): \n
 * Charge chemical potential \f$ \mu_Q \f$ used in case if
 * Use_Thermal_Multiplicities is true to compute thermal densities \f$ n_i \f$.
 *
 * \key Account_Resonance_Widths (bool, optional, default = true): \n
 * In case of thermal initialization: true -- account for resonance
 * spectral functions, while computing multiplicities and sampling masses,
 * false -- simply use pole masses.
 *
 * Normally, one wants this option true. For example for the detailed balance
 * studies it is better to account for spectral functions, because then at t =
 * 0 one has exactly the expected thermal grand-canonical multiplicities, that
 * can be compared to final ones.  However, by switching true/false one can
 * observe the effect of spectral function on the multiplicity. This is useful
 * for understanding the implications of different ways of sampling resonances
 * in hydrodynamics.
 *
 * \key Jet: \n
 * This subset of config values is used to put a single high energy particle
 * (a "jet") in the center of the box, on an trajectory along
 * the x axis; if no pdg is specified no jet is produced.
 *
 * \li \key Jet_PDG (int, optional):
 * The type of particle to be used as a jet, as given by its PDG code;
 * if none is provided no jet is initialized.
 *
 * \li \key Jet_Momentum (double, optional, default = 20.):
 * The initial momentum to give to the jet particle (in GeV)
 *
 * \n
 * Examples: Configuring a Box Simulation
 * --------------
 * The following example configures an infinite matter simulation in a Box with
 * 10 fm cube length at a temperature of 200 MeV. The particles are initialized
 * with thermal momenta at a start time of 10 fm. The particle numbers at
 * initialization are 100 \f$ \pi^+ \f$, 100 \f$ \pi^0 \f$, 100 \f$ \pi^- \f$,
 * 50 protons and 50 neutrons.
 *
 *\verbatim
 Modi:
     Box:
         Length: 10.0
         Temperature: 0.2
         Initial_Condition: "thermal momenta"
         Start_Time: 10.0
         Init_Multiplicities:
             211: 100
             111: 100
             -211: 100
             2212: 50
             2112: 50
 \endverbatim
 *
 * On the contrary, it is also possible to initialize a thermal box based on
 * thermal multiplicities. This is done via
 *\verbatim
 Modi:
     Box:
         Length: 10.0
         Temperature: 0.2
         Use_Thermal_Multiplicities: True
         Initial_Condition: "thermal momenta"
         Baryon_Chemical_Potential: 0.0
         Strange_Chemical_Potential: 0.0
         Charge_Chemical_Potential: 0.0
         Account_Resonance_Widths: True
 \endverbatim
 *
 * If one wants to simulate a jet in the hadronic medium, this can be done
 * by using the following configuration setup:
 *\verbatim
 Modi:
     Box:
         Length: 10.0
         Temperature: 0.2
         Use_Thermal_Multiplicities: True
         Initial_Condition: "thermal momenta"
         Jet:
             Jet_PDG: 211
             Jet_Momentum: 100.0
\endverbatim
 * \n
 *
 * \note
 * The box modus is most useful for infinite matter simulations
 * with thermal and chemical equilibration and detailed balance. Detailed
 * balance can however not be conserved if 3-body decays (or higher) are
 * performed. To yield useful results applying a SMASH box simulation, it is
 * therefore necessary to modify the provided default particles.txt and
 * decaymodes.txt by removing 3-body and higher order decays from
 * the decaymodes file and all corresponding particles that can no longer be
 * produced from the particles file. In addtion, strings need to be
 * turned off, since they also break detailed balance due to lacking
 * backreactions. \n
 * SMASH is shipped with example files (config.yaml, particles.txt,
 * decaymodes.txt) meeting the above mentioned requirements to set up an
 * infinite matter simulation. They are located in /input/box. To run SMASH
 * with the provided example files, execute \n
 * \n
 * \verbatim
    ./smash -i INPUT_DIR/box/config.yaml -p INPUT_DIR/box/particles.txt -d
 INPUT_DIR/box/decaymodes.txt \endverbatim
 * \n
 * Where 'INPUT_DIR' needs to be replaced by the path to the input directory
 * ('../input', if the build directory is located in the smash
 * folder).
 */
BoxModus::BoxModus(Configuration modus_config,
                   const ExperimentParameters &parameters)
    : initial_condition_(modus_config.take({"Box", "Initial_Condition"})),
      length_(modus_config.take({"Box", "Length"})),
      equilibration_time_(
          modus_config.take({"Box", "Equilibration_Time"}, -1.)),
      temperature_(modus_config.take({"Box", "Temperature"})),
      start_time_(modus_config.take({"Box", "Start_Time"}, 0.)),
      use_thermal_(
          modus_config.take({"Box", "Use_Thermal_Multiplicities"}, false)),
      mub_(modus_config.take({"Box", "Baryon_Chemical_Potential"}, 0.)),
      mus_(modus_config.take({"Box", "Strange_Chemical_Potential"}, 0.)),
      muq_(modus_config.take({"Box", "Charge_Chemical_Potential"}, 0.)),
      account_for_resonance_widths_(
          modus_config.take({"Box", "Account_Resonance_Widths"}, true)),
      init_multipl_(use_thermal_
                        ? std::map<PdgCode, int>()
                        : modus_config.take({"Box", "Init_Multiplicities"})
                              .convert_for(init_multipl_)),
      insert_jet_(modus_config.has_value({"Box", "Jet", "Jet_PDG"})),
      jet_pdg_(insert_jet_ ? modus_config.take({"Box", "Jet", "Jet_PDG"})
                                 .convert_for(jet_pdg_)
                           : pdg::p),  // dummy default; never used
      jet_mom_(modus_config.take({"Box", "Jet", "Jet_Momentum"}, 20.)) {
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

  /* Make total 3-momentum 0 */
  for (ParticleData &data : *particles) {
    data.set_4momentum(data.momentum().abs(),
                       data.momentum().threevec() -
                           momentum_total.threevec() / particles->size());
  }

  /* Add a single highly energetic particle in the center of the box (jet) */
  if (insert_jet_) {
    auto &jet_particle = particles->create(jet_pdg_);
    jet_particle.set_formation_time(start_time_);
    jet_particle.set_4position(FourVector(start_time_, 0., 0., 0.));
    jet_particle.set_4momentum(ParticleType::find(jet_pdg_).mass(),
                               ThreeVector(jet_mom_, 0., 0.));
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
