/*
 *
 *    Copyright (c) 2013-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
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
#include "smash/configuration.h"
#include "smash/constants.h"
#include "smash/distributions.h"
#include "smash/experimentparameters.h"
#include "smash/fourvector.h"
#include "smash/hadgas_eos.h"
#include "smash/logging.h"
#include "smash/macros.h"
#include "smash/particles.h"
#include "smash/random.h"
#include "smash/spheremodus.h"
#include "smash/threevector.h"

namespace smash {

/*!\Userguide
 * \page input_modi_sphere_ Sphere
 *
 * \key Radius (double, required): \n
 * Radius of the Sphere, in fm.
 *
 * \key Sphere_Temperature (double, required):\n
 * Temperature to sample momenta in the sphere, in GeV.
 *
 * \key Start_Time (double, required):\n
 * Starting time of sphere calculation.
 *
 * \key Init_Multiplicities (int int, required):\n
 * Initial multiplicities per particle species.
 * Map of PDG number and quantity of this PDG number.
 * Controls how many particles of each species will be initialized.
 *
 * \key Use_Thermal_Multiplicities (bool, optional, default = false): \n
 * If this option is set to true then Init_Multiplicities are ignored and the
 * box is initialized with all particle species of the particle table that
 * belong to the hadron gas equation of state (see
 * HadronGasEos::is_eos_particle()). The multiplicities are sampled from
 * Poisson distributions \f$ Poi(n_i V) \f$, where \f$ n_i \f$ are the
 * grand-canonical thermal densities of the corresponding species and \f$ V \f$
 * is the box volume. This option simulates the grand-canonical ensemble, where
 * the number of particles is not fixed from event to event.
 *
 * \key Baryon_Chemical_Potential (double, optional, default = 0.0): \n
 * Baryon chemical potential \f$ \mu_B \f$ only used if
 * Use_Thermal_Multiplicities is true to compute thermal densities \f$ n_i \f$.
 *
 * \key Strange_Chemical_Potential (double, optional, default = 0.0): \n
 * Strangeness chemical potential \f$ \mu_S \f$ only used if
 * Use_Thermal_Multiplicities is true to compute thermal densities \f$ n_i \f$.
 *
 * \key Initial_Condition (string, optional, default = "thermal momenta") \n
 * Initial distribution to use for momenta of particles. Mainly used in the
 * expanding universe scenario, options are:
 * \li \key thermal momenta - equilibrium Boltzmann distribution
 * \li \key IC_ES - off-equilibrium distribution
 * \li \key IC_1M - off-equilibrium distribution
 * \li \key IC_2M - off-equilibrium distribution
 * \li \key IC_Massive - off-equilibrium distribution
 *
 * See \iref{Bazow:2016oky} and \iref{Tindall:2016try} for further explanations
 * about the different distribution functions.
 *
 * \n
 * Examples: Configuring a Sphere Simulation
 * --------------
 * The following example configures an expanding sphere with a radius of 5 fm
 * at a temperature of 200 MeV. The particles are initialized
 * with thermal momenta at a start time of 0 fm. The particle numbers at
 * initialization are 100 \f$ \pi^+ \f$, 100 \f$ \pi^0 \f$, 100 \f$ \pi^- \f$,
 * 50 protons and 50 neutrons.
 *
 *\verbatim
 Modi:
     Sphere:
         Radius: 5.0
         Sphere_Temperature: 0.2
         Initial_Condition: "thermal momenta"
         Start_Time: 0.0
         Init_Multiplicities:
             211: 100
             111: 100
             -211: 100
             2212: 50
             2112: 50
 \endverbatim
 *
 * It is also possible to initialize a sphere based on
 * thermal multiplicities. This is done via
 *\verbatim
 Modi:
     Box:
         Length: 10.0
         Temperature: 0.2
         Use_Thermal_Multiplicities: True
 \endverbatim
 *
 * \n
 * \note
 * SMASH is shipped with an example configuration file to set up an expanding
 * sphere simulation initialized with predefined initial particle
 * multiplicities. This file is located in /input/sphere. To run SMASH
 * with the provided example configuration for the sphere, execute \n
 * \n
 * \verbatim
    ./smash -i INPUT_DIR/sphere/config.yaml
 \endverbatim
 * \n
 * Where 'INPUT_DIR' needs to be replaced by the path to the input directory
 * ('../input', if the build directory is located in the smash
 * folder).
 *
 */

SphereModus::SphereModus(Configuration modus_config,
                         const ExperimentParameters &)
    : radius_(modus_config.take({"Sphere", "Radius"})),
      sphere_temperature_(modus_config.take({"Sphere", "Sphere_Temperature"})),
      start_time_(modus_config.take({"Sphere", "Start_Time"})),
      use_thermal_(
          modus_config.take({"Sphere", "Use_Thermal_Multiplicities"}, false)),
      mub_(modus_config.take({"Sphere", "Baryon_Chemical_Potential"}, 0.)),
      mus_(modus_config.take({"Sphere", "Strange_Chemical_Potential"}, 0.)),
      init_multipl_(use_thermal_
                        ? std::map<PdgCode, int>()
                        : modus_config.take({"Sphere", "Init_Multiplicities"})
                              .convert_for(init_multipl_)),
      init_distr_(modus_config.take({"Sphere", "Initial_Condition"},
                                    SphereInitialCondition::ThermalMomenta)) {}

/* console output on startup of sphere specific parameters */
std::ostream &operator<<(std::ostream &out, const SphereModus &m) {
  out << "-- Sphere Modus:\nRadius of the sphere: " << m.radius_ << " [fm]"
      << "\nTemperature for momentum sampling: " << m.sphere_temperature_
      << "\nStarting time for Sphere calculation: " << m.start_time_ << '\n';
  if (m.use_thermal_) {
    out << "Thermal multiplicities\n";
  } else {
    for (const auto &p : m.init_multipl_) {
      out << "Particle " << p.first << " initial multiplicity " << p.second
          << '\n';
    }
  }
  return out;
}

/* initial_conditions - sets particle data for @particles */
double SphereModus::initial_conditions(Particles *particles,
                                       const ExperimentParameters &parameters) {
  const auto &log = logger<LogArea::Sphere>();
  FourVector momentum_total(0, 0, 0, 0);
  const double T = this->sphere_temperature_;
  /* Create NUMBER OF PARTICLES according to configuration */
  if (use_thermal_) {
    const double V = 4.0 / 3.0 * M_PI * radius_ * radius_ * radius_;
    if (average_multipl_.empty()) {
      for (const ParticleType &ptype : ParticleType::list_all()) {
        if (HadronGasEos::is_eos_particle(ptype)) {
          const double n = HadronGasEos::partial_density(ptype, T, mub_, mus_);
          average_multipl_[ptype.pdgcode()] = n * V * parameters.testparticles;
        }
      }
    }
    double nb_init = 0.0, ns_init = 0.0;
    for (const auto &mult : average_multipl_) {
      const int thermal_mult_int = random::poisson(mult.second);
      particles->create(thermal_mult_int, mult.first);
      nb_init += mult.second * mult.first.baryon_number();
      ns_init += mult.second * mult.first.strangeness();
      log.debug(mult.first, " initial multiplicity ", thermal_mult_int);
    }
    log.info("Initial hadron gas baryon density ", nb_init);
    log.info("Initial hadron gas strange density ", ns_init);
  } else {
    for (const auto &p : init_multipl_) {
      particles->create(p.second * parameters.testparticles, p.first);
      log.debug("Particle ", p.first, " initial multiplicity ", p.second);
    }
  }
  /* loop over particle data to fill in momentum and position information */
  for (ParticleData &data : *particles) {
    Angles phitheta;
    /* thermal momentum according Maxwell-Boltzmann distribution */
    double momentum_radial, mass = data.pole_mass();
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
      case (SphereInitialCondition::ThermalMomenta):
      default:
        mass = HadronGasEos::sample_mass_thermal(data.type(), 1.0 / T);
        momentum_radial = sample_momenta_from_thermal(T, mass);
        break;
    }
    phitheta.distribute_isotropically();
    log.debug(data.type().name(), "(id ", data.id(), ") radial momentum ",
              momentum_radial, ", direction", phitheta);
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
  /* Make total 3-momentum 0 */
  for (ParticleData &data : *particles) {
    data.set_4momentum(data.momentum().abs(),
                       data.momentum().threevec() -
                           momentum_total.threevec() / particles->size());
  }

  /* Recalculate total momentum */
  momentum_total = FourVector(0, 0, 0, 0);
  for (ParticleData &data : *particles) {
    momentum_total += data.momentum();
    /* IC: debug checks */
    log.debug() << data;
  }
  /* allows to check energy conservation */
  log.info() << "Sphere initial total 4-momentum [GeV]: " << momentum_total;
  return start_time_;
}
}  // namespace smash
