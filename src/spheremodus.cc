/*
 *
 *    Copyright (c) 2013-2014
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

#include "include/algorithms.h"
#include "include/angles.h"
#include "include/configuration.h"
#include "include/constants.h"
#include "include/distributions.h"
#include "include/experimentparameters.h"
#include "include/fourvector.h"
#include "include/hadgas_eos.h"
#include "include/logging.h"
#include "include/macros.h"
#include "include/particles.h"
#include "include/random.h"
#include "include/spheremodus.h"
#include "include/threevector.h"

namespace Smash {

/*!\Userguide
 * \page input_modi_sphere_ Sphere
 *
 * \key Radius (float, required): \n
 * Radius of the Sphere.
 *
 * \key Sphere_Temperature (float, required):\n
 * Temperature for the momentum sampling in the sphere in GeV.
 *
 * \key Start_Time (float, required):\n
 * Starting time of Sphere calculation.
 *
 * \key Init_Multiplicities (int int, required):\n
 * Initial multiplicities per particle species.
 * Map of PDG number and quantity of this PDG number.
 * Controls how many particles of each sort will be initialized. \n
 * Example:
 * \verbatim
 Init_Multiplicities:
 2112: 200
 -2112: 100
 \endverbatim
 * It means that 200 neutrons and 100 antineutrons will be initialized.
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
 * Baryon chemical potential \f$ \mu_B \f$ used in case if
 * Use_Thermal_Multiplicities is true to compute thermal densities \f$ n_i \f$.
 *
 * \key Strange_Chemical_Potential (double, optional, default = 0.0): \n
 * Strangeness chemical potential \f$ \mu_S \f$ used in case if
 * Use_Thermal_Multiplicities is true to compute thermal densities \f$ n_i \f$.
 */


SphereModus::SphereModus(Configuration modus_config,
                         const ExperimentParameters &)
    : radius_(modus_config.take({"Sphere", "Radius"})),
      sphere_temperature_(modus_config.take({"Sphere", "Sphere_Temperature"})),
      start_time_(modus_config.take({"Sphere", "Start_Time"})),
      use_thermal_(
        modus_config.take({"Sphere", "Use_Thermal_Multiplicities"}, false)),
      mub_(modus_config.take({"Sphere", "Baryon_Chemical_Potential"}, 0.0f)),
      mus_(modus_config.take({"Sphere", "Strange_Chemical_Potential"}, 0.0f)),
      init_multipl_(use_thermal_ ? std::map<PdgCode, int>() :
                    modus_config.take({"Sphere", "Init_Multiplicities"}).
                    convert_for(init_multipl_)),
      init_distr_(modus_config.take({"Sphere", "Initial_Condition"},
                    SphereInitialCondition::ThermalMomenta)) {
}

/* console output on startup of sphere specific parameters */
std::ostream &operator<<(std::ostream &out, const SphereModus &m) {
  out << "-- Sphere Modus:\nRadius of the sphere: " << m.radius_ << " [fm]"
      << "\nTemperature for momentum sampling: " << m.sphere_temperature_
      << "\nStarting time for Sphere calculation: " << m.start_time_ << '\n';
  if (m.use_thermal_) {
    out << "Thermal multiplicities\n";
  } else {
    for (const auto &p : m.init_multipl_) {
      out << "Particle " << p.first << " initial multiplicity "
                         << p.second << '\n';
    }
  }
  return out;
}

/* initial_conditions - sets particle data for @particles */
float SphereModus::initial_conditions(Particles *particles,
  const ExperimentParameters &parameters) {
  const auto &log = logger<LogArea::Sphere>();
  FourVector momentum_total(0, 0, 0, 0);
  /* Create NUMBER OF PARTICLES according to configuration */
  if (use_thermal_) {
    const double T = sphere_temperature_;
    const double V = 4.0/3.0 * M_PI * radius_*radius_*radius_;
    for (const ParticleType &ptype : ParticleType::list_all()) {
      if (HadronGasEos::is_eos_particle(ptype)) {
        const double n = HadronGasEos::partial_density(ptype, T, mub_, mus_);
        const double thermal_mult = n*V*parameters.testparticles;
        assert(thermal_mult > 0.0);
        const int thermal_mult_int = Random::poisson(thermal_mult);
        particles->create(thermal_mult_int, ptype.pdgcode());
        log.debug(ptype.name(), " initial multiplicity ", thermal_mult_int);
      }
    }
    log.info() << "Initial baryon density "
               << HadronGasEos::net_baryon_density(T, mub_, mus_);
    log.info() << "Initial strange density "
               << HadronGasEos::net_strange_density(T, mub_, mus_);
  } else {
    for (const auto &p : init_multipl_) {
      particles->create(p.second*parameters.testparticles, p.first);
      log.debug() << "Particle " << p.first
                  << " initial multiplicity " << p.second;
    }
  }
  /* loop over particle data to fill in momentum and position information */
  for (ParticleData &data : *particles) {
    Angles phitheta;
    /* thermal momentum according Maxwell-Boltzmann distribution */
    double momentum_radial;
    /* assign momentum_radial according to requested distribution */
    switch (init_distr_) {
      case (SphereInitialCondition::ThermalMomenta):
        momentum_radial = sample_momenta_from_thermal(this->sphere_temperature_,
                                                      data.pole_mass());
        break;
      case (SphereInitialCondition::IC_ES):
        momentum_radial = sample_momenta_IC_ES(this->sphere_temperature_);
        break;
      case (SphereInitialCondition::IC_1M):
        momentum_radial = sample_momenta_IC_1M(this->sphere_temperature_,
                                                      data.pole_mass());
        break;
      case (SphereInitialCondition::IC_2M):
        momentum_radial = sample_momenta_IC_2M(this->sphere_temperature_,
                                                      data.pole_mass());
        break;
      case (SphereInitialCondition::IC_Massive):
        momentum_radial = sample_momenta_non_eq_mass(this->sphere_temperature_,
                                                      data.pole_mass());
        break;
      default:
        momentum_radial = sample_momenta_from_thermal(this->sphere_temperature_,
                                                      data.pole_mass());
        break;
    }
    phitheta.distribute_isotropically();
    log.debug("Particle ", data.id(), " radial momenta ", momentum_radial, ' ',
              phitheta);
    data.set_4momentum(data.pole_mass(), phitheta.threevec() * momentum_radial);
    momentum_total += data.momentum();
    /* uniform sampling in a sphere with radius r */
    double position_radial;
    position_radial = std::cbrt(Random::canonical()) * radius_;
    Angles pos_phitheta;
    pos_phitheta.distribute_isotropically();
    data.set_4position(FourVector(start_time_,
                                  pos_phitheta.threevec() * position_radial));
    data.set_formation_time(start_time_);
  }
  /* Make total 3-momentum 0 */
  for (ParticleData &data : *particles) {
    data.set_4momentum(data.pole_mass(), data.momentum().threevec() -
                       momentum_total.threevec()/particles->size());
  }

  /* Recalculate total momentum */
  momentum_total = FourVector(0, 0, 0, 0);
  for (ParticleData &data : *particles) {
    momentum_total += data.momentum();
    /* IC: debug checks */
    log.debug() << data;
  }
  /* allows to check energy conservation */
  log.info() << "Sphere initial total 4-momentum [GeV]: "
             << momentum_total;
  return start_time_;
}
}  // namespace Smash
