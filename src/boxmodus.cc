/*
 *    Copyright (c) 2013-2014
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

#include "include/algorithms.h"
#include "include/angles.h"
#include "include/boxmodus.h"
#include "include/configuration.h"
#include "include/constants.h"
#include "include/distributions.h"
#include "include/experimentparameters.h"
#include "include/logging.h"
#include "include/macros.h"
#include "include/outputinterface.h"
#include "include/particles.h"
#include "include/random.h"
#include "include/threevector.h"

namespace Smash {

/* console output on startup of box specific parameters */
std::ostream &operator<<(std::ostream &out, const BoxModus &m) {
  return out << "-- Box Modus:\n"
                "Size of the box: (" << m.length_ << " fm)³"
             << "\nInitial temperature: " << m.temperature_ << " GeV"
             << "\nIC type " << m.initial_condition_ << '\n';
}

/*!\Userguide
 * \page input_modi_box_ Box
 *
 * \key INITIAL_CONDITION: \n
 * Controls whether particles are created with
 * thermal momenta (sampled from a Maxwell-Boltzmann distribution) or
 * with average momentum \f$p = 3 \cdot T\f$ with T the temperature. The
 * latter is chosen if \key INITIAL_CONDITION is 2.
 *
 * \key LENGTH: \n
 * Length of the cube's edge in fm/c
 *
 * \key TEMPERATURE: \n
 * Temperature in the box in GeV.
 */
BoxModus::BoxModus(Configuration modus_config, const ExperimentParameters &)
    : initial_condition_(modus_config.take({"Box", "INITIAL_CONDITION"})),
                 length_(modus_config.take({"Box", "LENGTH"})),
            temperature_(modus_config.take({"Box", "TEMPERATURE"})),
             start_time_(modus_config.take({"Box", "START_TIME"})),
           init_multipl_(modus_config.take({"Box", "INIT_MULTIPLICITIES"}).
                                                convert_for(init_multipl_)) {
}

/* initial_conditions - sets particle data for @particles */
float BoxModus::initial_conditions(Particles *particles,
                                  const ExperimentParameters &parameters) {
  const auto &log = logger<LogArea::Box>();
  double momentum_radial = 0;
  Angles phitheta;
  FourVector momentum_total(0, 0, 0, 0);
  auto uniform_length = Random::make_uniform_distribution(0.0,
                                         static_cast<double>(this->length_));

  /* Create NUMBER OF PARTICLES according to configuration */
  for (const auto &p : init_multipl_) {
    particles->create(p.second*parameters.testparticles, p.first);
    log.debug() << "Particle " << p.first << " init multiplicity " << p.second;
  }
  number_density_initial_ = particles->size()/(length_*length_*length_);

  for (ParticleData &data : particles->data()) {
    /* Set MOMENTUM SPACE distribution */
    if (this->initial_condition_ != 2) {
      /* thermal momentum according Maxwell-Boltzmann distribution */
      momentum_radial = sample_momenta(this->temperature_, data.pole_mass());
    } else {
      /* IC == 2 initial thermal momentum is the average 3T */
      momentum_radial = 3.0 * this->temperature_;
    }
    phitheta.distribute_isotropically();
    log.debug() << data << ", radial mom:" << field << momentum_radial << ", "
                << phitheta;
    data.set_4momentum(data.pole_mass(), phitheta.threevec() * momentum_radial);
    momentum_total += data.momentum();

    /* Set COORDINATE SPACE distribution */
    ThreeVector pos{uniform_length(), uniform_length(), uniform_length()};
    data.set_4position(FourVector(start_time_, pos));

    /* IC: debug checks */
    log.debug() << data;
  }

  /* Make total 3-momentum 0 */
  for (ParticleData &data : particles->data()) {
    data.set_4momentum(data.pole_mass(), data.momentum().threevec() -
                       momentum_total.threevec()/particles->size());
  }

  /* Recalculate total momentum */
  momentum_total = FourVector(0, 0, 0, 0);
  for (ParticleData &data : particles->data()) {
    momentum_total += data.momentum();
  }
  /* allows to check energy conservation */
  log.info() << "IC total 4-momentum [GeV]: "
             << momentum_total;  // FIXME: for an info message the "IC" string
                                 // is too cryptic.
  return start_time_;
}

/* evolve - the core of the box, stepping forward in time */
int BoxModus::sanity_check(Particles *particles) {
  int wraps = 0;
  /* fixup positions on startup, particles need to be *inside* the box */
  for (ParticleData &data : particles->data()) {
    FourVector p = data.position();
    if (enforce_periodic_boundaries(p.begin() + 1, p.end(), length_)) {
      ++wraps;
    }
    data.set_4position(p);
  }
  return wraps;
}


/* propagate all particles */
void BoxModus::propagate(Particles *particles,
                         const ExperimentParameters &parameters,
                         const OutputsList &output_list) {
  const auto &log = logger<LogArea::Box>();
  FourVector distance, position;
  for (ParticleData &data : particles->data()) {
    /* propagation for this time step */
    distance = FourVector(0.0,
                          data.velocity() * parameters.timestep_duration());
    log.debug() << data << " motion: " << distance;
    /* treat the box boundaries */
    position = data.position();
    position += distance;
    position.set_x0(parameters.new_particle_time());
    bool wall_hit = enforce_periodic_boundaries(position.begin() + 1,
                                                position.end(), length_);
    const ParticleList incoming_particle{1, data};
    data.set_4position(position);
    if (wall_hit) {
      for (const auto &output : output_list) {
        output->at_interaction(incoming_particle, {1, data});
      }
    }
    log.debug() << data;
  }
}

}  // namespace Smash
