/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/action.h"

#include <assert.h>
#include <algorithm>
#include <sstream>

#include "include/angles.h"
#include "include/constants.h"
#include "include/kinematics.h"
#include "include/logging.h"
#include "include/pauliblocking.h"
#include "include/processbranch.h"
#include "include/quantumnumbers.h"

namespace Smash {

Action::Action(const ParticleList &in_part, float time)
              : incoming_particles_(in_part),
                time_of_execution_(time+in_part[0].position().x0()) {}

Action::~Action() = default;


bool Action::is_valid(const Particles &particles) const {
  return std::all_of(
      incoming_particles_.begin(), incoming_particles_.end(),
      [&particles](const ParticleData &p) { return particles.is_valid(p); });
}

bool Action::is_pauli_blocked(const Particles & particles,
                             const PauliBlocker& p_bl) const {
  const auto &log = logger<LogArea::PauliBlocking>();
  for (const auto &p : outgoing_particles_) {
    if (p.is_baryon()) {
      const auto f = p_bl.phasespace_dens(p.position().threevec(),
                                           p.momentum().threevec(),
                                           particles, p.pdgcode());
      if (f >  Random::uniform(0.f, 1.f)) {
        log.debug("Action ", *this, " is pauli-blocked with f = ", f);
        return true;
      }
    }
  }
  return false;
}

const ParticleList& Action::incoming_particles() const {
  return incoming_particles_;
}

void Action::update_incoming(const Particles &particles) {
  for (auto &p : incoming_particles_) {
    p = particles.lookup(p);
  }
}

FourVector Action::get_interaction_point() {
  // Estimate for the interaction point in the calculational frame
  FourVector interaction_point = FourVector(0., 0., 0., 0.);
  for (const auto &part : incoming_particles_) {
    interaction_point += part.position();
  }
  interaction_point /= incoming_particles_.size();
  return interaction_point;
}

void Action::perform(Particles *particles, uint32_t id_process) {
  assert(id_process != 0);
  const auto &log = logger<LogArea::Action>();

  for (ParticleData &p : outgoing_particles_) {
    // store the history info
    if (process_type_ != ProcessType::Wall) {
      p.set_history(p.get_history().collisions_per_particle+1, id_process,
                  process_type_, time_of_execution_, incoming_particles_);
    }
  }

  particles->update(incoming_particles_, outgoing_particles_,
                    process_type_ != ProcessType::Elastic);

  log.debug("Particle map now has ", particles->size(), " elements.");

  check_conservation(id_process);
}


std::pair<double, double> Action::sample_masses() const {
  const ParticleType &t_a = outgoing_particles_[0].type();
  const ParticleType &t_b = outgoing_particles_[1].type();

  // start with pole masses
  std::pair<double, double> masses = {t_a.mass(), t_b.mass()};

  const double cms_energy = sqrt_s();

  if (cms_energy < t_a.minimum_mass() + t_b.minimum_mass()) {
    const std::string reaction = incoming_particles_[0].type().name() +
                                 incoming_particles_[1].type().name() + "→" +
                                 t_a.name() + t_b.name();
    throw InvalidResonanceFormation(reaction + ": not enough energy, " +
      std::to_string(cms_energy) + " < " +
      std::to_string(t_a.minimum_mass()) + " + " +
      std::to_string(t_b.minimum_mass()));
  }

  /* If one of the particles is a resonance, sample its mass. */
  if (!t_a.is_stable() && t_b.is_stable()) {
    masses.first = t_a.sample_resonance_mass(t_b.mass(), cms_energy);
  } else if (!t_b.is_stable() && t_a.is_stable()) {
    masses.second = t_b.sample_resonance_mass(t_a.mass(), cms_energy);
  } else if (!t_a.is_stable() && !t_b.is_stable()) {
    // two resonances in final state
    masses = t_a.sample_resonance_masses(t_b, cms_energy);
  }

  return masses;
}


void Action::sample_angles(std::pair<double, double> masses) {
  const auto &log = logger<LogArea::Action>();

  ParticleData *p_a = &outgoing_particles_[0];
  ParticleData *p_b = &outgoing_particles_[1];

  const double cms_energy = sqrt_s();

  const double pcm = pCM(cms_energy, masses.first, masses.second);
  if (!(pcm > 0.0)) {
    log.warn("Particle: ", p_a->pdgcode(), " radial momentum: ", pcm);
    log.warn("Etot: ", cms_energy, " m_a: ", masses.first,
                                   " m_b: ", masses.second);
  }
  /* Here we assume an isotropic angular distribution. */
  Angles phitheta;
  phitheta.distribute_isotropically();

  p_a->set_4momentum(masses.first,   phitheta.threevec() * pcm);
  p_b->set_4momentum(masses.second, -phitheta.threevec() * pcm);

  log.debug("p_a: ", *p_a, "\np_b: ", *p_b);
}


void Action::sample_2body_phasespace() {
  /* This function only operates on 2-particle final states. */
  assert(outgoing_particles_.size() == 2);
  // first sample the masses
  const std::pair<double, double> masses = sample_masses();
  // after the masses are fixed (and thus also pcm), sample the angles
  sample_angles(masses);
}


void Action::check_conservation(const uint32_t id_process) const {
  QuantumNumbers before(incoming_particles_);
  QuantumNumbers after(outgoing_particles_);
  if (before != after) {
    std::stringstream particle_names;
    for (const auto& p : incoming_particles_) {
      particle_names << p.type().name();
    }
    particle_names << " vs. ";
    for (const auto& p : outgoing_particles_) {
      particle_names << p.type().name();
    }
    particle_names << "\n";
    const auto &log = logger<LogArea::Action>();
    std::string err_msg = before.report_deviations(after);
    log.error() << particle_names.str() << err_msg;
    throw std::runtime_error("Conservation laws violated in process " +
                             std::to_string(id_process));
  }
}

std::ostream &operator<<(std::ostream &out, const ActionList &actions) {
  out << "ActionList {\n";
  for (const auto &a : actions) {
    out << "- " << a << '\n';
  }
  return out << '}';
}

}  // namespace Smash
