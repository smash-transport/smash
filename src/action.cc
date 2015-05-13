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
#include "include/resonances.h"

namespace Smash {

Action::Action(const ParticleList &in_part, float time_of_execution)
    : incoming_particles_(in_part), time_of_execution_(time_of_execution) {}

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

FourVector Action::get_interaction_point() {
  // Estimate for the interaction point in the calculational frame
  FourVector interaction_point = FourVector(0., 0., 0., 0.);
  for (const auto &part : incoming_particles_) {
    interaction_point += part.position();
  }
  interaction_point /= incoming_particles_.size();
  return interaction_point;
}

void Action::perform(Particles *particles, size_t &id_process) {
  const auto &log = logger<LogArea::Action>();

  for (ParticleData &p : outgoing_particles_) {
    p.set_id_process(id_process);  // store the process id
  }
  if (process_type_ == ProcessType::Elastic) {
    for (std::size_t i = 0; i < incoming_particles_.size(); ++i) {
      outgoing_particles_[i] =
          particles->update(incoming_particles_[i], outgoing_particles_[i]);
    }
  } else {
    particles->replace(incoming_particles_, outgoing_particles_);
  }

  log.debug("Particle map now has ", particles->size(), " elements.");

  check_conservation(id_process);
  id_process++;
}


void Action::sample_cms_momenta() {
  const auto &log = logger<LogArea::Action>();
  /* This function only operates on 2-particle final states. */
  assert(outgoing_particles_.size() == 2);

  ParticleData *p_a = &outgoing_particles_[0];
  ParticleData *p_b = &outgoing_particles_[1];

  const ParticleType &t_a = p_a->type();
  const ParticleType &t_b = p_b->type();

  double mass_a = t_a.mass();
  double mass_b = t_b.mass();

  const double cms_energy = sqrt_s();

  if (cms_energy < t_a.minimum_mass() + t_b.minimum_mass()) {
    throw InvalidResonanceFormation("resonance_formation: not enough energy! " +
      std::to_string(cms_energy) + " " + std::to_string(t_a.minimum_mass()) +
      " " + std::to_string(t_b.minimum_mass()) + " " +
      p_a->pdgcode().string() + " " + p_b->pdgcode().string());
  }

  /* If one of the particles is a resonance, sample its mass. */
  /* TODO: Other particle assumed stable! */
  if (!t_a.is_stable()) {
    mass_a = sample_resonance_mass(t_a, t_b, cms_energy);
  } else if (!t_b.is_stable()) {
    mass_b = sample_resonance_mass(t_b, t_a, cms_energy);
  }

  double momentum_radial = pCM(cms_energy, mass_a, mass_b);
  if (!(momentum_radial > 0.0)) {
    log.warn("Particle: ", t_a.pdgcode(),
             " radial momentum: ", momentum_radial);
    log.warn("Etot: ", cms_energy, " m_a: ", mass_a, " m_b: ", mass_b);
  }
  /* TODO : Angles should be sampled from differential cross section
   * of this process. */
  Angles phitheta;
  phitheta.distribute_isotropically();

  p_a->set_4momentum(mass_a,  phitheta.threevec() * momentum_radial);
  p_b->set_4momentum(mass_b, -phitheta.threevec() * momentum_radial);

  log.debug("p_a: ", *p_a, "\np_b: ", *p_b);
}


void Action::check_conservation(const size_t &id_process) const {
  const auto &log = logger<LogArea::Action>();
  QuantumNumbers before(incoming_particles_);
  QuantumNumbers after(outgoing_particles_);
  std::string err_msg = before.report_deviations(after);
  if (before != after) {
    log.error() << err_msg;
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
