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
#include <sstream>

#include "include/angles.h"
#include "include/constants.h"
#include "include/logging.h"
#include "include/pauliblocking.h"
#include "include/processbranch.h"
#include "include/random.h"
#include "include/resonances.h"
#include "include/width.h"


namespace Smash {

Action::Action(const ParticleList &in_part, float time_of_execution)
    : incoming_particles_(in_part), time_of_execution_(time_of_execution),
      total_weight_(0.) {}

Action::Action(Action &&a)
    : incoming_particles_(std::move(a.incoming_particles_)),
      outgoing_particles_(std::move(a.outgoing_particles_)),
      subprocesses_(std::move(a.subprocesses_)),
      time_of_execution_(a.time_of_execution_),
      total_weight_(a.total_weight_) {}

Action::~Action() = default;

float Action::weight() const {
  return total_weight_;
}

void Action::add_process(ProcessBranch *p) {
  if (p->weight() > really_small) {
    total_weight_ += p->weight();
    subprocesses_.emplace_back(std::move(p));
  }
}

void Action::add_processes(ProcessBranchList pv) {
  if (subprocesses_.empty()) {
    subprocesses_ = std::move(pv);
    for (auto &proc : subprocesses_) {
      total_weight_ += proc->weight();
    }
  } else {
    subprocesses_.reserve(subprocesses_.size() + pv.size());
    for (auto &proc : pv) {
      total_weight_ += proc->weight();
      subprocesses_.emplace_back(std::move(proc));
    }
  }
}

bool Action::is_valid(const Particles &particles) const {
  for (const auto &part : incoming_particles_) {
    // Check if the particles still exists. If it decayed or scattered
    // inelastically it is gone.
    if (!particles.has_data(part.id())) {
      return false;
    }
    // If the particle has scattered elastically, its id_process has changed and
    // we consider it invalid.
    if (particles.data(part.id()).id_process() != part.id_process()) {
      return false;
    }
  }
  return true;
}

bool Action::is_pauliblocked(const Particles & particles,
                             const PauliBlocker* p_bl) const {
  for (const auto &p : outgoing_particles_) {
    if (p.is_baryon() &&
        p_bl->phasespace_dens(p.position().threevec(),
                              p.momentum().threevec(),
                              particles, p.pdgcode()) >
              Random::uniform(0.f,1.f)) {
      return true;
    }
  }
  return false;
}

ParticleList Action::incoming_particles() const {
  ParticleList l;
  for (const auto &part : incoming_particles_) {
    l.emplace_back(part);
  }
  return std::move(l);
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


const ProcessBranch* Action::choose_channel() {
  const auto &log = logger<LogArea::Action>();
  float random_weight = Random::uniform(0.f,total_weight_);
  float weight_sum = 0.;
  /* Loop through all subprocesses and select one by Monte Carlo, based on
   * their weights.  */
  for (const auto &proc : subprocesses_) {
    /* All processes apart from strings should have a well-defined final state. */
    if (proc->particle_number() < 1
        && proc->get_type() != ProcessBranch::String) {
      continue;
    }
    weight_sum += proc->weight();
    if (random_weight <= weight_sum) {
      /* Return the full process information. */
       return proc.get();
    }
  }
  /* Should never get here. */
  log.fatal(source_location, "Problem in choose_channel: ",
            subprocesses_.size(), " ", weight_sum, " ", total_weight_, " ",
            random_weight, "\n", *this);
  throw std::runtime_error("problem in choose_channel");
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

  double momentum_radial = pCM (cms_energy, mass_a, mass_b);
  if (!(momentum_radial > 0.0)) {
    log.warn("Particle: ", t_a.pdgcode(), " radial momentum: ", momentum_radial);
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
  bool violation = false;

  /* Check momentum conservation */
  FourVector momentum_difference;
  for (const auto &part : incoming_particles_) {
    momentum_difference += part.momentum();
  }
  for (const auto &p : outgoing_particles_) {
    momentum_difference -= p.momentum();
  }

  if (fabs(momentum_difference.x0()) > really_small) {
    violation = true;
    log.error("E conservation violation ", momentum_difference.x0());
  }
  if (fabs(momentum_difference.x1()) > really_small) {
    violation = true;
    log.error("px conservation violation ", momentum_difference.x1());
  }
  if (fabs(momentum_difference.x2()) > really_small) {
    violation = true;
    log.error("py conservation violation ", momentum_difference.x2());
  }
  if (fabs(momentum_difference.x3()) > really_small) {
    violation = true;
    log.error("pz conservation violation ", momentum_difference.x3());
  }

  // TODO: check other conservation laws (charge, baryon number, etc)

  if (violation) {
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
