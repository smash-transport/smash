/*
 *
 *    Copyright (c) 2014
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
#include "include/random.h"
#include "include/resonances.h"


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

void Action::add_process(ProcessBranch p) {
  total_weight_ += p.weight();
  subprocesses_.emplace_back(std::move(p));
}

void Action::add_processes(ProcessBranchList pv) {
  if (subprocesses_.empty()) {
    subprocesses_ = std::move(pv);
    for (auto &proc : subprocesses_) {
      total_weight_ += proc.weight();
    }
  } else {
    subprocesses_.reserve(subprocesses_.size() + pv.size());
    for (auto &proc : pv) {
      total_weight_ += proc.weight();
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

ParticleList Action::incoming_particles() const {
  ParticleList l;
  for (const auto &part : incoming_particles_) {
    l.emplace_back(part);
  }
  return std::move(l);
}


ParticleList Action::choose_channel() {
  if (total_weight_ < really_small) {
    return ParticleList();
  }
  double random_interaction = Random::canonical();
  float interaction_probability = 0.0;
  /* Loop through all subprocesses and select one by Monte Carlo, based on
   * their weights.  */
  for (const auto &proc : subprocesses_) {
    if (proc.pdg_list().size() < 1
        || proc.pdg_list()[0] == PdgCode::invalid()) {
      continue;
    }
    interaction_probability += proc.weight() / total_weight_;
    if (random_interaction < interaction_probability) {
      return proc.particle_list();
    }
  }
  /* Should never get here. */
  std::stringstream ss;
  ss << "problem in choose_channel: " << subprocesses_.size() << " " <<
         interaction_probability << " " << total_weight_;
  throw std::runtime_error(ss.str());
}


void Action::sample_cms_momenta() {
  const auto &log = logger<LogArea::Action>();
  /* This function only operates on 2-particle final states. */
  assert(outgoing_particles_.size() == 2);

  ParticleData *p1 = &outgoing_particles_[0];
  ParticleData *p2 = &outgoing_particles_[1];

  const ParticleType &t1 = p1->type();
  const ParticleType &t2 = p2->type();

  double mass1 = t1.mass();
  double mass2 = t2.mass();

  const double cms_energy = sqrt_s();

  if (cms_energy < t1.minimum_mass() + t2.minimum_mass()) {
    throw InvalidResonanceFormation("resonance_formation: not enough energy! " +
      std::to_string(cms_energy) + " " + std::to_string(t1.minimum_mass()) +
      " " + std::to_string(t2.minimum_mass()) + " " +
      p1->pdgcode().string() + " " + p2->pdgcode().string());
  }

  /* If one of the particles is a resonance, sample its mass. */
  /* XXX: Other particle assumed stable! */
  if (!t1.is_stable()) {
    mass1 = sample_resonance_mass(t1, t2, cms_energy);
  } else if (!t2.is_stable()) {
    mass2 = sample_resonance_mass(t2, t1, cms_energy);
  }

  double energy1 = (cms_energy * cms_energy + mass1 * mass1 - mass2 * mass2) /
                   (2.0 * cms_energy);
  double momentum_radial = std::sqrt(energy1 * energy1 - mass1 * mass1);
  if (!(momentum_radial > 0.0)) {
    log.warn("radial momenta ", momentum_radial);
  }
  /* XXX: Angles should be sampled from differential cross section
   * of this process. */
  Angles phitheta;
  phitheta.distribute_isotropically();
  if (!(energy1 > mass1)) {
    log.info("Particle ", t1.pdgcode(), " radial momenta ", momentum_radial,
             phitheta);
    log.info("Etot: ", cms_energy, " m_a: ", mass1, " m_b: ", mass2, " E_a: ",
             energy1);
  }

  p1->set_4momentum(FourVector(energy1, phitheta.threevec() * momentum_radial));
  p2->set_4momentum(FourVector(cms_energy - energy1,
                               -phitheta.threevec() * momentum_radial));

  log.debug("p1: ", *p1, "\np2: ", *p2);
}


void Action::check_conservation(const size_t &id_process) const {
  const auto &log = logger<LogArea::Action>();
  /* Check momentum conservation */
  FourVector momentum_difference;
  for (const auto &part : incoming_particles_) {
    momentum_difference += part.momentum();
  }
  for (const auto &p : outgoing_particles_) {
    momentum_difference -= p.momentum();
  }

  /* TODO: throw an exception */
  if (fabs(momentum_difference.x0()) > really_small) {
    log.warn("Process ", id_process, "\nE conservation violation ",
             momentum_difference.x0());
  }
  if (fabs(momentum_difference.x1()) > really_small) {
    log.warn("px conservation violation ", momentum_difference.x1());
  }
  if (fabs(momentum_difference.x2()) > really_small) {
    log.warn("py conservation violation ", momentum_difference.x2());
  }
  if (fabs(momentum_difference.x3()) > really_small) {
    log.warn("pz conservation violation ", momentum_difference.x3());
  }

  // TODO: check other conservation laws (baryon number etc)
}

std::ostream &operator<<(std::ostream &out, const ActionList &actions) {
  out << "ActionList {\n";
  for (const auto &a : actions) {
    out << "- " << a << '\n';
  }
  return out << '}';
}

}  // namespace Smash
