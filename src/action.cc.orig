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
#include "include/outputroutines.h"
#include "include/random.h"
#include "include/resonances.h"


namespace Smash {

Action::Action(const ParticleList &in_part, float time_of_execution)
    : incoming_particles_(in_part), time_of_execution_(time_of_execution),
      total_weight_(0.) {}

Action::~Action() {}

float Action::weight() const {
  return total_weight_;
}

void Action::add_process(const ProcessBranch &p) {
  subprocesses_.push_back(p);
  total_weight_ += p.weight();
}

void Action::add_processes(const ProcessBranchList &pv) {
  for (const auto &proc : pv) {
    subprocesses_.push_back(proc);
    total_weight_ += proc.weight();
  }
}

bool Action::is_valid(const Particles &particles) const {
  for (const auto &part : incoming_particles_) {
    /* Check if the particles still exist. */
    if (!particles.has_data(part.id())) {
      return false;
    }
    /* Check if particles have scattered in the meantime
     * (by checking if their energy or momentum has changed). */
    if ((fabs(part.momentum().x0()
             - particles.data(part.id()).momentum().x0()) > really_small)
        || (fabs(part.momentum().x1()
                 - particles.data(part.id()).momentum().x1()) > really_small)
        || (fabs(part.momentum().x2()
                 - particles.data(part.id()).momentum().x2()) > really_small)
        || (fabs(part.momentum().x3()
                 - particles.data(part.id()).momentum().x3()) > really_small)) {
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
  if (!(momentum_radial > 0.0))
    printf("Warning: radial momenta %g \n", momentum_radial);
  /* XXX: Angles should be sampled from differential cross section
   * of this process. */
  Angles phitheta;
  phitheta.distribute_isotropically();
  if (!(energy1 > mass1)) {
    printf("Particle %s radial momenta %g phi %g cos_theta %g\n",
           t1.pdgcode().string().c_str(), momentum_radial,
           phitheta.phi(), phitheta.costheta());
    printf("Etot: %g m_a: %g m_b %g E_a: %g\n", cms_energy, mass1, mass2,
           energy1);
  }

  p1->set_momentum(FourVector(energy1, phitheta.threevec() * momentum_radial));
  p2->set_momentum(FourVector(cms_energy - energy1,
                              -phitheta.threevec() * momentum_radial));

  printd("p0: %g %g \n", p1->momentum().x0(), p2->momentum().x0());
  printd("p1: %g %g \n", p1->momentum().x1(), p2->momentum().x1());
  printd("p2: %g %g \n", p1->momentum().x2(), p2->momentum().x2());
  printd("p3: %g %g \n", p1->momentum().x3(), p2->momentum().x3());
}


void Action::check_conservation(const size_t &id_process) const {
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
    printf("Process %zu\n", id_process);
    printf("Warning: E conservation violation %g\n",
           momentum_difference.x0());
  }
  if (fabs(momentum_difference.x1()) > really_small)
    printf("Warning: px conservation violation %g\n",
           momentum_difference.x1());
  if (fabs(momentum_difference.x2()) > really_small)
    printf("Warning: py conservation violation %g\n",
           momentum_difference.x2());
  if (fabs(momentum_difference.x3()) > really_small)
    printf("Warning: pz conservation violation %g\n",
           momentum_difference.x3());

  // TODO: check other conservation laws (baryon number etc)
}

}  // namespace Smash
