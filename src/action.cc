/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/action.h"

#include <assert.h>
#include <algorithm>
#include <sstream>

#include "smash/angles.h"
#include "smash/constants.h"
#include "smash/kinematics.h"
#include "smash/logging.h"
#include "smash/pauliblocking.h"
#include "smash/processbranch.h"
#include "smash/quantumnumbers.h"

namespace smash {
/// Destructor
Action::~Action() = default;

bool Action::is_valid(const Particles &particles) const {
  return std::all_of(
      incoming_particles_.begin(), incoming_particles_.end(),
      [&particles](const ParticleData &p) { return particles.is_valid(p); });
}

bool Action::is_pauli_blocked(const Particles &particles,
                              const PauliBlocker &p_bl) const {
  // Wall-crossing actions should never be blocked: currently
  // if the action is blocked, a particle continues to propagate in a straight
  // line. This would simply bring it out of the box.
  if (process_type_ == ProcessType::Wall) {
    return false;
  }
  const auto &log = logger<LogArea::PauliBlocking>();
  for (const auto &p : outgoing_particles_) {
    if (p.is_baryon()) {
      const auto f =
          p_bl.phasespace_dens(p.position().threevec(), p.momentum().threevec(),
                               particles, p.pdgcode(), incoming_particles_);
      if (f > random::uniform(0., 1.)) {
        log.debug("Action ", *this, " is pauli-blocked with f = ", f);
        return true;
      }
    }
  }
  return false;
}

const ParticleList &Action::incoming_particles() const {
  return incoming_particles_;
}

void Action::update_incoming(const Particles &particles) {
  for (auto &p : incoming_particles_) {
    p = particles.lookup(p);
  }
}

FourVector Action::get_interaction_point() const {
  // Estimate for the interaction point in the calculational frame
  FourVector interaction_point = FourVector(0., 0., 0., 0.);
  for (const auto &part : incoming_particles_) {
    interaction_point += part.position();
  }
  interaction_point /= incoming_particles_.size();
  return interaction_point;
}

std::pair<double, double> Action::get_potential_at_interaction_point() const {
  const ThreeVector r = get_interaction_point().threevec();
  double UB = 0.;
  double UI3 = 0.;
  /* Check:
   * Lattice is turned on. */
  if (UB_lat_pointer != nullptr) {
    /** \todo TODO(fengli):
     * A Lorentz transformation from the local rest frame to the
     * center of mass frame of the incoming particles is missing here. Since all
     * the actions take place in the center of mass frame of the incoming
     * particles, particles should see potentials different from UB_lat_ or
     * UI3_lat_ which are obtained in the local rest frame. But I don't think
     * the Lorentz transformation is important in the low energy heavy-ion
     * collisions, and turning on potentials violates the Lorentz covariance in
     * the current SMASH version anyway, so I'd like to leave it to another
     * issue in the future.
     */
    UB_lat_pointer->value_at(r, UB);
  }
  if (UI3_lat_pointer != nullptr) {
    UI3_lat_pointer->value_at(r, UI3);
  }
  return std::make_pair(UB, UI3);
}

void Action::input_potential(RectangularLattice<double> *UB_lat,
                             RectangularLattice<double> *UI3_lat,
                             Potentials *pot) {
  UB_lat_pointer = UB_lat;
  UI3_lat_pointer = UI3_lat;
  pot_pointer = pot;
}

void Action::perform(Particles *particles, uint32_t id_process) {
  assert(id_process != 0);
  const auto &log = logger<LogArea::Action>();

  for (ParticleData &p : outgoing_particles_) {
    // store the history info
    if (process_type_ != ProcessType::Wall) {
      p.set_history(p.get_history().collisions_per_particle + 1, id_process,
                    process_type_, time_of_execution_, incoming_particles_);
    }
  }

  /* For elastic collisions and box wall crossings it is not necessary to remove
   * particles from the list and insert new ones, it is enough to update their
   * properties. */
  particles->update(incoming_particles_, outgoing_particles_,
                    (process_type_ != ProcessType::Elastic) &&
                        (process_type_ != ProcessType::Wall));

  log.debug("Particle map now has ", particles->size(), " elements.");

  /* Check the conservation laws if the modifications of the total kinetic
   * energy of the outgoing particles by the mean field potentials are not
   * taken into account. */
  if (UB_lat_pointer == nullptr && UI3_lat_pointer == nullptr) {
    check_conservation(id_process);
  }
}

double Action::kinetic_energy_cms() const {
  const auto potentials = get_potential_at_interaction_point();
  return kinetic_energy_cms<ParticleList>(potentials, outgoing_particles_);
}

std::pair<double, double> Action::sample_masses() const {
  const ParticleType &t_a = outgoing_particles_[0].type();
  const ParticleType &t_b = outgoing_particles_[1].type();
  // start with pole masses
  std::pair<double, double> masses = {t_a.mass(), t_b.mass()};

  const double cms_kin_energy = kinetic_energy_cms();

  if (cms_kin_energy < t_a.min_mass_kinematic() + t_b.min_mass_kinematic()) {
    const std::string reaction = incoming_particles_[0].type().name() +
                                 incoming_particles_[1].type().name() + "→" +
                                 t_a.name() + t_b.name();
    throw InvalidResonanceFormation(
        reaction + ": not enough energy, " + std::to_string(cms_kin_energy) +
        " < " + std::to_string(t_a.min_mass_kinematic()) + " + " +
        std::to_string(t_b.min_mass_kinematic()));
  }

  /* If one of the particles is a resonance, sample its mass. */
  if (!t_a.is_stable() && t_b.is_stable()) {
    masses.first = t_a.sample_resonance_mass(t_b.mass(), cms_kin_energy);
  } else if (!t_b.is_stable() && t_a.is_stable()) {
    masses.second = t_b.sample_resonance_mass(t_a.mass(), cms_kin_energy);
  } else if (!t_a.is_stable() && !t_b.is_stable()) {
    // two resonances in final state
    masses = t_a.sample_resonance_masses(t_b, cms_kin_energy);
  }
  return masses;
}

void Action::sample_angles(std::pair<double, double> masses) {
  const auto &log = logger<LogArea::Action>();

  ParticleData *p_a = &outgoing_particles_[0];
  ParticleData *p_b = &outgoing_particles_[1];

  const double cms_kin_energy = kinetic_energy_cms();

  const double pcm = pCM(cms_kin_energy, masses.first, masses.second);
  if (!(pcm > 0.0)) {
    log.warn("Particle: ", p_a->pdgcode(), " radial momentum: ", pcm);
    log.warn("Ektot: ", cms_kin_energy, " m_a: ", masses.first,
             " m_b: ", masses.second);
  }
  /* Here we assume an isotropic angular distribution. */
  Angles phitheta;
  phitheta.distribute_isotropically();

  p_a->set_4momentum(masses.first, phitheta.threevec() * pcm);
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
    for (const auto &p : incoming_particles_) {
      particle_names << p.type().name();
    }
    particle_names << " vs. ";
    for (const auto &p : outgoing_particles_) {
      particle_names << p.type().name();
    }
    particle_names << "\n";
    const auto &log = logger<LogArea::Action>();
    std::string err_msg = before.report_deviations(after);
    log.error() << particle_names.str() << err_msg;
    /* Pythia does not conserve energy and momentum at high energy, so we just
     * print the error and continue. */
    if ((is_string_soft_process(process_type_)) ||
        (process_type_ == ProcessType::StringHard)) {
      return;
    }
    if (id_process == ID_PROCESS_PHOTON) {
      abort();
      throw std::runtime_error("Conservation laws violated in photon process");
    } else {
      throw std::runtime_error("Conservation laws violated in process " +
                               std::to_string(id_process));
    }
  }
}

std::ostream &operator<<(std::ostream &out, const ActionList &actions) {
  out << "ActionList {\n";
  for (const auto &a : actions) {
    out << "- " << a << '\n';
  }
  return out << '}';
}

}  // namespace smash
