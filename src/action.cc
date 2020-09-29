/*
 *
 *    Copyright (c) 2014-2019
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
#include "smash/potential_globals.h"
#include "smash/processbranch.h"
#include "smash/quantumnumbers.h"

namespace smash {
/// Destructor
Action::~Action() = default;
static constexpr int LPauliBlocking = LogArea::PauliBlocking::id;

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
  for (const auto &p : outgoing_particles_) {
    if (p.is_baryon()) {
      const auto f =
          p_bl.phasespace_dens(p.position().threevec(), p.momentum().threevec(),
                               particles, p.pdgcode(), incoming_particles_);
      if (f > random::uniform(0., 1.)) {
        logg[LPauliBlocking].debug("Action ", *this,
                                   " is pauli-blocked with f = ", f);
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
  /*
   * In case of periodic boundaries interaction point is not necessarily
   * (x1 + x2)/2. Consider only one dimension, e.g. x, the rest are analogous.
   * Instead of x, there can be x + k * L, where k is any integer and L
   * is period.Interaction point is either. Therefore, interaction point is
   * (x1 + k * L + x2 + m * L) / 2  = (x1 + x2) / 2 + n * L / 2. We need
   * this interaction point to be with [0, L], so n can be {-1, 0, 1}.
   * Which n to choose? Our guiding principle is that n should be such that
   * interaction point is closest to interacting particles.
   */
  if (box_length_ > 0) {
    assert(incoming_particles_.size() == 2);
    const FourVector r1 = incoming_particles_[0].position(),
                     r2 = incoming_particles_[1].position(), r = r1 - r2;
    for (int i = 1; i < 4; i++) {
      const double d = std::abs(r[i]);
      if (d > 0.5 * box_length_) {
        if (interaction_point[i] >= 0.5 * box_length_) {
          interaction_point[i] -= 0.5 * box_length_;
        } else {
          interaction_point[i] += 0.5 * box_length_;
        }
      }
    }
  }
  return interaction_point;
}

std::pair<FourVector, FourVector> Action::get_potential_at_interaction_point()
    const {
  const ThreeVector r = get_interaction_point().threevec();
  FourVector UB = FourVector();
  FourVector UI3 = FourVector();
  /* Check:
   * Lattice is turned on. */
  if (UB_lat_pointer != nullptr) {
    UB_lat_pointer->value_at(r, UB);
  }
  if (UI3_lat_pointer != nullptr) {
    UI3_lat_pointer->value_at(r, UI3);
  }
  return std::make_pair(UB, UI3);
}

void Action::perform(Particles *particles, uint32_t id_process) {
  assert(id_process != 0);

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

  logg[LAction].debug("Particle map now has ", particles->size(), " elements.");

  /* Check the conservation laws if the modifications of the total kinetic
   * energy of the outgoing particles by the mean field potentials are not
   * taken into account. */
  if (UB_lat_pointer == nullptr && UI3_lat_pointer == nullptr) {
    check_conservation(id_process);
  }
}

FourVector Action::total_momentum_of_outgoing_particles() const {
  const auto potentials = get_potential_at_interaction_point();
  /* scale_B returns the difference of the total force scales of the skyrme
   * potential between the initial and final states. */
  double scale_B = 0.0;
  /* scale_I3 returns the difference of the total force scales of the symmetry
   * potential between the initial and final states. */
  double scale_I3 = 0.0;
  for (const auto &p_in : incoming_particles_) {
    // Get the force scale of the incoming particle.
    const auto scale =
        ((pot_pointer != nullptr) ? pot_pointer->force_scale(p_in.type())
                                  : std::make_pair(0.0, 0));
    scale_B += scale.first;
    scale_I3 += scale.second * p_in.type().isospin3_rel();
  }
  for (const auto &p_out : outgoing_particles_) {
    // Get the force scale of the outgoing particle.
    const auto scale = ((pot_pointer != nullptr)
                            ? pot_pointer->force_scale(type_of_pout(p_out))
                            : std::make_pair(0.0, 0));
    scale_B -= scale.first;
    scale_I3 -= scale.second * type_of_pout(p_out).isospin3_rel();
  }
  /* Rescale to get the potential difference between the
   * initial and final state, and thus get the total momentum
   * of the outgoing particles*/
  return total_momentum() + potentials.first * scale_B +
         potentials.second * scale_I3;
}

void Action::assign_formation_time_to_outgoing_particles() {
  /* Find incoming particle with largest formation time i.e. the last formed
   * incoming particle. */
  auto last_formed_in_part =
      std::max_element(incoming_particles_.begin(), incoming_particles_.end(),
                       [](const ParticleData &a, const ParticleData &b) {
                         return a.formation_time() < b.formation_time();
                       });

  const double form_time_begin = last_formed_in_part->begin_formation_time();
  const double sc = last_formed_in_part->initial_xsec_scaling_factor();

  if (last_formed_in_part->formation_time() > time_of_execution_) {
    for (ParticleData &new_particle : outgoing_particles_) {
      new_particle.set_slow_formation_times(
          form_time_begin, last_formed_in_part->formation_time());
      new_particle.set_cross_section_scaling_factor(sc);
    }
  } else {
    for (ParticleData &new_particle : outgoing_particles_) {
      new_particle.set_formation_time(time_of_execution_);
    }
  }
}

std::pair<double, double> Action::sample_masses(
    const double kinetic_energy_cm) const {
  const ParticleType &t_a = outgoing_particles_[0].type();
  const ParticleType &t_b = outgoing_particles_[1].type();
  // start with pole masses
  std::pair<double, double> masses = {t_a.mass(), t_b.mass()};

  if (kinetic_energy_cm < t_a.min_mass_kinematic() + t_b.min_mass_kinematic()) {
    const std::string reaction = incoming_particles_[0].type().name() +
                                 incoming_particles_[1].type().name() + "→" +
                                 t_a.name() + t_b.name();
    throw InvalidResonanceFormation(
        reaction + ": not enough energy, " + std::to_string(kinetic_energy_cm) +
        " < " + std::to_string(t_a.min_mass_kinematic()) + " + " +
        std::to_string(t_b.min_mass_kinematic()));
  }

  /* If one of the particles is a resonance, sample its mass. */
  if (!t_a.is_stable() && t_b.is_stable()) {
    masses.first = t_a.sample_resonance_mass(t_b.mass(), kinetic_energy_cm);
  } else if (!t_b.is_stable() && t_a.is_stable()) {
    masses.second = t_b.sample_resonance_mass(t_a.mass(), kinetic_energy_cm);
  } else if (!t_a.is_stable() && !t_b.is_stable()) {
    // two resonances in final state
    masses = t_a.sample_resonance_masses(t_b, kinetic_energy_cm);
  }
  return masses;
}

void Action::sample_angles(std::pair<double, double> masses,
                           const double kinetic_energy_cm) {
  ParticleData *p_a = &outgoing_particles_[0];
  ParticleData *p_b = &outgoing_particles_[1];

  const double pcm = pCM(kinetic_energy_cm, masses.first, masses.second);
  if (!(pcm > 0.0)) {
    logg[LAction].warn("Particle: ", p_a->pdgcode(), " radial momentum: ", pcm);
    logg[LAction].warn("Ektot: ", kinetic_energy_cm, " m_a: ", masses.first,
                       " m_b: ", masses.second);
  }
  /* Here we assume an isotropic angular distribution. */
  Angles phitheta;
  phitheta.distribute_isotropically();

  p_a->set_4momentum(masses.first, phitheta.threevec() * pcm);
  p_b->set_4momentum(masses.second, -phitheta.threevec() * pcm);
  /* Debug message is printed before boost, so that p_a and p_b are
   * the momenta in the center of mass frame and thus opposite to
   * each other.*/
  logg[LAction].debug("p_a: ", *p_a, "\np_b: ", *p_b);
}

void Action::sample_2body_phasespace() {
  /* This function only operates on 2-particle final states. */
  assert(outgoing_particles_.size() == 2);
  const FourVector p_tot = total_momentum_of_outgoing_particles();
  const double cm_kin_energy = p_tot.abs();
  // first sample the masses
  const std::pair<double, double> masses = sample_masses(cm_kin_energy);
  // after the masses are fixed (and thus also pcm), sample the angles
  sample_angles(masses, cm_kin_energy);
}

void Action::sample_3body_phasespace() {
  assert(outgoing_particles_.size() == 3);
  if (!outgoing_particles_[0].type().is_stable() ||
      !outgoing_particles_[1].type().is_stable() ||
      !outgoing_particles_[2].type().is_stable()) {
    throw std::invalid_argument(
        "sample_3body_phasespace: Found resonance in to be sampled outgoing "
        "particles, but assumes stable particles.");
  }

  const double m_a = outgoing_particles_[0].type().mass(),
               m_b = outgoing_particles_[1].type().mass(),
               m_c = outgoing_particles_[2].type().mass();
  const double sqrts = sqrt_s();

  // sample mab from pCM(sqrt, mab, mc) pCM (mab, ma, mb) <= sqrts^2/4
  double mab, r, probability, pcm_ab, pcm;
  do {
    mab = random::uniform(m_a + m_b, sqrts - m_c);
    r = random::canonical();
    pcm = pCM(sqrts, mab, m_c);
    pcm_ab = pCM(mab, m_a, m_b);
    probability = pcm * pcm_ab * 4 / (sqrts * sqrts);
  } while (r > probability);
  Angles phitheta;
  phitheta.distribute_isotropically();
  outgoing_particles_[2].set_4momentum(m_c, pcm * phitheta.threevec());
  const ThreeVector beta_cm =
      pcm * phitheta.threevec() / std::sqrt(pcm * pcm + mab * mab);

  phitheta.distribute_isotropically();
  outgoing_particles_[0].set_4momentum(m_a, pcm_ab * phitheta.threevec());
  outgoing_particles_[1].set_4momentum(m_b, -pcm_ab * phitheta.threevec());
  outgoing_particles_[0].boost_momentum(beta_cm);
  outgoing_particles_[1].boost_momentum(beta_cm);
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
    std::string err_msg = before.report_deviations(after);
    logg[LAction].error() << particle_names.str() << err_msg;
    /* Pythia does not conserve energy and momentum at high energy, so we just
     * print the error and continue. */
    if ((is_string_soft_process(process_type_)) ||
        (process_type_ == ProcessType::StringHard)) {
      return;
    }
    if (id_process == ID_PROCESS_PHOTON) {
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
