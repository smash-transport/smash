/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/decayaction.h"

#include "include/angles.h"
#include "include/decaymodes.h"
#include "include/kinematics.h"
#include "include/logging.h"
#include "include/pdgcode.h"
#include "include/resonances.h"

namespace Smash {

DecayAction::DecayAction(const ParticleData &p, float time)
    : Action({p}, time),
      total_width_(0.) {}

void DecayAction::add_decays(DecayBranchList pv) {
  add_processes<DecayBranch>(std::move(pv), decay_channels_, total_width_);
}

void DecayAction::add_decay(DecayBranchPtr p) {
  add_process<DecayBranch>(p, decay_channels_, total_width_);
}

double DecayAction::sqrt_s() const {
  return incoming_particles_[0].momentum().abs();
}


void DecayAction::one_to_two() {
  /* Sample the masses and momenta. */
  sample_cms_momenta();
}

void DecayAction::one_to_three() {
  const auto &log = logger<LogArea::DecayModes>();
  ParticleData &outgoing_a = outgoing_particles_[0];
  ParticleData &outgoing_b = outgoing_particles_[1];
  ParticleData &outgoing_c = outgoing_particles_[2];
  const ParticleType &outgoing_a_type = outgoing_a.type();
  const ParticleType &outgoing_b_type = outgoing_b.type();
  const ParticleType &outgoing_c_type = outgoing_c.type();

  log.debug("Note: Doing 1->3 decay!");

  const double mass_a = outgoing_a_type.mass();
  const double mass_b = outgoing_b_type.mass();
  const double mass_c = outgoing_c_type.mass();
  const double mass_resonance = incoming_particles_[0].effective_mass();

  /* mandelstam-s limits for pairs ab and bc */
  const double s_ab_max = (mass_resonance - mass_c) * (mass_resonance - mass_c);
  const double s_ab_min = (mass_a + mass_b) * (mass_a + mass_b);
  const double s_bc_max = (mass_resonance - mass_a) * (mass_resonance - mass_a);
  const double s_bc_min = (mass_b + mass_c) * (mass_b + mass_c);

  log.debug("s_ab limits: ", s_ab_min, " ", s_ab_max);
  log.debug("s_bc limits: ", s_bc_min, " ", s_bc_max);

  /* randomly pick values for s_ab and s_bc
   * until the pair is within the Dalitz plot */
  double dalitz_bc_max = 0.0, dalitz_bc_min = 1.0;
  double s_ab = 0.0, s_bc = 0.5;
  while (s_bc > dalitz_bc_max || s_bc < dalitz_bc_min) {
    s_ab = Random::uniform(s_ab_min, s_ab_max);
    s_bc = Random::uniform(s_bc_min, s_bc_max);
    const double e_b_rest =
      (s_ab - mass_a * mass_a + mass_b * mass_b) / (2 * std::sqrt(s_ab));
    const double e_c_rest =
      (mass_resonance * mass_resonance - s_ab - mass_c * mass_c) /
      (2 * std::sqrt(s_ab));
    dalitz_bc_max = (e_b_rest + e_c_rest) * (e_b_rest + e_c_rest) -
                    (std::sqrt(e_b_rest * e_b_rest - mass_b * mass_b) -
                     std::sqrt(e_c_rest * e_c_rest - mass_c * mass_c)) *
                     (std::sqrt(e_b_rest * e_b_rest - mass_b * mass_b) -
                      std::sqrt(e_c_rest * e_c_rest - mass_c * mass_c));
    dalitz_bc_min = (e_b_rest + e_c_rest) * (e_b_rest + e_c_rest) -
                    (std::sqrt(e_b_rest * e_b_rest - mass_b * mass_b) +
                     std::sqrt(e_c_rest * e_c_rest - mass_c * mass_c)) *
                     (std::sqrt(e_b_rest * e_b_rest - mass_b * mass_b) +
                      std::sqrt(e_c_rest * e_c_rest - mass_c * mass_c));
  }

  log.debug("s_ab: ", s_ab, " s_bc: ", s_bc, " min: ", dalitz_bc_min, " max: ",
            dalitz_bc_max);

  /* Compute energy and momentum magnitude */
  const double energy_a =
      (mass_resonance * mass_resonance + mass_a * mass_a - s_bc) /
      (2 * mass_resonance);
  const double energy_c =
      (mass_resonance * mass_resonance + mass_c * mass_c - s_ab) /
      (2 * mass_resonance);
  const double energy_b =
      (s_ab + s_bc - mass_a * mass_a - mass_c * mass_c) / (2 * mass_resonance);
  const double momentum_a = std::sqrt(energy_a * energy_a - mass_a * mass_a);
  const double momentum_c = std::sqrt(energy_c * energy_c - mass_c * mass_c);
  const double momentum_b = std::sqrt(energy_b * energy_b - mass_b * mass_b);

  const double total_energy = sqrt_s();
  if (std::abs(energy_a + energy_b + energy_c - total_energy) > really_small) {
    log.warn("1->3: Ea + Eb + Ec: ", energy_a + energy_b + energy_c,
             " Total E: ", total_energy);
  }
  log.debug("Calculating the angles...");

  /* momentum_a direction is random */
  Angles phitheta;
  phitheta.distribute_isotropically();
  /* This is the angle of the plane of the three decay particles */
  outgoing_a.set_4momentum(mass_a, phitheta.threevec() * momentum_a);

  /* Angle between a and b */
  double theta_ab = std::acos(
      (energy_a * energy_b - 0.5 * (s_ab - mass_a * mass_a - mass_b * mass_b)) /
      (momentum_a * momentum_b));
  log.debug("theta_ab: ", theta_ab, " Ea: ", energy_a, " Eb: ", energy_b,
            " sab: ", s_ab, " pa: ", momentum_a, " pb: ", momentum_b);
  bool phi_has_changed = phitheta.add_to_theta(theta_ab);
  outgoing_b.set_4momentum(mass_b, phitheta.threevec() * momentum_b);

  /* Angle between b and c */
  double theta_bc = std::acos(
      (energy_b * energy_c - 0.5 * (s_bc - mass_b * mass_b - mass_c * mass_c)) /
      (momentum_b * momentum_c));
  log.debug("theta_bc: ", theta_bc, " Eb: ", energy_b, " Ec: ", energy_c,
            " sbc: ", s_bc, " pb: ", momentum_b, " pc: ", momentum_c);
  // pass information on whether phi has changed during the last adding
  // on to add_to_theta:
  phitheta.add_to_theta(theta_bc, phi_has_changed);
  outgoing_c.set_4momentum(mass_c, phitheta.threevec() * momentum_c);

  /* Momentum check */
  FourVector ptot = outgoing_a.momentum() + outgoing_b.momentum() +
                    outgoing_c.momentum();

  if (std::abs(ptot.x0() - total_energy) > really_small) {
    log.warn("1->3 energy not conserved! Before: ", total_energy, " After: ",
             ptot.x0());
  }
  if (std::abs(ptot.x1()) > really_small ||
      std::abs(ptot.x2()) > really_small ||
      std::abs(ptot.x3()) > really_small) {
    log.warn("1->3 momentum check failed. Total momentum: ", ptot.threevec());
  }

  log.debug("outgoing_a: ", outgoing_a.momentum(),
            "\noutgoing_b: ", outgoing_b.momentum(),
            "\noutgoing_c: ", outgoing_c.momentum());
}


void DecayAction::generate_final_state() {
  const auto &log = logger<LogArea::DecayModes>();
  log.debug("Process: Resonance decay. ");
  /*
   * Execute a decay process for the selected particle.
   *
   * Randomly select one of the decay modes of the particle
   * according to their relative weights. Then decay the particle
   * by calling function one_to_two or one_to_three.
   */
  const DecayBranch* proc = choose_channel<DecayBranch>(
      decay_channels_, total_width_);
  outgoing_particles_ = proc->particle_list();
  process_type_ = proc->get_type();
  L_ = proc->angular_momentum();

  switch (outgoing_particles_.size()) {
  case 2:
    one_to_two();
    break;
  case 3:
    one_to_three();
    break;
  default:
    throw InvalidDecay(
        "DecayAction::perform: Only 1->2 or 1->3 processes are supported. "
        "Decay from 1->" + std::to_string(outgoing_particles_.size()) +
        " was requested. (PDGcode=" + incoming_particles_[0].pdgcode().string()
        + ", mass=" + std::to_string(incoming_particles_[0].effective_mass())
        + ")");
  }

  /* Set positions and boost back. */
  ThreeVector velocity_CM = incoming_particles_[0].velocity();
  for (auto &p : outgoing_particles_) {
    log.debug("particle momenta in lrf ", p);
    p.boost_momentum(-velocity_CM);
    p.set_4position(incoming_particles_[0].position());
    log.debug("particle momenta in comp ", p);
  }
}


void DecayAction::sample_cms_momenta() {
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
    mass_a = sample_resonance_mass(t_a, t_b.mass(), cms_energy, L_);
  } else if (!t_b.is_stable()) {
    mass_b = sample_resonance_mass(t_b, t_a.mass(), cms_energy, L_);
  }

  double momentum_radial = pCM(cms_energy, mass_a, mass_b);
  if (!(momentum_radial > 0.0)) {
    log.warn("Particle: ", t_a.pdgcode(),
             " radial momentum: ", momentum_radial);
    log.warn("Etot: ", cms_energy, " m_a: ", mass_a, " m_b: ", mass_b);
  }
  /* Here we assume an isotropic angular distribution. */
  Angles phitheta;
  phitheta.distribute_isotropically();

  p_a->set_4momentum(mass_a,  phitheta.threevec() * momentum_radial);
  p_b->set_4momentum(mass_b, -phitheta.threevec() * momentum_radial);

  log.debug("p_a: ", *p_a, "\np_b: ", *p_b);
}


float DecayAction::raw_weight_value() const {
  return total_width_;
}


void DecayAction::format_debug_output(std::ostream &out) const {
  out << "Decay of " << incoming_particles_ << " to " << outgoing_particles_ <<
  ", sqrt(s)=" << format(sqrt_s(), "GeV", 11, 9);
}

}  // namespace Smash
