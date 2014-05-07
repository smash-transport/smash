/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/action.h"

#include "include/angles.h"
#include "include/outputroutines.h"
#include "include/pdgcode.h"
#include "include/resonances.h"

namespace Smash {


DecayAction::DecayAction(const std::vector<int> &in_part,
                         float time_of_execution, int interaction_type)
    : Action (in_part, time_of_execution, interaction_type) {}



/**
 * Kinematics of a 1-to-2 decay process.
 *
 * Given a resonance ID and the PDG codes of decay product particles,
 * sample the momenta and position of the products and add them
 * to the active particles data structure.
 *
 * \param[in,out] particles Particles in the simulation.
 *
 * \return The ID of the first new particle.
 */
int DecayAction::one_to_two (Particles *particles) {
  int resonance_id = incoming_particles_[0];
  PdgCode type_a = outgoing_particles_[0].pdgcode();
  PdgCode type_b = outgoing_particles_[1].pdgcode();

  //TODO(weil) This may be checked already when reading in the possible
  //decay channels.
  if (! type_a.is_hadron() || ! type_b.is_hadron()) {
    printf("Warning: decay products A: %s B: %s\n", type_a.string().c_str()
                                                  , type_b.string().c_str());
  }

  /* Add two new particles */
  ParticleData new_particle_a, new_particle_b;
  new_particle_a.set_pdgcode(type_a);
  new_particle_b.set_pdgcode(type_b);

  double mass_a = particles->particle_type(type_a).mass(),
    mass_b = particles->particle_type(type_b).mass();
  const double total_energy = particles->data(resonance_id).momentum().x0();

  /* If one of the particles is resonance, sample its mass */
  /* XXX: Other particle assumed stable! */
  if (particles->particle_type(type_a).width() > 0) {
    mass_a = sample_resonance_mass(particles, type_a, type_b, total_energy);
  } else if (particles->particle_type(type_b).width() > 0) {
    mass_b = sample_resonance_mass(particles, type_b, type_a, total_energy);
  }

  /* Sample the momenta */
  sample_cms_momenta(&new_particle_a, &new_particle_b, total_energy,
                     mass_a, mass_b);

  /* Both decay products begin from the same point */
  FourVector decay_point = particles->data(resonance_id).position();
  new_particle_a.set_position(decay_point);
  new_particle_b.set_position(decay_point);

  int id_first_new = particles->add_data(new_particle_a);
  particles->add_data(new_particle_b);

  /* 2 new particles created; return the id of the first one */
  return id_first_new;
}

/**
 * Kinematics of a 1-to-3 decay process.
 *
 * Given a resonance ID and the PDG codes of decay product particles,
 * sample the momenta and position of the products and add them
 * to the active particles data structure.
 *
 * \param[in,out] particles Particles in the simulation.
 *
 * \return The ID of the first new particle.
 */
int DecayAction::one_to_three (Particles *particles) {
  int resonance_id = incoming_particles_[0];
  PdgCode type_a = outgoing_particles_[0].pdgcode();
  PdgCode type_b = outgoing_particles_[1].pdgcode();
  PdgCode type_c = outgoing_particles_[2].pdgcode();

  //TODO(weil) This may be checked already when reading in the possible
  //decay channels.
  if (! type_a.is_hadron() || ! type_b.is_hadron() || ! type_c.is_hadron()) {
    printf("Warning: decay products A: %s B: %s C: %s\n",
           type_a.string().c_str(), type_b.string().c_str(),
           type_c.string().c_str());
  }
  printd("Note: Doing 1->3 decay!\n");

  /* Add three new particles */
  ParticleData new_particle_a, new_particle_b, new_particle_c;
  new_particle_a.set_pdgcode(type_a);
  new_particle_b.set_pdgcode(type_b);
  new_particle_c.set_pdgcode(type_c);

  FourVector momentum_resonance = particles->data(resonance_id).momentum();
  const double mass_a = particles->particle_type(type_a).mass(),
    mass_b = particles->particle_type(type_b).mass(),
    mass_c = particles->particle_type(type_c).mass(),
    mass_resonance = sqrt(momentum_resonance.Dot(momentum_resonance));

  /* mandelstam-s limits for pairs ab and bc */
  double s_ab_max = (mass_resonance - mass_c) * (mass_resonance - mass_c);
  double s_ab_min = (mass_a + mass_b) * (mass_a + mass_b);
  double s_bc_max = (mass_resonance - mass_a) * (mass_resonance - mass_a);
  double s_bc_min = (mass_b + mass_c) * (mass_b + mass_c);

  printd("s_ab limits: %g %g \n", s_ab_min, s_ab_max);
  printd("s_bc limits: %g %g \n", s_bc_min, s_bc_max);

  /* randomly pick values for s_ab and s_bc
   * until the pair is within the Dalitz plot */
  double dalitz_bc_max = 0.0, dalitz_bc_min = 1.0;
  double s_ab = 0.0, s_bc = 0.5;
  while (s_bc > dalitz_bc_max || s_bc < dalitz_bc_min) {
    s_ab = Random::uniform(s_ab_min, s_ab_max);
    s_bc = Random::uniform(s_bc_min, s_bc_max);
    double e_b_rest = (s_ab - mass_a * mass_a + mass_b * mass_b)
                           / (2 * sqrt(s_ab));
    double e_c_rest = (mass_resonance * mass_resonance - s_ab
                            - mass_c * mass_c) / (2 * sqrt(s_ab));
    dalitz_bc_max = (e_b_rest + e_c_rest) * (e_b_rest + e_c_rest)
      - (sqrt(e_b_rest * e_b_rest - mass_b * mass_b)
         - sqrt(e_c_rest * e_c_rest - mass_c * mass_c))
      * (sqrt(e_b_rest * e_b_rest - mass_b * mass_b)
         - sqrt(e_c_rest * e_c_rest - mass_c * mass_c));
    dalitz_bc_min = (e_b_rest + e_c_rest) * (e_b_rest + e_c_rest)
      - (sqrt(e_b_rest * e_b_rest - mass_b * mass_b)
         + sqrt(e_c_rest * e_c_rest - mass_c * mass_c))
      * (sqrt(e_b_rest * e_b_rest - mass_b * mass_b)
         + sqrt(e_c_rest * e_c_rest - mass_c * mass_c));
  }

  printd("s_ab: %g s_bc: %g min: %g max: %g\n",
         s_ab, s_bc, dalitz_bc_min, dalitz_bc_max);

  /* Compute energy and momentum magnitude */
  const double energy_a = (mass_resonance * mass_resonance + mass_a * mass_a
                           - s_bc) / (2 * mass_resonance);
  const double energy_c = (mass_resonance * mass_resonance + mass_c * mass_c
                           - s_ab) / (2 * mass_resonance);
  const double energy_b = (s_ab + s_bc - mass_a * mass_a - mass_c * mass_c)
                           / (2 * mass_resonance);
  const double momentum_a = sqrt(energy_a * energy_a - mass_a * mass_a);
  const double momentum_c = sqrt(energy_c * energy_c - mass_c * mass_c);
  const double momentum_b = sqrt(energy_b * energy_b - mass_b * mass_b);

  const double total_energy = momentum_resonance.x0();
  if (fabs(energy_a + energy_b + energy_c - total_energy) > really_small)
    printf("1->3 warning: Ea + Eb + Ec: %g Total E: %g",
           energy_a + energy_b + energy_c, total_energy);
  printd("Calculating the angles...\n");

  /* momentum_a direction is random */
  Angles phitheta;
  phitheta.distribute_isotropically();
  /* This is the angle of the plane of the three decay particles */
  new_particle_a.set_momentum(mass_a,
                              momentum_a * phitheta.x(),
                              momentum_a * phitheta.y(),
                              momentum_a * phitheta.z());

  /* Angle between a and b */
  double theta_ab = acos((energy_a * energy_b - 0.5 * (s_ab - mass_a * mass_a
                          - mass_b * mass_b)) / (momentum_a * momentum_b));
  printd("theta_ab: %g Ea: %g Eb: %g sab: %g pa: %g pb: %g\n",
         theta_ab, energy_a, energy_b, s_ab, momentum_a, momentum_b);
  bool phi_has_changed = phitheta.add_to_theta(theta_ab);
  new_particle_b.set_momentum(mass_b,
                              momentum_b * phitheta.x(),
                              momentum_b * phitheta.y(),
                              momentum_b * phitheta.z());

  /* Angle between b and c */
  double theta_bc = acos((energy_b * energy_c - 0.5 *(s_bc - mass_b * mass_b
                         - mass_c * mass_c)) / (momentum_b * momentum_c));
  printd("theta_bc: %g Eb: %g Ec: %g sbc: %g pb: %g pc: %g\n",
         theta_bc, energy_b, energy_c, s_bc, momentum_b, momentum_c);
  // pass information on whether phi has changed during the last adding
  // on to add_to_theta:
  phitheta.add_to_theta(theta_bc, phi_has_changed);
  new_particle_c.set_momentum(mass_c,
                              momentum_c * phitheta.x(),
                              momentum_c * phitheta.y(),
                              momentum_c * phitheta.z());

  /* Momentum check */
  double energy = new_particle_a.momentum().x0()
    + new_particle_b.momentum().x0() + new_particle_c.momentum().x0();
  double px = new_particle_a.momentum().x1() + new_particle_b.momentum().x1()
    + new_particle_c.momentum().x1();
  double py = new_particle_a.momentum().x2() + new_particle_b.momentum().x2()
    + new_particle_c.momentum().x2();
  double pz = new_particle_a.momentum().x3() + new_particle_b.momentum().x3()
    + new_particle_c.momentum().x3();

  if (fabs(energy - total_energy) > really_small)
    printf("1->3 energy not conserved! Before: %g After: %g\n",
           total_energy, energy);

  if (fabs(px) > really_small || fabs(py) > really_small
      || fabs(pz) > really_small)
    printf("1->3 momentum check failed. Total momentum: %g %g %g\n",
           px, py, pz);

  /* All decay products begin from the same point */
  FourVector decay_point = particles->data(resonance_id).position();
  new_particle_a.set_position(decay_point);
  new_particle_b.set_position(decay_point);
  new_particle_c.set_position(decay_point);

  int id_first_new = particles->add_data(new_particle_a);
  particles->add_data(new_particle_b);
  particles->add_data(new_particle_c);

  printd("p0: %g %g %g \n", new_particle_a.momentum().x0(),
         new_particle_b.momentum().x0(), new_particle_c.momentum().x0());
  printd("p1: %g %g %g \n", new_particle_a.momentum().x1(),
         new_particle_b.momentum().x1(), new_particle_c.momentum().x1());
  printd("p2: %g %g %g \n", new_particle_a.momentum().x2(),
         new_particle_b.momentum().x2(), new_particle_c.momentum().x2());
  printd("p3: %g %g %g \n", new_particle_a.momentum().x3(),
         new_particle_b.momentum().x3(), new_particle_c.momentum().x3());

  /* 3 new particles created; return the id of the first one */
  return id_first_new;
}


void DecayAction::choose_channel (Particles *particles) {
  const PdgCode pdgcode = particles->type(incoming_particles_[0]).pdgcode();

  /* Get the decay modes of this resonance */
  const std::vector<ProcessBranch> decaymodes
    = particles->decay_modes(pdgcode).decay_mode_list();
  /* Get the first decay mode and its branching ratio */
  std::vector<ProcessBranch>::const_iterator mode = decaymodes.begin();
  float cumulated_probability = mode->weight();
  /* Ratios of decay channels should add to 1; pick a random number
   * between 0 and 1 to select the decay mode to be used
   */
  double random_mode = Random::canonical();
  /* Keep adding to the probability until it exceeds the random value */
  while (random_mode > cumulated_probability &&  mode != decaymodes.end()) {
    cumulated_probability += mode->weight();
    ++mode;
  }

  outgoing_particles_ = mode->particle_list();
}


/**
 * Execute a decay process for the selected particle.
 *
 * Randomly select one of the decay modes of the particle
 * according to their relative weights. Then decay the particle
 * by calling function one_to_two or one_to_three.
 *
 * \param[in,out] particles Particles in the simulation.
 *
 * \return The ID of the first decay product.
 */
int DecayAction::resonance_decay (Particles *particles) {

  /* Decide for a particular decay channel. */
  choose_channel (particles);

  /* We found our decay branch, get the decay product pdgs and do the decay */
  size_t decay_particles = outgoing_particles_.size();
  int new_id_a = -1;
  switch (decay_particles) {
  case 2:
    new_id_a = one_to_two (particles);
    break;
  case 3:
    new_id_a = one_to_three (particles);
    break;
  default:
    printf("Warning: Not a 1->2 or 1->3 process!\n");
    printf("Number of decay particles: %zu \n", decay_particles);
    printf("Decay particles: ");
    for (size_t i = 0; i < decay_particles; i++) {
      printf("%s ", outgoing_particles_[i].pdgcode().string().c_str());
    }
    printf("\n");
  }
  return new_id_a;
}

void DecayAction::perform (Particles *particles, size_t &id_process) {
  FourVector velocity_CM;
  int id_a = incoming_particles_[0];

  /* Check if particle still exists. */
  if (is_valid(*particles)) {
    printf("ScatterAction::perform: ID %i not found!\n", id_a);
    return;
  }

  if (interaction_type_ != 2)
    printf("Decays warning: ID %i (%s) has process type %i.\n",
           id_a, particles->type(id_a).name().c_str(), interaction_type_);

  /* Save a copy of the initial state */
  ParticleData initial_data = particles->data(id_a);
  ParticleType initial_type = particles->type(id_a);

  printd("Process: Resonance decay. ");
  printd_momenta("Resonance momenta before decay", particles->data(id_a));

  /* boost to rest frame */
  velocity_CM.set_x0(1.0);
  velocity_CM.set_x1(particles->data(id_a).momentum().x1()
                      / particles->data(id_a).momentum().x0());
  velocity_CM.set_x2(particles->data(id_a).momentum().x2()
                      / particles->data(id_a).momentum().x0());
  velocity_CM.set_x3(particles->data(id_a).momentum().x3()
                      / particles->data(id_a).momentum().x0());
  particles->data_pointer(id_a)->set_momentum(
      particles->data(id_a).momentum().LorentzBoost(velocity_CM));
  particles->data_pointer(id_a)->set_position(
      particles->data(id_a).position().LorentzBoost(velocity_CM));

  printd_momenta("Boosted resonance momenta before decay",
                 particles->data(id_a));

  /* Save the highest id before decay */
  size_t old_max_id = particles->id_max();
  /* Do the decay; this returns the smallest new id */
  size_t id_new_a = resonance_decay (particles);
  /* There's going to be at least 2 new particles */
  size_t id_new_b = id_new_a + 1;

  printd_momenta("particle 1 momenta in lrf", particles->data(id_new_a));
  printd_momenta("particle 2 momenta in lrf", particles->data(id_new_b));

  boost_back_CM(particles->data_pointer(id_new_a),
                particles->data_pointer(id_new_b), &velocity_CM);

  /* How many new particles we have exactly */
  size_t new_particles = particles->id_max() - old_max_id;

  /* Check for the possible third final state particle */
  int id_new_c = -1;
  if (new_particles == 3) {
    id_new_c = id_new_b + 1;
    FourVector velocity = velocity_CM;
    velocity *= -1;
    velocity.set_x0(1.0);
    FourVector momentum_c(particles->data(id_new_c).momentum());
    FourVector position_c(particles->data(id_new_c).position());
    /* Boost the momenta back to lab frame */
    momentum_c = momentum_c.LorentzBoost(velocity);
    /* Boost the position back to lab frame */
    position_c = position_c.LorentzBoost(velocity);
    /* Write the oscar output for this particle */
    particles->data_pointer(id_new_c)->set_momentum(momentum_c);
    particles->data_pointer(id_new_c)->set_position(position_c);
  }

  printd_momenta("particle 1 momenta in comp", particles->data(id_new_a));
  printd_momenta("particle 2 momenta in comp", particles->data(id_new_b));

  FourVector final_momentum(particles->data(id_new_a).momentum()
    + particles->data(id_new_b).momentum());

  /* unset collision time for both particles + keep id + unset partner */
  particles->data_pointer(id_new_a)->set_collision_past(id_process);
  particles->data_pointer(id_new_b)->set_collision_past(id_process);
  if (new_particles == 3) {
    final_momentum += particles->data(id_new_c).momentum();
    particles->data_pointer(id_new_c)->set_collision_past(id_process);
  }
  printd("Particle map has now %zu elements. \n", particles->size());

  id_process++;

  /* Check momentum conservation */
  FourVector momentum_difference;
  momentum_difference += initial_data.momentum();
  momentum_difference -= final_momentum;
  if (fabs(momentum_difference.x0()) > really_small) {
    printf("Process %zu type %i particle %s decay to %s and %s ",
           id_process, interaction_type_, particles->type(id_a).name().c_str(),
           particles->type(id_new_a).name().c_str(),
           particles->type(id_new_b).name().c_str());
    if (new_particles == 3) {
      printf("and %s ", particles->type(id_new_c).name().c_str());
    }
    printf("time %g\n", initial_data.position().x0());
    printf("Warning: Interaction type %i E conservation violation %g\n",
           interaction_type_, momentum_difference.x0());
  }
  if (fabs(momentum_difference.x1()) > really_small)
    printf("Warning: Interaction type %i px conservation violation %g\n",
           interaction_type_, momentum_difference.x1());
  if (fabs(momentum_difference.x2()) > really_small)
    printf("Warning: Interaction type %i py conservation violation %g\n",
           interaction_type_, momentum_difference.x2());
  if (fabs(momentum_difference.x3()) > really_small)
    printf("Warning: Interaction type %i pz conservation violation %g\n",
           interaction_type_, momentum_difference.x3());

  /* Remove decayed particle */
  particles->remove(id_a);
  printd("ID %i has decayed and removed from the list.\n", id_a);
}



}
