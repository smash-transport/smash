/*
 *
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */
#include "include/particles.h"

#include <cstdio>
#include <cstring>

#include "include/constants.h"
#include "include/distributions.h"
#include "include/FourVector.h"
#include "include/outputroutines.h"
#include "include/ParticleData.h"
#include "include/ParticleType.h"

/* boost_CM - boost to center of momentum */
void boost_CM(ParticleData *particle1, ParticleData *particle2,
  FourVector *velocity) {
  FourVector momentum1(particle1->momentum()), momentum2(particle2->momentum());
  FourVector position1(particle1->position()), position2(particle2->position());
  double cms_energy = momentum1.x0() + momentum2.x0();

  // CMS 4-velocity
  velocity->set_x0(1.0);
  velocity->set_x1((momentum1.x1() + momentum2.x1()) / cms_energy);
  velocity->set_x2((momentum1.x2() + momentum2.x2()) / cms_energy);
  velocity->set_x3((momentum1.x3() + momentum2.x3()) / cms_energy);

  // Boost the momenta into CMS frame
  momentum1 = momentum1.LorentzBoost(*velocity);
  momentum2 = momentum2.LorentzBoost(*velocity);

  // Boost the positions into CMS frame
  position1 = position1.LorentzBoost(*velocity);
  position2 = position2.LorentzBoost(*velocity);

  particle1->set_momentum(momentum1);
  particle1->set_position(position1);
  particle2->set_momentum(momentum2);
  particle2->set_position(position2);
}

/* boost_back_CM - boost back from center of momentum */
void boost_back_CM(ParticleData *particle1, ParticleData *particle2,
  FourVector *velocity_orig) {
  FourVector momentum1(particle1->momentum()), momentum2(particle2->momentum());
  FourVector position1(particle1->position()), position2(particle2->position());
  FourVector velocity = *velocity_orig;

  /* To boost back set 1 + velocity */
  velocity *= -1;
  velocity.set_x0(1.0);

  /* Boost the momenta back to lab frame */
  momentum1 = momentum1.LorentzBoost(velocity);
  momentum2 = momentum2.LorentzBoost(velocity);

  /* Boost the positions back to lab frame */
  position1 = position1.LorentzBoost(velocity);
  position2 = position2.LorentzBoost(velocity);

  particle1->set_momentum(momentum1);
  particle1->set_position(position1);
  particle2->set_momentum(momentum2);
  particle2->set_position(position2);
}

/* particle_distance - measure distance between two particles
 *                     in center of momentum
 */
double particle_distance(ParticleData *particle_orig1,
  ParticleData *particle_orig2) {
  /* Copy the particles in order to boost them and to forget the copy */
  ParticleData particle1 = *particle_orig1, particle2 = *particle_orig2;
  FourVector velocity_CM;

  /* boost particles in center of momenta frame */
  boost_CM(&particle1, &particle2, &velocity_CM);
  FourVector position_difference = particle1.position() - particle2.position();
  printd("Particle %d<->%d position difference: %g %g %g %g [fm]\n",
    particle1.id(), particle2.id(), position_difference.x0(),
    position_difference.x1(), position_difference.x2(),
    position_difference.x3());

  FourVector momentum_difference = particle1.momentum() - particle2.momentum();
  printd("Particle %d<->%d momentum difference: %g %g %g %g [fm]\n",
    particle1.id(), particle2.id(), momentum_difference.x0(),
    momentum_difference.x1(), momentum_difference.x2(),
    momentum_difference.x3());
  /* zero momentum leads to infite distance */
  if (momentum_difference.x1() == 0 || momentum_difference.x2() == 0
      || momentum_difference.x3() == 0)
    return  - position_difference.DotThree(position_difference);

  /* UrQMD squared distance criteria:
   * arXiv:nucl-th/9803035 (3.27): in center of momemtum frame
   * position of particle a: x_a
   * position of particle b: x_b
   * velocity of particle a: v_a
   * velocity of particle b: v_b
   * d^2_{coll} = (x_a - x_b)^2 - ((x_a - x_a) . (v_a - v_b))^2 / (v_a - v_b)^2
   */
  return - position_difference.DotThree(position_difference)
    + position_difference.DotThree(momentum_difference)
      * position_difference.DotThree(momentum_difference)
      / momentum_difference.DotThree(momentum_difference);
}

/* time_collision - measure collision time of two particles */
double collision_time(ParticleData *particle1, ParticleData *particle2) {
  /* UrQMD collision time
   * arXiv:1203.4418 (5.15): in computational frame
   * position of particle a: x_a
   * position of particle b: x_b
   * momentum of particle a: p_a
   * momentum of particle b: p_b
   * t_{coll} = - (x_a - x_b) . (p_a - p_b) / (p_a - p_b)^2
   */
  FourVector position_difference = particle1->position()
    - particle2->position();
  printd("Particle %d<->%d position difference: %g %g %g %g [fm]\n",
    particle1->id(), particle2->id(), position_difference.x0(),
    position_difference.x1(), position_difference.x2(),
    position_difference.x3());
  FourVector velocity_difference = particle1->momentum()
    / particle1->momentum().x0()
    - particle2->momentum() / particle2->momentum().x0();
  printd("Particle %d<->%d velocity difference: %g %g %g %g [fm]\n",
    particle1->id(), particle2->id(), velocity_difference.x0(),
    velocity_difference.x1(), velocity_difference.x2(),
    velocity_difference.x3());
  /* zero momentum leads to infite distance, particles are not approaching */
  if (velocity_difference.x1() == 0.0 || velocity_difference.x2() == 0.0
      || velocity_difference.x3() == 0.0)
    return -1.0;
  return - position_difference.DotThree(velocity_difference)
           / velocity_difference.DotThree(velocity_difference);
}

/* momenta_exchange - soft scattering */
void momenta_exchange(ParticleData *particle1, ParticleData *particle2,
  const float &particle1_mass, const float &particle2_mass) {
  /* debug output */
  printd("center of momenta 1: %g %g %g %g \n", particle1->momentum().x0(),
    particle1->momentum().x1(), particle1->momentum().x2(),
    particle1->momentum().x3());
  printd("center of momenta 2: %g %g %g %g \n", particle2->momentum().x0(),
    particle2->momentum().x1(), particle2->momentum().x2(),
    particle2->momentum().x3());

  /* center of momentum hence this is equal for both particles */
  const double momentum_radial = sqrt(particle1->momentum().x0()
    * particle1->momentum().x0() - particle1_mass * particle1_mass);
  printd("Particle 1: momentum %g mass %g \n", particle1->momentum().x0(),
    particle1_mass);
  /* particle exchange momenta and scatter to random direction */
  const double phi =  2.0 * M_PI * drand48();
  const double cos_theta = -1.0 + 2.0 * drand48();
  const double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
  printd("Random momentum: %g %g %g %g \n", momentum_radial, phi, cos_theta,
    sin_theta);
  const FourVector momentum1(sqrt(particle1_mass * particle1_mass
    + momentum_radial * momentum_radial),
     momentum_radial * cos(phi) * sin_theta,
     momentum_radial * sin(phi) * sin_theta, momentum_radial * cos_theta);
  particle1->set_momentum(momentum1);
  const FourVector momentum2(sqrt(particle2_mass * particle2_mass
    + momentum_radial * momentum_radial),
    - momentum_radial * cos(phi) * sin_theta,
    - momentum_radial * sin(phi) * sin_theta, -momentum_radial * cos_theta);
  particle2->set_momentum(momentum2);

  /* debug output */
  printd("exchanged momenta 1: %g %g %g %g \n", particle1->momentum().x0(),
    particle1->momentum().x1(), particle1->momentum().x2(),
    particle1->momentum().x3());
  printd("exchanged momenta 2: %g %g %g %g \n", particle2->momentum().x0(),
    particle2->momentum().x1(), particle2->momentum().x2(),
    particle2->momentum().x3());
}

/* resonance_cross_section - energy-dependent cross section
 * for producing a resonance
 */

double resonance_cross_section(ParticleData *particle1, ParticleData *particle2,
  ParticleType *type_particle1, ParticleType *type_particle2,
  std::vector<ParticleType> *type_list) {
  const int charge1 = (*type_particle1).charge(),
    charge2 = (*type_particle2).charge();

  /* Total charge defines the type of resonance */
  /* We have no resonances with charge > 1 */
  if (abs(charge1 + charge2) > 1)
    return 0.0;

  std::string resonance_name;
  if (charge1 + charge2 == 1)
    resonance_name = "rho+";
  else if (charge1 + charge2 == -1)
    resonance_name = "rho-";
  else
    resonance_name = "rho0";

  /* Find the width and mass of the desired resonance */
  float resonance_width = -1.0, resonance_mass = 0.0;
  size_t type_index = 0;
  while (resonance_width < 0 && type_index < (*type_list).size()) {
    if (strcmp((*type_list)[type_index].name().c_str(),
               resonance_name.c_str()) == 0) {
      resonance_width = (*type_list)[type_index].width();
      resonance_mass = (*type_list)[type_index].mass();
      printd("Found resonance %s with mass %f and width %f.\n",
             resonance_name.c_str(), resonance_mass, resonance_width);
      printd("Original pions: %s %s Charges: %i %i \n",
             (*type_particle1).name().c_str(), (*type_particle2).name().c_str(),
             (*type_particle1).charge(), (*type_particle2).charge());
    }
    type_index++;
  }

  /* If there was no such resonance in the list, return 0 */
  if (resonance_width < 0.0)
    return 0.0;

  /* Mandelstam s = (p_a + p_b)^2 = square of CMS energy */
  const double mandelstam_s =
       ( (*particle1).momentum() + (*particle2).momentum()).Dot(
         (*particle1).momentum() + (*particle2).momentum() );

  /* XXX: Arbitrary proportionality constant */
  const double proportionality = 10.0;

  /* Calculate resonance production cross section
   * using the Breit-Wigner distribution as probability amplitude
   */
  return proportionality * breit_wigner(mandelstam_s,
          resonance_mass, resonance_width) / mandelstam_s;
}

/* 1->2 resonance decay process */
size_t resonance_decay(std::vector<ParticleData> *particles,
  std::vector<ParticleType> *types, std::map<int, int> *map_type,
  int *particle_id) {
  /* Add two new particles */
  size_t old_size = (*particles).size();
  size_t new_id_a = old_size, new_id_b = old_size + 1;
  (*particles).resize(old_size + 2);
  (*particles)[new_id_a].set_id(new_id_a);
  (*particles)[new_id_b].set_id(new_id_b);

  const int charge = (*types)[(*map_type)[*particle_id]].charge();

  std::string name_a, name_b;
  if (charge == 0) {
    name_a = "pi+";
    name_b = "pi-";
  } else if ( charge == 1 ) {
    name_a = "pi+";
    name_b = "pi0";
  } else if ( charge == -1 ) {
    name_a = "pi-";
    name_b = "pi0";
  }

  /* Find the desired particle types */
  bool not_found_a = true, not_found_b = true;
  size_t type_index = 0;
  while ( (not_found_a || not_found_b) && type_index < (*types).size() ) {
    if (strcmp((*types)[type_index].name().c_str(),
               name_a.c_str()) == 0) {
      printd("Found particle %s.\n", name_a.c_str());
      (*map_type)[new_id_a] = type_index;
      not_found_a = false;
    }
    if (strcmp((*types)[type_index].name().c_str(),
               name_b.c_str()) == 0) {
      printd("Found particle %s.\n", name_b.c_str());
      (*map_type)[new_id_b] = type_index;
      not_found_b = false;
    }
    type_index++;
  }

  const double total_energy = ((*particles)[*particle_id]).momentum().x0();
  double mass_a = (*types)[(*map_type)[new_id_a]].mass(),
    mass_b = (*types)[(*map_type)[new_id_b]].mass();
  double energy_a = (total_energy * total_energy
                     + mass_a * mass_a - mass_b * mass_b)
                    / (2.0 * total_energy);

  double momentum_radial = sqrt(energy_a * energy_a - mass_a * mass_a);
  /* phi in the range from [0, 2 * pi) */
  double phi = 2.0 * M_PI * drand48();
  /* cos(theta) in the range from [-1.0, 1.0) */
  double cos_theta = -1.0 + 2.0 * drand48();
  double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
  printd("Particle %lu radial momenta %g phi %g cos_theta %g\n", new_id_a,
         momentum_radial, phi, cos_theta);
  (*particles)[new_id_a].set_momentum(mass_a,
      momentum_radial * cos(phi) * sin_theta,
      momentum_radial * sin(phi) * sin_theta,
      momentum_radial * cos_theta);
  (*particles)[new_id_b].set_momentum(mass_b,
    - (*particles)[new_id_a].momentum().x1(),
    - (*particles)[new_id_a].momentum().x2(),
    - (*particles)[new_id_a].momentum().x3());

  /* Both decay products begin from the same point */
  FourVector decay_point = (*particles)[*particle_id].position();
  (*particles)[new_id_a].set_position(decay_point);
  (*particles)[new_id_b].set_position(decay_point);

  /* No collision yet */
  (*particles)[new_id_a].set_collision(-1, 0, -1);
  (*particles)[new_id_b].set_collision(-1, 0, -1);

  printd("Created %s and %s with IDs %lu and %lu \n",
  (*types)[(*map_type)[new_id_a]].name().c_str(),
  (*types)[(*map_type)[new_id_b]].name().c_str(), new_id_a, new_id_b);

  return new_id_a;
}

/* 2->1 resonance formation process */
size_t resonance_formation(std::vector<ParticleData> *particles,
  std::vector<ParticleType> *types, std::map<int, int> *map_type,
                         int *particle_id, int *other_id) {
  /* Add a new particle */
  size_t old_size = (*particles).size();
  size_t new_id = old_size;
  (*particles).resize(old_size + 1);
  (*particles)[new_id].set_id(new_id);

  /* Which resonance is formed */
  const int charge1 = (*types)[(*map_type)[*particle_id]].charge(),
    charge2 = (*types)[(*map_type)[*other_id]].charge();

  std::string resonance_name;
  if (charge1 + charge2 == 1)
    resonance_name = "rho+";
  else if (charge1 + charge2 == -1)
    resonance_name = "rho-";
  else
    resonance_name = "rho0";

  /* Find the desired resonance */
  bool not_found = true;
  size_t type_index = 0;
  while (not_found && type_index < (*types).size()) {
    if (strcmp((*types)[type_index].name().c_str(),
               resonance_name.c_str()) == 0) {
      printd("Found resonance %s.\n", resonance_name.c_str());
      printd("Parent particles: %s %s \n",
             (*types)[(*map_type)[*particle_id]].name().c_str(),
             (*types)[(*map_type)[*other_id]].name().c_str());
      (*map_type)[new_id] = type_index;
      not_found = false;
    }
    type_index++;
  }

  /* Initialize the timer for unstable particle */
  (*particles)[new_id].set_lifetime(0.0);

  /* Center-of-momentum frame of initial particles
   * is the rest frame of the resonance
   */
  const double energy = (*particles)[*particle_id].momentum().x0()
    + (*particles)[*other_id].momentum().x0();
  /* We use fourvector to set 4-momentum, as setting it
   * with doubles requires that particle is on
   * mass shell, which is not generally true for resonances
   */
  FourVector resonance_momentum(energy, 0.0, 0.0, 0.0);
  (*particles)[new_id].set_momentum(resonance_momentum);

  printd("Momentum of the new particle: %g %g %g %g \n",
    (*particles)[new_id].momentum().x0(), (*particles)[new_id].momentum().x1(),
    (*particles)[new_id].momentum().x2(), (*particles)[new_id].momentum().x3());

  /* The real position should be between parents in the computational frame! */
  (*particles)[new_id].set_position(1.0, 0.0, 0.0, 0.0);

  /* No collision yet */
  (*particles)[new_id].set_collision(-1, 0, -1);

  printd("Created %s with ID %i \n",
         (*types)[(*map_type)[new_id]].name().c_str(),
         (*particles)[new_id].id());

  return new_id;
}
