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

  /* We have no resonances with charge > 1 */
  if (abs(charge1 + charge2) > 0)
    return 0.0;

  /* Otherwise, total charge defines the type of resonance */
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

  double cross_section = 1.0;

  return cross_section;
}
