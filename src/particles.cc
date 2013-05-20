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

#include "include/FourVector.h"
#include "include/ParticleData.h"
#include "include/ParticleType.h"
#include "include/constants.h"
#include "include/box.h"
#include "include/outputroutines.h"

/* boost_COM - boost to center of momentum */
void boost_COM(ParticleData *particle1, ParticleData *particle2,
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

/* boost_back_COM - boost back from center of momentum */
void boost_back_COM(ParticleData *particle1, ParticleData *particle2,
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
  FourVector velocity_com;

  /* boost particles in center of momenta frame */
  boost_COM(&particle1, &particle2, &velocity_com);
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
  FourVector velocity_difference = particle1->momentum()
    / particle1->momentum().x0()
    - particle2->momentum() / particle2->momentum().x0();
  return - position_difference.DotThree(velocity_difference)
           / velocity_difference.DotThree(velocity_difference);
}

/* momenta_exchange - soft scattering */
void momenta_exchange(ParticleData *particle1, ParticleData *particle2,
  const float &particle1_mass, const float &particle2_mass) {
  double phi, theta;
  /* center of momentum hence equal for both particles */
  double momentum_radial = sqrt(particle1->momentum().x0()
    * particle1->momentum().x0() - particle1_mass * particle1_mass);
  printd("center of momenta 1: %g %g %g %g \n", particle1->momentum().x0(),
    particle1->momentum().x1(), particle1->momentum().x2(),
    particle1->momentum().x3());
  printd("center of momenta 2: %g %g %g %g \n", particle2->momentum().x0(),
    particle2->momentum().x1(), particle2->momentum().x2(),
    particle2->momentum().x3());

  /* particle exchange momenta and scatter to random direction */
  phi =  2 * M_PI * drand48();
  theta = M_PI * drand48();
  particle1->set_momentum(particle1_mass,
     momentum_radial * cos(phi) * sin(theta),
     momentum_radial * sin(phi) * sin(theta), momentum_radial * cos(theta));
  particle2->set_momentum(particle2_mass,
     -momentum_radial * cos(phi) * sin(theta),
     -momentum_radial * sin(phi) * sin(theta), -momentum_radial * cos(theta));

  printd("exchanged momenta 1: %g %g %g %g \n", particle1->momentum().x0(),
    particle1->momentum().x1(), particle1->momentum().x2(),
    particle1->momentum().x3());
  printd("exchanged momenta 2: %g %g %g %g \n", particle2->momentum().x0(),
    particle2->momentum().x1(), particle2->momentum().x2(),
    particle2->momentum().x3());
}
