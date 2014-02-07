/*
 *
 *    Copyright (c) 2013
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include <cmath>
#include <cstdio>

#include "include/angles.h"
#include "include/constants.h"
#include "include/distributions.h"
#include "include/fourvector.h"
#include "include/macros.h"
#include "include/outputroutines.h"
#include "include/particledata.h"
#include "include/particles.h"

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
  if (fabs(momentum_difference.x1()) < really_small
      && fabs(momentum_difference.x2()) < really_small
      && fabs(momentum_difference.x3()) < really_small)
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
double collision_time(const ParticleData &particle1,
  const ParticleData &particle2) {
  /* UrQMD collision time
   * arXiv:1203.4418 (5.15): in computational frame
   * position of particle a: x_a
   * position of particle b: x_b
   * momentum of particle a: p_a
   * momentum of particle b: p_b
   * t_{coll} = - (x_a - x_b) . (p_a - p_b) / (p_a - p_b)^2
   */
  FourVector position_difference = particle1.position()
    - particle2.position();
  printd("Particle %d<->%d position difference: %g %g %g %g [fm]\n",
    particle1.id(), particle2.id(), position_difference.x0(),
    position_difference.x1(), position_difference.x2(),
    position_difference.x3());
  FourVector velocity_difference = particle1.momentum()
    / particle1.momentum().x0()
    - particle2.momentum() / particle2.momentum().x0();
  printd("Particle %d<->%d velocity difference: %g %g %g %g [fm]\n",
    particle1.id(), particle2.id(), velocity_difference.x0(),
    velocity_difference.x1(), velocity_difference.x2(),
    velocity_difference.x3());
  /* zero momentum leads to infite distance, particles are not approaching */
  if (fabs(velocity_difference.x1()) < really_small
      && fabs(velocity_difference.x2()) < really_small
      && fabs(velocity_difference.x3()) < really_small)
    return -1.0;
  return - position_difference.DotThree(velocity_difference)
           / velocity_difference.DotThree(velocity_difference);
}

/* momenta_exchange - soft scattering */
void momenta_exchange(ParticleData *particle1, ParticleData *particle2) {
  /* debug output */
  printd_momenta("center of momenta 1", *particle1);
  printd_momenta("center of momenta 2", *particle2);

  /* center of momentum hence this is equal for both particles */
  const double momentum_radial = sqrt(particle1->momentum().x1()
    * particle1->momentum().x1() + particle1->momentum().x2() *
    particle1->momentum().x2() + particle1->momentum().x3() *
    particle1->momentum().x3());

  /* particle exchange momenta and scatter to random direction */
  /* XXX: Angles should be sampled from differential cross section
   * of this process
   */
  Angles phitheta;
  phitheta.distribute_isotropically();
  printd("Random momentum: %g %g %g %g \n", momentum_radial, phitheta.phi(),
        phitheta.costheta(), phitheta.sintheta());

  /* Only direction of 3-momentum, not magnitude, changes in CM frame,
   * thus particle energies remain the same (Lorentz boost will change them for
   * computational frame, however)
   */
  const FourVector momentum1(particle1->momentum().x0(),
     momentum_radial * phitheta.x(),
     momentum_radial * phitheta.y(),
     momentum_radial * phitheta.z());
  particle1->set_momentum(momentum1);
  const FourVector momentum2(particle2->momentum().x0(),
    - momentum_radial * phitheta.x(),
    - momentum_radial * phitheta.y(),
    - momentum_radial * phitheta.z());
  particle2->set_momentum(momentum2);

  /* debug output */
  printd_momenta("exchanged momenta 1", *particle1);
  printd_momenta("exchanged momenta 2", *particle2);
}

void sample_cms_momenta(ParticleData *particle1, ParticleData *particle2,
                        const double cms_energy, const double mass1,
                        const double mass2) {
  double energy1 = (cms_energy * cms_energy + mass1 * mass1 - mass2 * mass2) /
                   (2.0 * cms_energy);
  double momentum_radial = sqrt(energy1 * energy1 - mass1 * mass1);
  if (!(momentum_radial > 0.0))
    printf("Warning: radial momenta %g \n", momentum_radial);
  /* XXX: Angles should be sampled from differential cross section
   * of this process
   */
  Angles phitheta;
  phitheta.distribute_isotropically();
  if (!(energy1 > mass1)) {
    printf("Particle %d radial momenta %g phi %g cos_theta %g\n",
           particle1->pdgcode(), momentum_radial, phitheta.phi(),
           phitheta.costheta());
    printf("Etot: %g m_a: %g m_b %g E_a: %g\n", cms_energy, mass1, mass2,
           energy1);
  }
  /* We use fourvector to set 4-momentum, as setting it
   * with doubles requires that particle uses its
   * pole mass, which is not generally true for resonances
   */
  FourVector momentum1(energy1, momentum_radial * phitheta.x(),
                       momentum_radial * phitheta.y(),
                       momentum_radial * phitheta.z());
  particle1->set_momentum(momentum1);

  FourVector momentum2(cms_energy - energy1, -momentum_radial * phitheta.x(),
                       -momentum_radial * phitheta.y(),
                       -momentum_radial * phitheta.z());
  particle2->set_momentum(momentum2);

  printd("p0: %g %g \n", momentum1.x0(), momentum2.x0());
  printd("p1: %g %g \n", momentum1.x1(), momentum2.x1());
  printd("p2: %g %g \n", momentum1.x2(), momentum2.x2());
  printd("p3: %g %g \n", momentum1.x3(), momentum2.x3());
}
