/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_PARTICLES_H_
#define SRC_INCLUDE_PARTICLES_H_

/* necessary forward declarations */
class ParticleData;
class FourVector;

/* boost_COM - boost to center of momentum */
void boost_COM(ParticleData *particle1, ParticleData *particle2,
  FourVector *velocity);

/* boost_from_COM - boost back from center of momentum */
void boost_from_COM(ParticleData *particle1, ParticleData *particle2,
  FourVector *velocity_orig);

/* particle_distance - measure distance between two particles */
double particle_distance(ParticleData *particle_orig1,
  ParticleData *particle_orig2);

/* time_collision - measure collision time of two particles */
double collision_time(ParticleData *particle1, ParticleData *particle2);

/* momenta_exchange - soft scattering */
void momenta_exchange(ParticleData *particle1, ParticleData *particle2,
  const float &particle1_mass, const float &particle2_mass);
#endif  // SRC_INCLUDE_PARTICLES_H_
