/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_RESONANCES_H_
#define SRC_INCLUDE_RESONANCES_H_

#include <map>

/* necessary forward declarations */
class Particles;
class ParticleData;
class ParticleType;

/* resonance_cross_section - energy-dependent cross section
 * for producing a resonance
 */
std::map<int, double> resonance_cross_section(
  const ParticleData &particle1, const ParticleData &particle2,
  const ParticleType &type_particle1, const ParticleType &type_particle2,
  const Particles &particles);

/* 1->2 process kinematics */
int one_to_two(Particles *particles, int resonance_id, int type_a, int type_b);

/* 1->3 process kinematics */
int one_to_three(Particles *particles, int resonance_id,
                 int type_a, int type_b, int type_c);

/* resonance decay process */
int resonance_decay(Particles *particles, int particle_id);

/* 2->1 resonance formation process */
int resonance_formation(Particles *particles, int particle_id, int other_id,
  int resonance_type);

#endif  // SRC_INCLUDE_RESONANCES_H_
