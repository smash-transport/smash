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

/* calculate_minimum_mass
 * - calculate the minimum rest energy the resonance must have
 * to be able to decay through any of its decay channels
 * NB: This function assumes stable decay products!
 */
float calculate_minimum_mass(Particles *particles, int pdgcode);

/* resonance_cross_section - energy-dependent cross section
 * for producing a resonance
 */
std::map<int, double> resonance_cross_section(
  const ParticleData &particle1, const ParticleData &particle2,
  const ParticleType &type_particle1, const ParticleType &type_particle2,
  Particles *particles);

/* 2->1 resonance formation process */
int resonance_formation(Particles *particles, int particle_id, int other_id,
  int resonance_type);

#endif  // SRC_INCLUDE_RESONANCES_H_
