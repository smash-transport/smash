/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_RESONANCES_H_
#define SRC_INCLUDE_RESONANCES_H_

#include <cstdio>

#include <map>
#include <vector>

/* necessary forward declarations */
class ParticleData;
class ParticleType;

/* resonance_cross_section - energy-dependent cross section
 * for producing a resonance
 */
std::map<int, double> resonance_cross_section(
  const ParticleData &particle1, const ParticleData &particle2,
  const ParticleType &type_particle1, const ParticleType &type_particle2,
  std::vector<ParticleType> *type_list);

/* 1->2 resonance decay process */
size_t resonance_decay(std::map<int, ParticleData> *particles,
  std::vector<ParticleType> *types, std::map<int, int> *map_type,
  int *particle_id, int *id_max);

/* 2->1 resonance formation process */
size_t resonance_formation(std::map<int, ParticleData> *particles,
  std::vector<ParticleType> *types, std::map<int, int> *map_type,
  int *particle_id, int *other_id, int resonance_type, int *id_max);

#endif  // SRC_INCLUDE_RESONANCES_H_
