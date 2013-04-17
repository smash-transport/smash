/*
 *    Copyright (c) 2012 maximilian attems <attems@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_INITIAL_CONDITIONS_H_
#define SRC_INCLUDE_INITIAL_CONDITIONS_H_

/* forward declarations */
class ParticleData;
class ParticleType;
class box;

ParticleData* initial_conditions(ParticleData *particles,
  ParticleType *particle_type, int &number, box *box);

#endif  // SRC_INCLUDE_INITIAL_CONDITIONS_H_
