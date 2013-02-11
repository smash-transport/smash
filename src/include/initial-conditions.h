/*
 *    Copyright (c) 2012 maximilian attems <attems@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_INITIAL_CONDITIONS_H_
#define SRC_INCLUDE_INITIAL_CONDITIONS_H_

/* forward declaration */
class ParticleData;

ParticleData* initial_conditions(ParticleData *particles, int &number);

#endif  // SRC_INCLUDE_INITIAL_CONDITIONS_H_
