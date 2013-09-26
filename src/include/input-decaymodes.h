/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_INPUT_DECAYMODES_H_
#define SRC_INCLUDE_INPUT_DECAYMODES_H_

#include <vector>

/* forward declarations */
class Particles;

extern const char *separator;

/* read input file particle types */
void input_decaymodes(Particles *particles, char *path);

#endif  // SRC_INCLUDE_INPUT_DECAYMODES_H_
