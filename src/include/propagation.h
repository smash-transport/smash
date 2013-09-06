/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_PROPAGATION_H_
#define SRC_INCLUDE_PROPAGATION_H_


class Box;
class Parameters;
class Particles;

void propagate_particles(Particles *particles, Parameters const &parameters,
                         Box const &box);

#endif  // SRC_INCLUDE_PROPAGATION_H_
