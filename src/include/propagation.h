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
class Laboratory;
class Particles;
class FourVector;

void propagate_particles(Particles *particles, Laboratory const &parameters,
                         Box const &box);

/* enforce periodic boundary conditions */
FourVector boundary_condition(FourVector position, const Box &box,
                              bool *boundary_hit);

#endif  // SRC_INCLUDE_PROPAGATION_H_
