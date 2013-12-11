/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_PROPAGATION_H_
#define SRC_INCLUDE_PROPAGATION_H_


class BoxBoundaryConditions;
class Particles;
class FourVector;

void propagate_particles(Particles *particles);

/* enforce periodic boundary conditions */
FourVector boundary_condition(FourVector position, bool *boundary_hit);

#endif  // SRC_INCLUDE_PROPAGATION_H_
