/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_DISTRIBUTIONS_H_
#define SRC_INCLUDE_DISTRIBUTIONS_H_

class Box;
class ParticleType;

/* Breit-Wigner distribution for calculating resonance
 * production probability
 */
double breit_wigner(double mandelstam_s, float resonance_mass,
                          float resonance_width);

/* density_integrand - Maxwell-Boltzmann distribution */
double inline density_integrand(double momentum, double temp, double mass);

/* sample_momenta - return thermal momenta */
double sample_momenta(Box *box, ParticleType type);

#endif  // SRC_INCLUDE_DISTRIBUTIONS_H_
