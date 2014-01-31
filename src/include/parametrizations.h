/*
 *
 *    Copyright (c) 2013
 *      SMASH Team
 * 
 *    GNU General Public License (GPLv3 or later)
 *
 */
#ifndef SRC_INCLUDE_PARAMETRIZATIONS_H_
#define SRC_INCLUDE_PARAMETRIZATIONS_H_

/* pp elastic cross section parametrization */
float pp_elastic(double p_lab, double mandelstam_s, float nucleon_mass);

/* pp total cross section parametrization */
float pp_total(double p_lab);

/* np elastic cross section parametrization */
float np_elastic(double p_lab, double mandelstam_s, float nucleon_mass);

/* np total cross section parametrization */
float np_total(double p_lab);

/* ppbar elastic cross section parametrization */
float ppbar_elastic(double p_lab);

/* ppbar total cross section parametrization */
float ppbar_total(double p_lab);

#endif  // SRC_INCLUDE_PARAMETRIZATIONS_H_
