/*
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_EXPERIMENTPARAMETERS_H_
#define SRC_INCLUDE_EXPERIMENTPARAMETERS_H_

struct ExperimentParameters {
  /* number of test particle */
  int testparticles = 1;
  /* temporal time step */
  float eps = 0.001f;
  /* cross section of the elastic scattering */
  float cross_section = 10.0f;
};

#endif  // SRC_INCLUDE_EXPERIMENTPARAMETERS_H_
