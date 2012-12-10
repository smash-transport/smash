/*
 *    Copyright (c) 2012 maximilian attems <attems@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_PARTICLEDATA_H_
#define SRC_INCLUDE_PARTICLEDATA_H_

class ParticleData {
  private:
    /* Each particle has a unique identifier */
    int id;
    /* longitudinal momenta */
    double momenta_l;
    /* tranverse momenta */
    double momenta_t;
};

#endif  // SRC_INCLUDE_PARTICLEDATA_H_
