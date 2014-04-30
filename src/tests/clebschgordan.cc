/*
 *    Copyright (c) 2013-2014
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */

#include <cmath>

#include "include/clebschgordan.h"

using namespace Smash;

int main() {
  /* spins are two times the actual values,
   * so j = 1 for spin-1/2 particle, 2 for spin-1 particle, etc.
   * Ordering of spins in array:
   *  0: j1, 1: j2, 2: j3, 3: m1, 4: m2, 5: m3
  */
  int spins[7][6];
  double correct_coefficient[7];
  const double tolerance = 1e-6;
  ClebschGordan cgfactor;

  spins[0][0] = 1;
  spins[0][1] = 1;
  spins[0][2] = 2;
  spins[0][3] = 1;
  spins[0][4] = 1;
  spins[0][5] = 2;
  correct_coefficient[0] = 1.0;

  spins[1][0] = 1;
  spins[1][1] = 1;
  spins[1][2] = 2;
  spins[1][3] = 1;
  spins[1][4] = -1;
  spins[1][5] = 0;
  correct_coefficient[1] = 1 / sqrt(2.0);

  spins[2][0] = 2;
  spins[2][1] = 1;
  spins[2][2] = 1;
  spins[2][3] = 2;
  spins[2][4] = -1;
  spins[2][5] = 1;
  correct_coefficient[2] = sqrt(2.0 / 3.0);

  spins[3][0] = 2;
  spins[3][1] = 1;
  spins[3][2] = 3;
  spins[3][3] = -2;
  spins[3][4] = 1;
  spins[3][5] = -1;
  correct_coefficient[3] = sqrt(1.0 / 3.0);

  spins[4][0] = 2;
  spins[4][1] = 2;
  spins[4][2] = 2;
  spins[4][3] = 0;
  spins[4][4] = 2;
  spins[4][5] = 2;
  correct_coefficient[4] = -1 / sqrt(2.0);

  spins[5][0] = 2;
  spins[5][1] = 2;
  spins[5][2] = 2;
  spins[5][3] = 0;
  spins[5][4] = 0;
  spins[5][5] = 0;
  correct_coefficient[5] = 0.0;

  spins[6][0] = 2;
  spins[6][1] = 2;
  spins[6][2] = 4;
  spins[6][3] = 2;
  spins[6][4] = -2;
  spins[6][5] = 0;
  correct_coefficient[6] = 1 / sqrt(6.0);

  for (int i = 0; i < 7; i++) {
    double cgresult = cgfactor(spins[i][0], spins[i][1], spins[i][2],
                        spins[i][3], spins[i][4], spins[i][5]);
    if (fabs(cgresult - correct_coefficient[i]) > tolerance)
      return -i;
  }
  return 0;
}
