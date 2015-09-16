/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "../include/clebschgordan.h"

/* spins are two times the actual values,
 * so j = 1 for spin-1/2 particle, 2 for spin-1 particle, etc.
 * Ordering of spins in array:
 *  0: j1, 1: j2, 2: j3, 3: m1, 4: m2, 5: m3
 */
int spin[7][3];
int spinz[7][3];
float correct_coefficient[7];

TEST(coefficient) {
  spin[0][0] = 1;
  spin[0][1] = 1;
  spin[0][2] = 2;
  spinz[0][0] = 1;
  spinz[0][1] = 1;
  spinz[0][2] = 2;
  correct_coefficient[0] = 1.0;

  spin[1][0] = 1;
  spin[1][1] = 1;
  spin[1][2] = 2;
  spinz[1][0] = 1;
  spinz[1][1] = -1;
  spinz[1][2] = 0;
  correct_coefficient[1] = 1 / sqrt(2.0);

  spin[2][0] = 2;
  spin[2][1] = 1;
  spin[2][2] = 1;
  spinz[2][0] = 2;
  spinz[2][1] = -1;
  spinz[2][2] = 1;
  correct_coefficient[2] = sqrt(2.0 / 3.0);

  spin[3][0] = 2;
  spin[3][1] = 1;
  spin[3][2] = 3;
  spinz[3][0] = -2;
  spinz[3][1] = 1;
  spinz[3][2] = -1;
  correct_coefficient[3] = sqrt(1.0 / 3.0);

  spin[4][0] = 2;
  spin[4][1] = 2;
  spin[4][2] = 2;
  spinz[4][0] = 0;
  spinz[4][1] = 2;
  spinz[4][2] = 2;
  correct_coefficient[4] = -1 / sqrt(2.0);

  spin[5][0] = 2;
  spin[5][1] = 2;
  spin[5][2] = 2;
  spinz[5][0] = 0;
  spinz[5][1] = 0;
  spinz[5][2] = 0;
  correct_coefficient[5] = 0.0;

  spin[6][0] = 2;
  spin[6][1] = 2;
  spin[6][2] = 4;
  spinz[6][0] = 2;
  spinz[6][1] = -2;
  spinz[6][2] = 0;
  correct_coefficient[6] = 1 / sqrt(6.0);
  for (int i = 0; i < 7; i++) {
    float cg = Smash::clebsch_gordan(spin[i][0], spin[i][1], spin[i][2],
                                     spinz[i][0], spinz[i][1], spinz[i][2]);
    COMPARE(cg, correct_coefficient[i])
      << '\n' // Using double quotes here produces an error(?!)
      << "J1: " << spin[i][0] << " Jz1: " << spinz[i][0] << "\n"
      << "J2: " << spin[i][1] << " Jz2: " << spinz[i][1] << "\n"
      << "J3: " << spin[i][2] << " Jz3: " << spinz[i][2] << "\n"
      << "CG: " << cg
      << " Correct: " << correct_coefficient[i];
  }
}
