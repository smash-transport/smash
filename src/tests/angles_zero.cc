/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "../include/smash/angles.h"

// tests if the angles add up to zero.

using namespace smash;

int main() {
  const int NUMBER = 100000000;
  double sumx = 0;
  double sumy = 0;
  double sumz = 0;
  for (int c = 0; c < NUMBER; c++) {
    Angles dir;
    dir.distribute_isotropically();
    sumx += dir.x();
    sumy += dir.y();
    sumz += dir.z();
  }
  printf("%g %g %g\n", sumx / (NUMBER + 0.0), sumy / (NUMBER + 0.0),
         sumz / (NUMBER + 0.0));
  // FIXME This test does not do anything
  return 0;
}
