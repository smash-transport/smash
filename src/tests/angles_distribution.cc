/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "../include/smash/angles.h"

// tests distributions for the angles.

using namespace smash;

int main() {
  // histograms in 1 degree steps:
  int phi_histo[360] = {0};
  int theta_histo[180] = {0};
  // ... and .01 steps in cosine.
  int cos_histo[200] = {0};
  const double dangle = M_PI / 180.0;
  const double dcosine = .01;

  const int NUMBER = 1000000;
  for (int c = 0; c < NUMBER; c++) {
    Angles dir;
    dir.distribute_isotropically();
    phi_histo[static_cast<int>(floor(dir.phi() / dangle))]++;
    theta_histo[static_cast<int>(floor(dir.theta() / dangle))]++;
    cos_histo[static_cast<int>(floor((dir.costheta() + 1.0) / dcosine))]++;

    // print some of the vectors for plotting and visual checking:
    if (c % 500 == 0)
      fprintf(stderr, "%8.5f %8.5f %8.5f\n", dir.x(), dir.y(), dir.z());
  }
  printf("#%7s %8s %8s\n", "angle", "dN/dphi", "dN/dtheta");

  for (int p = 0; p < 180; p++) {
    printf("%3d %8.6f %5.2f %8.5f %8.5f %8.5f\n", p, p * dangle,
           p * dcosine - 1.0, (phi_histo[p] + 0.0) / dangle / (NUMBER + 0.0),
           (cos_histo[p] + 0.0) / dcosine / (NUMBER + 0.0),
           (theta_histo[p] + 0.0) / dangle / (NUMBER + 0.0));
  }
  for (int p = 180; p < 200; p++) {
    printf("%3d %8.6f %5.2f %8.5f %8.5f\n", p, p * dangle, p * dcosine - 1.0,
           (phi_histo[p] + 0.0) / dangle / (NUMBER + 0.0),
           (cos_histo[p] + 0.0) / dcosine / (NUMBER + 0.0));
  }
  for (int p = 200; p < 360; p++) {
    printf("%3d %8.6f  0.00 %8.5f\n", p, p * dangle,
           (phi_histo[p] + 0.0) / dangle / (NUMBER + 0.0));
  }
  // FIXME this test does not do anything

  return 0;
}
