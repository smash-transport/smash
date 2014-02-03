/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "../include/Angles.h"

int main() {
  static const double EPS = 1e-7;
  static const int NUMBER_OF_TRIES = 1000000;
  static const int MIN_M = -6;
  static const int MAX_M = 12;

  Angles dir;

  // first, usual set:
  dir.set_phi(.5);
  if (fabs(dir.phi()-.5) >= EPS)
    return 1;
  dir.set_costheta(.3);
  if (fabs(dir.costheta()-.3) >= EPS)
    return 1;
  dir.set_theta(.3);
  if (fabs(dir.theta()-.3) >= EPS)
    return 1;

  // second: check accessors and the relations between them:
  for (int c = 0; c < NUMBER_OF_TRIES; c++) {
    dir.distribute_isotropically();
    // sintheta**2 + costheta**2 = 1
    double theta_one = dir.sintheta()*dir.sintheta()
                     + dir.costheta()*dir.costheta();
    if (fabs(theta_one - 1.0) >= EPS)
      return 2;
    // x**2 + y**2 + z**2 = 1
    double xyz_one = dir.x()*dir.x()
                   + dir.y()*dir.y()
                   + dir.z()*dir.z();
    if (fabs(xyz_one - 1.0) >= EPS)
      return 3;

    // compare cos(theta) and costheta:
    if (fabs(cos(dir.theta()) - dir.costheta()) >= EPS)
      return 4;
    if (fabs(dir.theta() - acos(dir.costheta())) >= EPS)
      return 5;
  }

  // third: Unusual set
  for (int m = MIN_M; m < MAX_M; m++) {
    // set phi outside [0..2pi]
    dir.set_phi(2.0 * M_PI * m + .5);
    if (fabs(dir.phi() - .5) >= EPS)
      return 6;
    // set theta in [2*n*pi .. (2*n+1)*pi]
    dir.set_theta(2.0 * M_PI * m + .7);
    if (fabs(dir.theta() - .7) >= EPS) {
      return 6;
    }
    // set theta in [(2*n-1)*pi .. 2*n*pi]
    // dir.set_theta(2.0 * M_PI * m - .7);
    // if (fabs(dir.theta() - .7) >= EPS)
    //  return 6;
  }
  for (double newcos = -8.0; newcos < 8.0; newcos += .2) {
    // Did I catch an exception?
    bool caught = false;
    // this should not work:
    try {
      dir.set_costheta(newcos);
    } catch(...) {
      caught = true;
    }
    // if I don't catch an expression although newcos is outside range:
    if (!caught && (newcos < -1 || newcos > 1))
      return 7;
    // if I catch an exception although newcos is inside rande:
    if (caught && (newcos >= -1 && newcos <= 1))
      return 8;
  }

  // now, set and check if the other field(s) are effected.
  double old_phi = dir.phi();
  dir.set_costheta(.2);
  if (fabs(old_phi - dir.phi()) >= EPS)
    return 9;
  dir.set_theta(.2);
  if (fabs(old_phi - dir.phi()) >= EPS)
    return 9;
  double old_costheta = dir.costheta();
  double old_z = dir.z();
  dir.set_phi(.4);
  if (fabs(old_costheta - dir.costheta()) >= EPS)
    return 10;
  if (fabs(old_z - dir.z()) >= EPS)
    return 11;

  // and now, the climax: add_to_theta.
  for (double current_phi = 0.0; current_phi < 2 * M_PI;
                                 current_phi += M_PI / 180.0) {
    dir.set_phi(current_phi);
    dir.set_theta(0.0);
    bool sign = dir.add_to_theta(M_PI / 2.0);
    // phi shouldn't have changed:
    if (fabs(current_phi - dir.phi()) >= EPS)
      return 33;
    if (sign)
      return 33;
    // theta, though, should have changed.
    if (fabs(dir.theta() - M_PI / 2.0) >= EPS)
      return 34;
    // 90+120 degrees: phi changed, theta is at 150 degrees
    sign = dir.add_to_theta(2.0 * M_PI / 3.0);
    if (current_phi < M_PI) {
      if (fabs(current_phi + M_PI - dir.phi()) >= EPS)
        return 35;
    } else {
      if (fabs(current_phi - M_PI - dir.phi()) >= EPS)
        return 35;
    }
    if (!sign)
      return 35;
    if (fabs(dir.theta() - 5.0 * M_PI / 6.0) >= EPS)
      return 36;
    // go back over the pole: +60 degrees
    sign = dir.add_to_theta(M_PI / 3.0);
    // phi is back to original;
    if (fabs(current_phi - dir.phi()) >= EPS)
      return 37;
    if (!sign)
      return 37;
    // theta is the same as before (150+60 = 210 equiv 150)
    if (fabs(dir.theta() - 5.0 * M_PI / 6.0) >= EPS)
      return 38;
    // go over two poles: +120 + 120 (going over more than 180 degrees
    // in one step is forbidden; we'll check that later)
    sign = dir.add_to_theta(2.0 * M_PI / 3.0);
    // other direction because of pole-crossing:
    bool sign2 = dir.add_to_theta(2.0 * M_PI / 3.0, sign);
    // phi stays the same (two changes)
    if (fabs(current_phi - dir.phi()) >= EPS)
      return 39;
    if (sign2 && !sign)
      return 39;
    // theta is now 30 degrees:
    if (fabs(dir.theta() - M_PI / 6.0) >= EPS)
      return 40;

    // add_to_theta should bail out for too big additions:
    bool caught = false;
    // this should not work:
    try {
      dir.add_to_theta(M_PI * 1.1);
    } catch(...) {
      caught = true;
    }
    if (!caught)
      return 41;
  }

  // we can also wildly add things to theta and see if the boolean
  // really corresponds to a change in phi:
  for (int c = 0; c < NUMBER_OF_TRIES; c++) {
    double current_phi = dir.phi();
    if (dir.add_to_theta(M_PI * drand48())
        && fabs(current_phi - dir.phi()) < EPS)
      return 42;
  }

  // that's all, folks!
  return 0;
}
