/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "unittest.h"
#include "../include/fourvector.h"
#include "../include/angles.h"

using namespace Smash;

constexpr double accuracy = 4e-9;
Angles dir;
auto cos_like = Random::make_uniform_distribution(-1.0, +1.0);

ThreeVector random_velocity();
ThreeVector random_velocity() {
  dir.distribute_isotropically();
  double beta = Random::canonical();
  return dir.threevec()*beta;
}

// here, we boost a velocity vector with itself.
TEST(self_boost) {
  for(int i = 0; i < 1000000; i++) {
    ThreeVector velocity = random_velocity();
    // u_mu is a real four-vector.
    double gamma = 1./sqrt(1. - velocity.sqr());
    FourVector u_mu = FourVector (gamma, velocity*gamma);
    FourVector boosted = u_mu.LorentzBoost(velocity);
    COMPARE_ABSOLUTE_ERROR(boosted.x0(), 1.0, accuracy) << " at loop " << i;
    COMPARE_ABSOLUTE_ERROR(boosted.x1(), 0.0, accuracy) << " at loop " << i;
    COMPARE_ABSOLUTE_ERROR(boosted.x2(), 0.0, accuracy) << " at loop " << i;
    COMPARE_ABSOLUTE_ERROR(boosted.x3(), 0.0, accuracy) << " at loop " << i;
  }
}

// try to keep the invariants invariant
// 1. "length" of a four-vector
TEST(keep_invariant_length) {
  for(int i = 0; i < 1000; i++) {
    ThreeVector velocity = random_velocity();
    for(int j = 0; j < 1000; j++) {
      FourVector a( cos_like()
                  , cos_like()
                  , cos_like()
                  , cos_like());
      FourVector A = a.LorentzBoost(velocity);
      COMPARE_RELATIVE_ERROR(a.sqr(), A.sqr(), accuracy) << " at loop " << i
                                                                << "*" << j;
    }
  }
}

// 2. scalar product between two four-vectors
TEST(keep_invariant_angle) {
  for(int i = 0; i < 1000; i++) {
    ThreeVector velocity = random_velocity();
    for(int j = 0; j < 1000; j++) {
      FourVector a( cos_like()
                  , cos_like()
                  , cos_like()
                  , cos_like());
      FourVector b( cos_like()
                  , cos_like()
                  , cos_like()
                  , cos_like());
      FourVector A = a.LorentzBoost(velocity);
      FourVector B = b.LorentzBoost(velocity);
      COMPARE_RELATIVE_ERROR(a.Dot(b), A.Dot(B), accuracy) << " at loop " << i
                                                                << "*" << j;
    }
  }
}

// Lorentz transformation and back should get the same vector:
TEST(back_and_forth) {
  for(int i = 0; i < 1000; i++) {
    ThreeVector velocity = random_velocity();
    for(int j = 0; j < 1000; j++) {
      FourVector a( cos_like()
                  , cos_like()
                  , cos_like()
                  , cos_like());
      FourVector forward  = a.LorentzBoost(velocity);
      FourVector backward = forward.LorentzBoost(-velocity);
      COMPARE_RELATIVE_ERROR(backward.x0(), a.x0(), accuracy) << " at loop "
                                                        << i << "*" << j;
      COMPARE_RELATIVE_ERROR(backward.x1(), a.x1(), accuracy) << " at loop "
                                                        << i << "*" << j;
      COMPARE_RELATIVE_ERROR(backward.x2(), a.x2(), accuracy) << " at loop "
                                                        << i << "*" << j;
      COMPARE_RELATIVE_ERROR(backward.x3(), a.x3(), accuracy) << " at loop "
                                                        << i << "*" << j;
    }
  }
}
