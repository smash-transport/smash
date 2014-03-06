/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "tests/unittest.h"
#include "include/fourvector.h"
#include "include/angles.h"

using namespace Smash;

static const double kEPS = 1e-10;
Angles dir;

FourVector random_velocity();
FourVector random_velocity() {
  dir.distribute_isotropically();
  double beta = rng.uniform(0.0,1.0);
  // velocity-"vector" is not normalized (that's how LorentzBoost
  // works):
  return FourVector(1.0, beta*dir.x(), beta*dir.y(), beta*dir.z());
}

// here, we boost a velocity vector with itself.
TEST(self_boost) {
  for(int i = 0; i < 1000000; i++) {
    FourVector velocity = random_velocity();
    // u_mu is a real four-vector.
    double gamma = 1.0/sqrt(velocity.Dot());
    FourVector u_mu = velocity*gamma;
    FourVector boosted = u_mu.LorentzBoost(velocity);
    FUZZY_COMPARE(boosted.x0(),1.0) << " at loop " << i;
    VERIFY(std::abs(boosted.x1() - 0.0) < kEPS) << " at loop " << i;
    VERIFY(std::abs(boosted.x2() - 0.0) < kEPS) << " at loop " << i;
    VERIFY(std::abs(boosted.x3() - 0.0) < kEPS) << " at loop " << i;
  }
}

// try to keep the invariants invariant
// 1. "length" of a four-vector
TEST(keep_invariant_length) {
  UnitTest::setFuzzyness<double>(32);
  for(int i = 0; i < 1000; i++) {
    FourVector velocity = random_velocity();
    for(int j = 0; j < 1000; j++) {
      FourVector a( rng.cos_like()
                  , rng.cos_like()
                  , rng.cos_like()
                  , rng.cos_like());
      FourVector A = a.LorentzBoost(velocity);
      FUZZY_COMPARE(a.Dot(), A.Dot());
    }
  }
}

// 2. scalar product between two four-vectors
TEST(keep_invariant_angle) {
  UnitTest::setFuzzyness<double>(64);
  for(int i = 0; i < 1000; i++) {
    FourVector velocity = random_velocity();
    for(int j = 0; j < 1000; j++) {
      FourVector a( rng.cos_like()
                  , rng.cos_like()
                  , rng.cos_like()
                  , rng.cos_like());
      FourVector b( rng.cos_like()
                  , rng.cos_like()
                  , rng.cos_like()
                  , rng.cos_like());
      FourVector A = a.LorentzBoost(velocity);
      FourVector B = b.LorentzBoost(velocity);
      FUZZY_COMPARE(a.Dot(b), A.Dot(B));
    }
  }
}

// Lorentz transformation and back should get the same vector:
TEST(back_and_forth) {
  UnitTest::setFuzzyness<double>(64);
  for(int i = 0; i < 1000; i++) {
    FourVector velocity = random_velocity();
    FourVector back(1, -velocity.x1(), -velocity.x2(), -velocity.x3());
    for(int j = 0; j < 1000; j++) {
      FourVector a( rng.cos_like()
                  , rng.cos_like()
                  , rng.cos_like()
                  , rng.cos_like());
      FourVector forward  = a.LorentzBoost(velocity);
      FourVector backward = forward.LorentzBoost(back);
      FUZZY_COMPARE(backward.x0(),a.x0());
      FUZZY_COMPARE(backward.x1(),a.x1());
      FUZZY_COMPARE(backward.x2(),a.x2());
      FUZZY_COMPARE(backward.x3(),a.x3());
    }
  }
}
