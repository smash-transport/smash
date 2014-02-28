/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "tests/unittest.h"
#include "include/fourvector.h"

using namespace Smash;

FourVector p(4, 1, 2, 3);

TEST(boost_to_rest) {
  FourVector v(1.0, 0.25, 0.5, 0.75);
  FourVector p_boosted = p.LorentzBoost(v);
  FUZZY_COMPARE(p_boosted.Dot(),p.Dot());
  // space-like components should have been boosted to zero:
  FUZZY_COMPARE(p_boosted.x1(),0.0);
  FUZZY_COMPARE(p_boosted.x2(),0.0);
  FUZZY_COMPARE(p_boosted.x3(),0.0);
  UnitTest::setFuzzyness<double>(2);
  // boost back and compare with original vector:
  FourVector v_neg(1, -v.x1(), -v.x2(), -v.x3());
  FourVector p_returned = p_boosted.LorentzBoost(v_neg);
  FUZZY_COMPARE(p.x1(), p_returned.x1());
  FUZZY_COMPARE(p.x2(), p_returned.x2());
  FUZZY_COMPARE(p.x3(), p_returned.x3());
}
