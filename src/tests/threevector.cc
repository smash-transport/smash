/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "../include/threevector.h"
#include "unittest.h"

using namespace Smash;

TEST(assign) {
  ThreeVector A;
  COMPARE(A.x1(), 0.0);
  COMPARE(A.x2(), 0.0);
  COMPARE(A.x3(), 0.0);
  // exact assignment with binary-representable numbers
  ThreeVector B(0.25, -2.5, 3.75);
  COMPARE(B.x1(), 0.25);
  COMPARE(B.x2(), -2.5);
  COMPARE(B.x3(), 3.75);
  // those numbers cannot be represented exactly in binary:
  ThreeVector C(0.33, -2.1, 3.85);
  FUZZY_COMPARE(B.x1(), 0.25);
  FUZZY_COMPARE(B.x2(), -2.5);
  FUZZY_COMPARE(B.x3(), 3.75);
}

TEST(set) {
  ThreeVector A;
  A.set_x1(0.4);
  FUZZY_COMPARE(A.x1(), 0.4);
  // make sure the other components are unaffected:
  FUZZY_COMPARE(A.x2(), 0.0);
  FUZZY_COMPARE(A.x3(), 0.0);
  A.set_x2(0.6);
  FUZZY_COMPARE(A.x1(), 0.4);
  FUZZY_COMPARE(A.x2(), 0.6);
  FUZZY_COMPARE(A.x3(), 0.0);
  A.set_x3(0.8);
  FUZZY_COMPARE(A.x1(), 0.4);
  FUZZY_COMPARE(A.x2(), 0.6);
  FUZZY_COMPARE(A.x3(), 0.8);
}

TEST(sqr_abs) {
  ThreeVector A(3.0, 4.0, 0.0);
  FUZZY_COMPARE(A.sqr(), 25.0);
  FUZZY_COMPARE(A.abs(), 5.0);
  ThreeVector B(0.3, 0.0, 0.4);
  FUZZY_COMPARE(B.sqr(), 0.25);
  FUZZY_COMPARE(B.abs(), 0.5);
  ThreeVector C(0.123456, 1.9e7, -1.9e7);
  FUZZY_COMPARE(C.sqr(), C.abs() * C.abs());
  FUZZY_COMPARE(C.sqr(),      0.123456*0.123456 + 2 * 1.9e7 * 1.9e7);
  FUZZY_COMPARE(C.abs(), sqrt(0.123456*0.123456 + 2 * 1.9e7 * 1.9e7));
}

TEST(arithmetic) {
  ThreeVector A(3.0, 4.0, 0.0);
  ThreeVector B(0.3, 0.0, 0.4);
  ThreeVector C = A - B;
  FUZZY_COMPARE(C.x1(), 2.7);
  FUZZY_COMPARE(C.x2(), 4.0);
  FUZZY_COMPARE(C.x3(), -.4);
  C += B;
  FUZZY_COMPARE(C.x1(), A.x1());
  FUZZY_COMPARE(C.x2(), A.x2());
  FUZZY_COMPARE(C.x3(), A.x3());
  C -= A;
  FUZZY_COMPARE(C.x1(), 0.0);
  FUZZY_COMPARE(C.x2(), 0.0);
  FUZZY_COMPARE(C.x3(), 0.0);
  ThreeVector D = A;
  D *= 2;
  FUZZY_COMPARE(D.x1(), 2.0 * A.x1());
  FUZZY_COMPARE(D.x2(), 2.0 * A.x2());
  FUZZY_COMPARE(D.x3(), 2.0 * A.x3());
  D /= 2;
  FUZZY_COMPARE(D.x1(), A.x1());
  FUZZY_COMPARE(D.x2(), A.x2());
  FUZZY_COMPARE(D.x3(), A.x3());

}
