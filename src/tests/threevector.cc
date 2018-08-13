/*
 *
 *    Copyright (c) 2013-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "../include/smash/threevector.h"

using namespace smash;

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
  FUZZY_COMPARE(C.sqr(), 0.123456 * 0.123456 + 2 * 1.9e7 * 1.9e7);
  FUZZY_COMPARE(C.abs(), sqrt(0.123456 * 0.123456 + 2 * 1.9e7 * 1.9e7));
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

TEST(get_angles) {
  // zero vector
  ThreeVector Z(0., 0., 0.);
  FUZZY_COMPARE(Z.get_phi(), 0.);
  FUZZY_COMPARE(Z.get_theta(), 0.);
  // unit vectors
  ThreeVector A(1., 0., 0.);
  FUZZY_COMPARE(A.get_phi(), 0.);
  FUZZY_COMPARE(A.get_theta(), M_PI / 2.);
  ThreeVector B(0., 1., 0.);
  FUZZY_COMPARE(B.get_phi(), M_PI / 2.);
  FUZZY_COMPARE(B.get_theta(), M_PI / 2.);
  ThreeVector C(-1., 0., 0.);
  FUZZY_COMPARE(C.get_phi(), M_PI);
  FUZZY_COMPARE(C.get_theta(), M_PI / 2.);
  ThreeVector D(0., -1., 0.);
  FUZZY_COMPARE(D.get_phi(), -M_PI / 2.);
  FUZZY_COMPARE(D.get_theta(), M_PI / 2.);
  ThreeVector E(0., 0., 1.);
  FUZZY_COMPARE(E.get_phi(), 0.);
  FUZZY_COMPARE(E.get_theta(), 0.);
  ThreeVector F(0., 0., -1.);
  FUZZY_COMPARE(F.get_phi(), 0.);
  FUZZY_COMPARE(F.get_theta(), M_PI);
  // some other cases
  ThreeVector G(1., 1., 0.);
  FUZZY_COMPARE(G.get_phi(), M_PI / 4.);
  FUZZY_COMPARE(G.get_theta(), M_PI / 2.);
  ThreeVector H(1., -1., 0.);
  FUZZY_COMPARE(H.get_phi(), -M_PI / 4.);
  FUZZY_COMPARE(H.get_theta(), M_PI / 2.);
  ThreeVector I(1., 0., 1.);
  FUZZY_COMPARE(I.get_phi(), 0.);
  FUZZY_COMPARE(I.get_theta(), M_PI / 4.);
  ThreeVector J(1., 0., -1.);
  FUZZY_COMPARE(J.get_phi(), 0.);
  FUZZY_COMPARE(J.get_theta(), M_PI * 3. / 4.);
}

TEST(rotations) {
  // rotate around y
  ThreeVector A(1., 0., 0.);
  A.rotate_around_y(M_PI / 2);
  COMPARE_ABSOLUTE_ERROR(A.x1(), 0., 1.e-15);
  COMPARE_ABSOLUTE_ERROR(A.x2(), 0., 1.e-15);
  COMPARE_ABSOLUTE_ERROR(A.x3(), -1., 1.e-15);
  A.rotate_around_y(M_PI / 2);
  COMPARE_ABSOLUTE_ERROR(A.x1(), -1., 1.e-15);
  COMPARE_ABSOLUTE_ERROR(A.x2(), 0., 1.e-15);
  COMPARE_ABSOLUTE_ERROR(A.x3(), 0., 1.e-15);
  A.rotate_around_y(M_PI / 2);
  COMPARE_ABSOLUTE_ERROR(A.x1(), 0., 1.e-15);
  COMPARE_ABSOLUTE_ERROR(A.x2(), 0., 1.e-15);
  COMPARE_ABSOLUTE_ERROR(A.x3(), 1., 1.e-15);
  A.rotate_around_y(M_PI / 2);
  COMPARE_ABSOLUTE_ERROR(A.x1(), 1., 1.e-15);
  COMPARE_ABSOLUTE_ERROR(A.x2(), 0., 1.e-15);
  COMPARE_ABSOLUTE_ERROR(A.x3(), 0., 1.e-15);
  // rotate around z
  ThreeVector B(1., 0., 0.);
  B.rotate_around_z(M_PI / 2);
  COMPARE_ABSOLUTE_ERROR(B.x1(), 0., 1.e-15);
  COMPARE_ABSOLUTE_ERROR(B.x2(), 1., 1.e-15);
  COMPARE_ABSOLUTE_ERROR(B.x3(), 0., 1.e-15);
  B.rotate_around_z(M_PI / 2);
  COMPARE_ABSOLUTE_ERROR(B.x1(), -1., 1.e-15);
  COMPARE_ABSOLUTE_ERROR(B.x2(), 0., 1.e-15);
  COMPARE_ABSOLUTE_ERROR(B.x3(), 0., 1.e-15);
  B.rotate_around_z(M_PI / 2);
  COMPARE_ABSOLUTE_ERROR(B.x1(), 0., 1.e-15);
  COMPARE_ABSOLUTE_ERROR(B.x2(), -1., 1.e-15);
  COMPARE_ABSOLUTE_ERROR(B.x3(), 0., 1.e-15);
  B.rotate_around_z(M_PI / 2);
  COMPARE_ABSOLUTE_ERROR(B.x1(), 1., 1.e-15);
  COMPARE_ABSOLUTE_ERROR(B.x2(), 0., 1.e-15);
  COMPARE_ABSOLUTE_ERROR(B.x3(), 0., 1.e-15);
  // rotate_z_axis_to
  ThreeVector C(0., 0., 1.);
  ThreeVector R(1., 1., 1.);
  C.rotate_z_axis_to(R);
  COMPARE_ABSOLUTE_ERROR(C.x1(), R.x1() / R.abs(), 1.e-15);
  COMPARE_ABSOLUTE_ERROR(C.x2(), R.x2() / R.abs(), 1.e-15);
  COMPARE_ABSOLUTE_ERROR(C.x3(), R.x3() / R.abs(), 1.e-15);
}

TEST(compares) {
  ThreeVector a = {};
  const auto b = a;
  VERIFY(a == b);
  VERIFY(!(a != b));
  a[2] += 1;
  VERIFY(!(a == b));
  VERIFY(a != b);
  a = b;
  a[1] += 1;
  VERIFY(!(a == b));
  VERIFY(a != b);
  a = b;
  a[0] += 1;
  VERIFY(!(a == b));
  VERIFY(a != b);
}
