/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#include "../include/FourVector.h"

int main() {
  FourVector p(4, 1, 2, 3);
  FourVector v(1.0, 0.25, 0.5, 0.75);
  FourVector p_boosted = p.LorentzBoost(v);
  const double tolerance = 1.0e-6;

  /* Check that the scalar product remains invariant */
  if (fabs(p_boosted.Dot(p_boosted) - p.Dot(p)) > tolerance)
    return -1;

  /* Check that the boost gives stationary vector */
  if (p_boosted.x1() > tolerance || p_boosted.x2() > tolerance
      || p_boosted.x3() > tolerance)
    return -2;

  /* Check that the return boost returns the original vector */
  FourVector v_neg(1, -v.x1(), -v.x2(), -v.x3());
  FourVector p_returned = p_boosted.LorentzBoost(v_neg);

  if (fabs(p_returned.x0() - p.x0()) > tolerance
      || fabs(p_returned.x1() - p.x1()) > tolerance
      || fabs(p_returned.x2() - p.x2()) > tolerance
      || fabs(p_returned.x3() - p.x3()) > tolerance)
    return -3;

  return 0;
}
