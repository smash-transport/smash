/*
 *
 *    Copyright (c) 2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <vir/test.h>  // This include has to be first

#include "setup.h"

#include "../include/smash/constants.h"
#include "../include/smash/formfactors.h"

using namespace smash;

TEST(blatt_weisskopf) {
  vir::test::setFuzzyness<double>(2);
  const double x = 0.5 * hbarc;
  COMPARE(blatt_weisskopf_sqr(x, 0), 1.);
  const double x2 = 0.5 * 0.5;
  COMPARE(blatt_weisskopf_sqr(x, 1), x2 / 1.25);
  const double x4 = x2 * x2;
  COMPARE(blatt_weisskopf_sqr(x, 2), x4 / 9.8125);
  const double x6 = x4 * x2;
  COMPARE(blatt_weisskopf_sqr(x, 3), x6 / 236.640625);
  const double x8 = x4 * x4;
  COMPARE(blatt_weisskopf_sqr(x, 4), x8 / 11427.34765625);
  const double x10 = x8 * x2;
  FUZZY_COMPARE(blatt_weisskopf_sqr(x, 5), x10 / 918229.981445313);
}
