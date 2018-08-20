/*
 *
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "../include/smash/tabulation.h"

using namespace smash;

TEST(constant) {
  // tabulate a constant function
  const Tabulation tab(0., 10., 10, [](double) { return 1.; });
  // check closest-point values
  FUZZY_COMPARE(tab.get_value_step(-1.), 0.);
  FUZZY_COMPARE(tab.get_value_step(0.), 1.);
  FUZZY_COMPARE(tab.get_value_step(5.), 1.);
  FUZZY_COMPARE(tab.get_value_step(7.5), 1.);
  FUZZY_COMPARE(tab.get_value_step(10.), 1.);
  FUZZY_COMPARE(tab.get_value_step(20.), 1.);
  // check interpolated values
  FUZZY_COMPARE(tab.get_value_linear(-1.), 0.);
  FUZZY_COMPARE(tab.get_value_linear(0.), 1.);
  FUZZY_COMPARE(tab.get_value_linear(5.), 1.);
  FUZZY_COMPARE(tab.get_value_linear(7.5), 1.);
  FUZZY_COMPARE(tab.get_value_linear(10.), 1.);
  // check extrapolated values
  FUZZY_COMPARE(tab.get_value_linear(20.), 1.);
}

TEST(linear) {
  // tabulate a linear function
  const Tabulation tab(0., 10., 10, [](double x) { return x; });
  // check closest-point values
  FUZZY_COMPARE(tab.get_value_step(-1.), 0.);
  FUZZY_COMPARE(tab.get_value_step(0.), 0.);
  FUZZY_COMPARE(tab.get_value_step(1.), 1.);
  FUZZY_COMPARE(tab.get_value_step(5.), 5.);
  FUZZY_COMPARE(tab.get_value_step(7.5), 8.);
  FUZZY_COMPARE(tab.get_value_step(8.4), 8.);
  FUZZY_COMPARE(tab.get_value_step(10.), 10.);
  FUZZY_COMPARE(tab.get_value_step(11.), 10.);
  // check interpolated values
  FUZZY_COMPARE(tab.get_value_linear(0.), 0.);
  FUZZY_COMPARE(tab.get_value_linear(1.), 1.);
  FUZZY_COMPARE(tab.get_value_linear(5.), 5.);
  FUZZY_COMPARE(tab.get_value_linear(7.5), 7.5);
  FUZZY_COMPARE(tab.get_value_linear(8.4), 8.4);
  FUZZY_COMPARE(tab.get_value_linear(10.), 10.);
  // check extrapolated values
  FUZZY_COMPARE(tab.get_value_linear(11.), 11.);
  FUZZY_COMPARE(tab.get_value_linear(20.), 20.);
}

TEST(quadratic) {
  // tabulate a quadratic function
  const Tabulation tab(-2., 4., 20, [](double x) { return x * x; });
  const double error = 1E-5;
  UnitTest::setFuzzyness<double>(2);
  // check closest-point values
  FUZZY_COMPARE(tab.get_value_step(-3.), 0.);
  FUZZY_COMPARE(tab.get_value_step(-2.), 4.);
  FUZZY_COMPARE(tab.get_value_step(-1.), 1.);
  COMPARE_ABSOLUTE_ERROR(tab.get_value_step(0.), 0., error);
  COMPARE_ABSOLUTE_ERROR(tab.get_value_step(0.5), 0.36, error);
  FUZZY_COMPARE(tab.get_value_step(1.), 1.);
  FUZZY_COMPARE(tab.get_value_step(1.2), 1.44);
  FUZZY_COMPARE(tab.get_value_step(2.), 4.);
  FUZZY_COMPARE(tab.get_value_step(3.), 4.);
  // check interpolated values
  FUZZY_COMPARE(tab.get_value_linear(-3.), 0.);
  FUZZY_COMPARE(tab.get_value_linear(-2.), 4.);
  FUZZY_COMPARE(tab.get_value_linear(-1.), 1.);
  COMPARE_ABSOLUTE_ERROR(tab.get_value_linear(0.), 0., error);
  COMPARE_ABSOLUTE_ERROR(tab.get_value_linear(0.5), 0.26, error);
  FUZZY_COMPARE(tab.get_value_linear(1.), 1.);
  FUZZY_COMPARE(tab.get_value_linear(1.2), 1.44);
  FUZZY_COMPARE(tab.get_value_linear(2.), 4.);
  // check extrapolated values
  COMPARE_ABSOLUTE_ERROR(tab.get_value_linear(3.), 7.8, error);
}
