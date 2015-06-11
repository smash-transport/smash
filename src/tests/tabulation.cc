/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "../include/tabulation.h" 

using namespace Smash;

TEST(constant) {
  // tabulate a constant function
  Tabulation tab(0., 10., 10, [](double) { return 1.; });
  // check closest-point values
  FUZZY_COMPARE(tab.get_value_step(-1.), 0.f);
  FUZZY_COMPARE(tab.get_value_step( 0.), 1.f);
  FUZZY_COMPARE(tab.get_value_step( 5.), 1.f);
  FUZZY_COMPARE(tab.get_value_step(7.5), 1.f);
  FUZZY_COMPARE(tab.get_value_step(10.), 1.f);
  FUZZY_COMPARE(tab.get_value_step(20.), 1.f);
  // check interpolated values
  FUZZY_COMPARE(tab.get_value_linear(-1.), 0.f);
  FUZZY_COMPARE(tab.get_value_linear( 0.), 1.f);
  FUZZY_COMPARE(tab.get_value_linear( 5.), 1.f);
  FUZZY_COMPARE(tab.get_value_linear(7.5), 1.f);
  FUZZY_COMPARE(tab.get_value_linear(10.), 1.f);
  FUZZY_COMPARE(tab.get_value_linear(20.), 1.f);
}

TEST(linear) {
  // tabulate a linear function
  Tabulation tab(0., 10., 10, [](double x) { return x; });
  // check closest-point values
  FUZZY_COMPARE(tab.get_value_step(-1.),  0.f);
  FUZZY_COMPARE(tab.get_value_step( 0.),  0.f);
  FUZZY_COMPARE(tab.get_value_step( 1.),  1.f);
  FUZZY_COMPARE(tab.get_value_step( 5.),  5.f);
  FUZZY_COMPARE(tab.get_value_step(7.5),  8.f);
  FUZZY_COMPARE(tab.get_value_step(8.4),  8.f);
  FUZZY_COMPARE(tab.get_value_step(10.), 10.f);
  FUZZY_COMPARE(tab.get_value_step(11.), 10.f);
  // check interpolated values
  FUZZY_COMPARE(tab.get_value_linear( 0.),  0.f);
  FUZZY_COMPARE(tab.get_value_linear( 1.),  1.f);
  FUZZY_COMPARE(tab.get_value_linear( 5.),  5.f);
  FUZZY_COMPARE(tab.get_value_linear(7.5), 7.5f);
  FUZZY_COMPARE(tab.get_value_linear(8.4), 8.4f);
  FUZZY_COMPARE(tab.get_value_linear(10.), 10.f);
  FUZZY_COMPARE(tab.get_value_linear(11.), 10.f);
}

TEST(quadratic) {
  // tabulate a quadratic function
  Tabulation tab(-2., 4., 20, [](double x) { return x*x; });
  // check closest-point values
  FUZZY_COMPARE(tab.get_value_step(-3.), 0.f);
  FUZZY_COMPARE(tab.get_value_step(-2.), 4.f);
  FUZZY_COMPARE(tab.get_value_step(-1.), 1.f);
  FUZZY_COMPARE(tab.get_value_step( 0.), 0.f);
  COMPARE_ABSOLUTE_ERROR(tab.get_value_step(0.5), 0.36f, 1E-6f); 
  FUZZY_COMPARE(tab.get_value_step( 1.), 1.f);
  FUZZY_COMPARE(tab.get_value_step(1.2), 1.44f);
  FUZZY_COMPARE(tab.get_value_step( 2.), 4.f);
  FUZZY_COMPARE(tab.get_value_step( 3.), 4.f);
  // check interpolated values
  FUZZY_COMPARE(tab.get_value_linear(-3.), 0.f);
  FUZZY_COMPARE(tab.get_value_linear(-2.), 4.f);
  FUZZY_COMPARE(tab.get_value_linear(-1.), 1.f);
  FUZZY_COMPARE(tab.get_value_linear( 0.), 0.f);
  COMPARE_ABSOLUTE_ERROR(tab.get_value_linear(0.5), 0.26f, 1E-6f);
  FUZZY_COMPARE(tab.get_value_linear( 1.), 1.f);
  FUZZY_COMPARE(tab.get_value_linear(1.2), 1.44f);
  FUZZY_COMPARE(tab.get_value_linear( 2.), 4.f);
  FUZZY_COMPARE(tab.get_value_linear( 3.), 4.f);
}