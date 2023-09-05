/*
 *
 *    Copyright (c) 2023
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/rootsolver.h"
using namespace smash;
TEST(root_of_cosine) {
  std::function<double(double)> cosinus = [](double x) { return std::cos(x); };
  auto rootsolver = RootSolver1D(cosinus);
  const auto hopefully_pi_half = rootsolver.try_find_root(0.0, M_PI, 10000);
  VERIFY(hopefully_pi_half.has_value());
  COMPARE_ABSOLUTE_ERROR((*hopefully_pi_half), M_PI / 2.0, 1e-7);
}
