#include "unittest.h"

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include "setup.h"

#include "../include/tabulationnd.h"


using namespace Smash;

//TEST(compile) { PhotonCrossSection<ComputationMethod::Lookup> xs_object; }

TEST(linear_3d) {
  const double x0 = 0.0, x1 = 4.0, y0 = 0.0, y1 = 4.0, z0 = 0.0, z1 = 4.0;
  const double dx = 0.1, dy = 0.1, dz = 0.1;
  TabulationND<3> tab{x0, x1, y0, y1, z0, z1, dx, dy, dz,
                     [](double x, double y, double z) { return x + y + z; }};
  
  FUZZY_COMPARE(tab.get_linear(0.0, 0.0, 0.0), 0.0);
  FUZZY_COMPARE(tab.get_linear(1.0, 0.0, 0.0), 1.0);
  FUZZY_COMPARE(tab.get_linear(0.0, 1.0, 0.0), 1.0);
  FUZZY_COMPARE(tab.get_linear(0.0, 0.0, 1.0), 1.0);
  FUZZY_COMPARE(tab.get_linear(1.0, 1.0, 1.0), 3.0);
  FUZZY_COMPARE(tab.get_linear(4.0, 1.0, 1.0), 6.0);
  FUZZY_COMPARE(tab.get_linear(4.0, 4.0, 4.0), 12.0);
  FUZZY_COMPARE(tab.get_linear(4.0, 0.0, 0.0), 4.0);
  FUZZY_COMPARE(tab.get_linear(0.0, 4.0, 0.0), 4.0);
  FUZZY_COMPARE(tab.get_linear(0.0, 0.0, 4.0), 4.0);




  COMPARE_RELATIVE_ERROR(tab.get_linear(0.05, 0.0, 0.0), 0.05, 0.01);
  COMPARE_RELATIVE_ERROR(tab.get_linear(0.05, 0.05, 0.05), 0.15, 0.01);
  
  
  
}
