#include "unittest.h"

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include "../include/kinematics.h"
#include "../include/photoncrosssections.h"
#include "../include/tabulationnd.h"
#include "setup.h"

TEST(compile)
{
  PhotonCrossSection<ComputationMethod::Lookup> xs_object;
  
}