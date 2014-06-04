/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include <boost/filesystem.hpp>
#include "../include/outputinterface.h"
#include "../include/oscaroutput.h"
#include "../include/particles.h"

using namespace Smash;

const bf::path testoutputpath = bf::absolute(SMASH_TEST_OUTPUT_PATH);

TEST(directory_is_created) {
  bf::create_directory(testoutputpath);
  VERIFY(bf::exists(testoutputpath));
}

TEST(file_is_created) {
  OscarOutput *oscout = new OscarOutput(testoutputpath);
  VERIFY(bf::exists(testoutputpath / "collision.dat"));
}
