/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "tests/unittest.h"

#include "include/experiment.h"

TEST(create_box) {
  VERIFY(!!ExperimentBase::create("Box", 0));
}

TEST(create_collider) {
  VERIFY(!!ExperimentBase::create("Collider", 0));
}

TEST_CATCH(create_invalid, ExperimentBase::InvalidModusRequest) {
  ExperimentBase::create("Invalid", 0);
}
