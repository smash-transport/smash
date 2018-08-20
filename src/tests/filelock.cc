/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "setup.h"

#include "../include/smash/filelock.h"

using namespace smash;

static const bf::path testoutputpath = bf::absolute(SMASH_TEST_OUTPUT_PATH);

TEST(directory_is_created) {
  bf::create_directories(testoutputpath);
  VERIFY(bf::exists(testoutputpath));
}

TEST(lock_is_created) {
  const bf::path lockpath = testoutputpath / "lock";
  {
    FileLock lock(lockpath);
    VERIFY(lock.acquire());
    VERIFY(bf::exists(lockpath));
  }
  VERIFY(!bf::exists(lockpath));
}

TEST(lock_is_denied) {
  const bf::path lockpath = testoutputpath / "lock";
  {
    FileLock lock1(lockpath);
    FileLock lock2(lockpath);
    VERIFY(lock1.acquire());
    VERIFY(!lock2.acquire());
    VERIFY(bf::exists(lockpath));
  }
  VERIFY(!bf::exists(lockpath));
}
