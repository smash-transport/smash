/*
 *
 *    Copyright (C) 2009-2014 Matthias Kretz <kretz@kde.org>
 *    Copyright (c) 2014,2017,2019,2023
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_TESTS_ASSERT_H_
#define SRC_TESTS_ASSERT_H_

/* When deprecating C functions with the  DEPERCATE_C_FNS option the assert test
 * has to be disabled, since due to the -include option used, the system assert
 * is included first. You still get warnings about it, though. */
#ifndef DONT_ASSERT_TEST
#ifdef assert
#error \
    "The system assert.h (or some other header that defines assert) was included before our own assert.h"
#endif
#endif  // DONT_ASSERT_TEST

namespace UnitTest {
void unittest_assert(bool cond, const char *code, const char *file, int line);
}  // namespace UnitTest
#define assert(cond) UnitTest::unittest_assert(cond, #cond, __FILE__, __LINE__);

#endif  // SRC_TESTS_ASSERT_H_
