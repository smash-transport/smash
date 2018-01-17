/*  TODO(mkretz): license?

    Copyright (C) 2009-2014 Matthias Kretz <kretz@kde.org>

    Vc is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    Vc is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with Vc.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SRC_TESTS_ASSERT_H_
#define SRC_TESTS_ASSERT_H_

#ifdef assert
#error \
    "The system assert.h (or some other header that defines assert) was included before our own assert.h"
#endif

namespace UnitTest {
void unittest_assert(bool cond, const char *code, const char *file, int line);
}  // namespace UnitTest
#define assert(cond) UnitTest::unittest_assert(cond, #cond, __FILE__, __LINE__);

#endif  // SRC_TESTS_ASSERT_H_
