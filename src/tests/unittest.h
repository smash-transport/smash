/*
 *    Copyright (c) 2014-2018,2020-2022,2025
 *      SMASH Team
 */

/*!
 * \page doxypage_unit_testing
 * \tableofcontents
 *
 * See also: \ref unittest
 *
 * \section unittest_intro Introduction
 * (The “Introduction” text is adapted from the <a
 * href="https://en.wikipedia.org/wiki/Unit_testing">Wikipedia article on Unit
 * testing</a>, for a stronger focus on SMASH and C++.)
 *
 * Unit testing is a software testing method by which
 * individual units of source code
 * are tested to determine if they are fit for use.
 * View a unit as the smallest testable part of an application.
 * Thus, a unit is typically an entire class.
 *
 * Ideally, each test case is independent from the others. Substitutes such as
 * <a href="https://en.wikipedia.org/wiki/Method_stub">method stubs</a> and <a
 * href="https://en.wikipedia.org/wiki/Mock_object">mock objects</a> can be used
 * to assist testing a module in isolation. Unit tests are written and
 * run by software developers to ensure that code meets its design and behaves
 * as intended.
 *
 * \subsection unittest_intro_benefits Benefits
 * The goal of unit testing is to isolate each part of the program and show that
 * the individual parts are correct. A unit test provides a strict, written
 * contract that the piece of code must satisfy. As a result, it affords several
 * benefits.
 *
 * \subsubsection unittest_intro_benefits_find_problems Finds problems early
 * Unit testing finds problems early in the development cycle.
 *
 * In test-driven development, unit tests are created before the code itself is
 * written. When the tests pass, that code is considered complete. The same unit
 * tests are run against that function frequently as the larger code base is
 * developed either as the code is changed or via an automated process with the
 * build. If the unit tests fail, it is considered to be a bug either in the
 * changed code or the tests themselves. The unit tests then allow the location
 * of the fault or failure to be easily traced. Since the unit tests alert the
 * development team of the problem before handing the code off to testers or
 * clients, it is still early in the development process.
 *
 * \subsubsection unittest_intro_benefits_change Facilitates change
 * Unit testing allows the programmer to refactor code at a later date, and make
 * sure the module still works correctly (e.g., in regression testing). The
 * procedure is to write test cases for all functions and methods so that
 * whenever a change causes a fault, it can be quickly identified.
 *
 * Readily available unit tests make it easy for the programmer to check whether
 * a piece of code is still working properly.
 *
 * \subsubsection unittest_intro_benefits_integration Simplifies integration
 * Unit testing may reduce uncertainty in the units themselves and can be used
 * in a bottom-up testing style approach. By testing the parts of a program
 * first and then testing the sum of its parts, integration testing becomes much
 * easier.
 *
 * \subsubsection unittest_intro_benefits_design Design
 * When software is developed using a test-driven approach, the combination of
 * writing the unit test to specify the interface plus the refactoring
 * activities performed after the test is passing, may take the place of formal
 * design. Each unit test can be seen as a design element specifying classes,
 * methods, and observable behaviour.
 *
 ******************************************************************************
 *
 * \section unittest_smash SMASH Specific Hints and Rules
 * The typical unit test in SMASH is centered around one C++ class.
 * But many of the classes in SMASH rely on specific data to do any useful
 * operations.
 * The obvious candidates are
 * - \ref smash::ParticleData
 * - \ref smash::ParticleType
 * - \ref smash::Particles
 * - \ref smash::Configuration
 * - \ref smash::ProcessBranch
 *
 * Each of these are interfaces to data that most classes in SMASH read,
 * modify, or create. For example, consider testing smash::DecayAction. The
 * class is created with a const-ref to a smash::ParticleData object. This
 * class in turn requires a smash::ParticleType object for its constructor. To
 * make things worse, the smash::DecayAction::perform function requires a
 * pointer to smash::Particles (which contains a map of all existing
 * smash::ParticleData objects). The \c perform function further calls
 * smash::Action::choose_channel which requires a std::vector of
 * smash::ProcessBranch to determine the the final state particles.
 *
 * \subsection unittest_smash_compromise Compromise
 * We see that testing \c DecayAction in isolation will be hard. If we'd want to
 * follow the purist rules for unit testing we'd have to mock all those classes.
 * Up to now we have not used mocking, as it would create even more work when
 * the design of SMASH changes. We should consider mocking on a case by case
 * basis, but feel free to just use the real thing for now.
 *
 * Instead of creating complete mock classes we can use the BUILD_TESTS macro in
 * actual SMASH classes to easily construct mock objects
 *
 * See \ref doxypage_unit_testing_mocking.
 *
 * \subsection unittest_smash_good_example Good Example
 * While implementing the initial conditions (see
 * \c src/test/initial_conditions.cc at the very end), the test to verify
 * momentum conservation was written first and then the code has been adjusted,
 * until the test passed successfully. This can serve as a positive example of
 * test-driven development.
 *
 ******************************************************************************
 * \section unittest_run Running tests & Test-driven development
 * Tests are built with cmake. They can be disabled via the BUILD_TESTING option
 * of cmake; per default tests are enabled.
 *
 * Once a test is built you will find an executable in the top-level build
 * directory. Every test has a name that is used for the \c .cc file, the
 * executable name, and the \c make target names. Take for example the \c
 * pdgcode test:
 * - The source for the test is in \c src/tests/pdgcode.cc.
 * - The executable will be in \c ${CMAKE_BINARY_DIR}/pdgcode.
 * - The make file will have two targets: \c pdgcode and \c run_pdgcode. The
 *   former will only build the executable, the latter will build and run it.
 *
 * For test-driven development, the run target can be quite handy, since it
 * requires a single command to compile and run a unit test. A vim user, for
 * example, will be able to call `:make run_pdgcode` (possibly mapped to a key,
 * such as F10) and get compiler output and test output into the Quickfix
 * buffer.
 *
 * Finally, if you want to run all tests in the test suite you can first build
 * all tests with
 *
 *     make -j4 all
 *
 * and then run all tests with
 *
 *     ctest -j4
 *
 * . Using \c ctest instead of `make test` allows you to run tests in parallel
 * with the \c -j flag. (Use a number that corresponds to the number of cores on
 * the machine where you're working on, instead of \c -j4.)
 * If you run \c ctest with the \c -V flag you will see more output from the
 * tests. The test output is always captured in a file, so that you don't
 * necessarily have to rerun a failed test to see what happened. You can find it
 * in the \c Testing directory in the build dir.
 *
 ******************************************************************************
 *
 * \section unittest_getstarted Get started
 * In SMASH we use a unit testing framework that was originally developed for
 * the [Vc library](http://code.compeng.uni-frankfurt.de/projects/vc). It
 * simplifies test creation to the bare minimum. The following code suffices to
 * run a test:
 * \code
 * #include "unittest.h"
 *
 * TEST(test_name) {
 *   int test = 1 + 1;
 *   COMPARE(test, 2) << "more details";
 *   VERIFY(1 > 0);
 * }
 * \endcode
 * This creates one test function (called "test_name"). This function is called
 * without any further code and executes two checks. If, for some reason, the
 * compiler would determine that test needs to have the value 3, then the output
 * would be:
   \verbatim
    FAIL: ┍ at /home/mkretz/src/smash/src/tests/testfile.cc:5 (0x40451f):
    FAIL: │ test (3) == 2 (2) -> false more details
    FAIL: ┕ test_name

    Testing done. 0 tests passed. 1 tests failed.
   \endverbatim
 * Let's take a look at what this tells us.
 * 1. The test macro that failed was in testfile.cc in line 5.
 * 2. If you want to look at the disassembly, the failure was at 0x40451f.
 * 3. The \ref COMPARE macro compared the expression `test` against the
 expression
 *    `2`. It shows that `test` had a value of `3` while `2` had a value of `2`
 *    (what a surprise). Since the values are not equal `test == 2` returns \c
 *    false.
 * 4. The \ref COMPARE, \ref FUZZY_COMPARE, \ref VERIFY, and \ref FAIL macros
 can be used as
 *    streams. The output will only appear on failure and will be printed right
 *    after the normal output of the macro.
 * 5. Finally the name of the failed test (the name specified inside the \ref
 TEST()
 *    macro) is printed.
 * 6. At the end of the run, a summary of the test results is shown. This may be
 *    important when there are many \ref TEST functions.
 *
 * If the test passed you'll see:
   \verbatim
    PASS: test_name

    Testing done. 1 tests passed. 0 tests failed.
   \endverbatim
 *
 * You can compile tests with the \c smash_add_unittest macro. You only need to
 * pass it the name of the \c .cc file (without the file extension). So, if your
 * test code above was saved in tests/testfile.cc, then you'd add the line
 * \code
 * smash_add_unittest(testfile)
 * \endcode
 * to the \c CMakeLists.txt .
 * You will then get two new targets that you can build with make: \c testfile
 * and \c run_testtest . The latter can be used to build and run a test quickly
 * in "code - compile - test" cycles in test-driven development.
 */

/**
 * \page doxypage_unit_testing_mocking
 *
 * This is a list of functions and classes that can be useful in unit tests for
 * creating objects that are necessary for testing a class in (more or less)
 * isolation.
 */

/**
 * \addtogroup unittest
 * @{
 */

/**
 * \brief Defines a test function.
 *
 * Consider this to expand to `void
 * function_name()`. The function_name will also be the name that appears in the
 * output after PASS/FAIL.
 */
#define TEST(function_name)

/**
 * \brief Same as above, but expects the code to throw an exception of type \p
 * ExceptionType.
 *
 * If the code does not throw (or throws a different exception),
 * the test is considered failed.
 */
#define TEST_CATCH(function_name, ExceptionType)

/**
 * \brief Define a test function template, with type parameter T, which is
 * specialized for all types in the \p typelist.
 */
#define TEST_TYPES(T, test_name, typelist)

/**
 * \brief Verifies that \p condition is \c true.
 */
#define VERIFY(condition)

/**
 * \brief Verifies that \p test_value is equal to \p reference.
 */
#define COMPARE(test_value, reference)

/**
 * \brief Verifies that the difference between \p test_value and \p reference is
 * smaller than \p allowed_difference.
 *
 * If the test fails the output will show the actual difference between \p
 * test_value and \p reference. If this value is positive \p test_value is too
 * large. If it is negative \p test_value is too small.
 */
#define COMPARE_ABSOLUTE_ERROR(test_value, reference, allowed_difference)

/**
 * \brief Verifies that the difference between \p test_value and \p reference is
 * smaller than `allowed_relative_difference * reference`.
 *
 * If the test fails the output will show the actual difference between \p
 * test_value and \p reference. If this value is positive \p test_value is too
 * large. If it is negative \p test_value is too small.
 *
 * The following example tests that `a` is no more than 1% different from `b`:
 * \code
 * COMPARE_ABSOLUTE_ERROR(a, b, 0.01);
 * \endcode
 *
 * \note This test macro still works even if \p reference is set to 0. It will
 * then compare the difference against `allowed_relative_difference * <smallest
 * positive normalized value of reference type>`.
 */
#define COMPARE_RELATIVE_ERROR(test_value, reference, \
                               allowed_relative_difference)

/**
 * \brief Verifies that \p test_value is equal to \p reference within a
 * pre-defined distance in units of least precision (ulp).
 *
 * If the test fails it will print the distance in ulp between \p test_value and
 * \p reference as well as the maximum allowed distance. Often this difference
 * is not visible in the value because the conversion of a double/float to a
 * string needs to round the value to a sensible length.
 *
 * The allowed distance can be modified by calling:
 * \code
 * vir::test::setFuzzyness<float>(4);
 * vir::test::setFuzzyness<double>(7);
 * \endcode
 *
 * <h3> ulp </h3>
 * Unit of least precision is a unit that is derived from the the least
 * significant bit in the mantissa of a floating-point value. Consider a
 * single-precision number (23 mantissa bits) with exponent \f$e\f$. Then 1
 * ulp is \f$2^{e-23}\f$. Thus, \f$\log_2(u)\f$ signifies the the number
 * incorrect mantissa bits (with \f$u\f$ the distance in ulp).
 *
 * If \p test_value and \p reference have a different exponent the meaning of
 * ulp depends on the variable you look at. The FUZZY_COMPARE code always uses
 * \p reference to determine the magnitude of 1 ulp.
 *
 * Example:
 * The value `1.` is `0x3f800000` in binary. The value
 * `1.00000011920928955078125` with binary representation `0x3f800001`
 * therefore has a distance of 1 ulp.
 * A positive distance means the \p test_value is larger than the \p reference.
 * A negative distance means the \p test_value is smaller than the \p reference.
 * * `FUZZY_COMPARE(1.00000011920928955078125, 1.)` will show a distance of 1
 * * `FUZZY_COMPARE(1., 1.00000011920928955078125)` will show a distance of -1
 *
 * The value `0.999999940395355224609375` with binary representation
 * `0x3f7fffff` has a smaller exponent than `1.`:
 * * `FUZZY_COMPARE(0.999999940395355224609375, 1.)` will show a distance of
 * -0.5
 * * `FUZZY_COMPARE(1., 0.999999940395355224609375)` will show a distance of 1
 *
 * <h3> Comparing to 0 </h3>
 * Distance to 0 is implemented as comparing to
 * <tt>std::numeric_limits<T>::min()</tt>
 * instead and adding 1 to the resulting distance.
 */
#define FUZZY_COMPARE(test_value, reference)

/**
 * \brief Call this to fail a test.
 */
#define FAIL()

namespace vir {
namespace test {

/**
 * \brief Pass code that should fail an assertion to this function.
 */
template <class F>
inline void expect_assert_failure(F &&f);

/**
 * \brief Use this to mark that the failure of a following test is
 * expected.
 *
 * \code
 * TEST(something) {
 *   // this needs to pass
 *   VERIFY(true);
 *   vir::test::expect_failure();
 *   // if this fails, the test will be marked "XFAIL" and the failure
 *   // will not be counted
 *   VERIFY(false);
 *   // this will not be checked anymore, because the TEST stops at the
 *   // (expectedly) failing VERIFY.
 *   VERIFY(true);
 * }
 * \endcode
 */
inline void expect_failure();

}  // namespace test
}  // namespace vir

/**
 * @}
 */
