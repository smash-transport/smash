# Vir's Unit Test Framework

[![license](https://img.shields.io/github/license/mattkretz/virtest.svg)](https://github.com/mattkretz/virtest/blob/master/LICENSE)
[![language](https://img.shields.io/badge/language-C%2B%2B11-blue.svg)](https://isocpp.org/)

[![Build Status](https://travis-ci.org/mattkretz/virtest.svg)](https://travis-ci.org/mattkretz/virtest)
[![Build status](https://ci.appveyor.com/api/projects/status/lxqk5tqs4og6dr3e?svg=true)](https://ci.appveyor.com/project/mattkretz/virtest)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/77236083639c44a19c2d73609349f5a3)](https://www.codacy.com/manual/mattkretz/virtest?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=mattkretz/virtest&amp;utm_campaign=Badge_Grade)

## Why another test framework?

The test framework was developed inside the [Vc](https://github.com/VcDevel/Vc) repository.
The goal was to build a test framework supporting:

* Minimal / no test registration or setup. Just write the test function and you're good.

* Simple way to disable compilation of tests, without having to comment out sections of the source
  file.

* Simple instantiation of test code with types of a given list.

* Support for fuzzy compares of floating point results with fine control over the ULP specification.

* Assertion testing (i.e. verify that assertions fail on violated preconditions).

* Simple but effective output (no XML, JSON, whatever; outputs a recognizable source location for more
  effective test driven development)

All the test frameworks I looked at in 2009 (and 2010) did not even come close to supporting the
above requirements.

Since I'm now very familiar with this test framework I want to use it in my other projects. The only
sensible choice is to release the test framework on its own.

## Usage

### Creating an executable
To write a test executable all you need is to include the test header:
```cpp
#include <vir/test.h>
```

This defines a main function, but at this point there are no tests, so it'll pass with the following
output:
```
 Testing done. 0 tests passed. 0 tests failed. 0 tests skipped.
```

### Creating a test function
Simple test functions are created with the `TEST` macro. Checks inside the test are done with
macros. The need for macros is due to the requirement to output the source location on failure. (The
macros `__FILE__` and `__LINE__` only yield the right value when expanded at the location of the
test.) You can use the following macros:

* `COMPARE(value, reference)`
  Compares `value` against `reference`, requiring the two to be equal. The comparison is done via
  equality operator, if it is usable. If no equality operator is defined for the type a fallback to
  memcmp is done. Note that this may yield incorrect failures if the types contain uninitialized
  padding bits. Also, if the equality operator does not return a boolean, the implementation will
  try to reduce the result to a boolean via calling `all_of(value == reference)`.

* `COMPARE_TYPES(T1, T2)`
  Test whether `T1` and `T2` are the same type (including value category). The 
  comparison is done via `std::is_same`. On failure this test macro prints the 
  `typeid` name wrapped inside a `vir::test::type<T>` template, to not lose 
  information about cv- and ref-qualifiers (`typeid` drops them).

* `FUZZY_COMPARE(value, reference)`
  Verifies that `value` is equal to `reference` within a pre-defined distance 
  in units of least precision (ulp). If the test fails it will print the 
  distance in ulp between `value` and `reference` as well as the maximum 
  allowed distance. Often this difference is not visible in the value because 
  the conversion of a double/float to a string needs to round the value to a 
  sensible length. The allowed distance can be modified by calling:
  ```cpp
  vir::test::setFuzzyness<float>(4);
  vir::test::setFuzzyness<double>(7);
  ```
 
  * ulp: Unit of least precision is a unit that is derived from the least 
    significant bit in the mantissa of a floating-point value. Consider a 
    single-precision number (23 mantissa bits) with exponent *e*. Then 1 ulp is 
    *2ᵉ⁻²³*. Thus, *log₂(u)* signifies the number of incorrect mantissa bits 
    (with *u* the distance in ulp). If `value` and `reference` have a different 
    exponent the meaning of ulp depends on the variable you look at. The 
    `FUZZY_COMPARE` implementation always uses `reference` to determine the 
    magnitude of 1 ulp. Example: The value `1.f` is `0x3f800000` in binary. The 
    value `1.00000011920928955078125f` with binary representation `0x3f800001` 
    therefore has a distance of 1 ulp. A positive distance means the `value` is 
    larger than the `reference`. A negative distance means the `value` is 
    smaller than the `reference`.
    * `FUZZY_COMPARE(1.00000011920928955078125f, 1.f)` will show a distance of 
      1 ulp

    * `FUZZY_COMPARE(1.f, 1.00000011920928955078125f)` will show a distance of 
      -1 ulp
      
    The value `0.999999940395355224609375f` with binary representation 
    `0x3f7fffff` has a smaller exponent than `1.f`:
    * `FUZZY_COMPARE(0.999999940395355224609375f, 1.f)` will show a distance of
      -0.5 ulp

    * `FUZZY_COMPARE(1.f, 0.999999940395355224609375f)` will show a distance of 
      1 ulp
 
  * Comparing to 0: Distance to 0 is implemented as comparing to 
    `std::numeric_limits<T>::min()` instead and adding 1 to the resulting 
    distance.

* `COMPARE_ABSOLUTE_ERROR(value, reference, error)`
  As above, but allowing an absoluted difference between `value` and `reference`.
  Verifies that the difference between `value` and `reference` is smaller than 
  `error`. If the test fails, the output will show the actual difference 
  between `value` and `reference`. If this value is positive `value` is too 
  large. If it is negative `value` is too small.

* `COMPARE_RELATIVE_ERROR(value, reference, error)`
  Verifies that the difference between `value` and `reference` is smaller than 
  `error * reference`. If the test fails, the output will show the actual 
  difference between `value` and `reference`. If this value is positive `value` 
  is too large. If it is negative `value` is too small. The following example 
  tests that `a` is no more than 1% different from `b`:
  ```cpp
  COMPARE_ABSOLUTE_ERROR(a, b, 0.01);
  ```
  If `reference` is set to 0, this macro compares the difference against `error *
  <smallest positive normalized value of reference type>`.

* `MEMCOMPARE(value, reference)`
  Executes a memcmp over the storage bytes of `value` and `reference`. The 
  number of bytes compared is determined via `sizeof`.

* `VERIFY(boolean)`
  Passes if the argument converted to `bool` is `true`. Fails otherwise.

* `FAIL()`
  Immediately fails a test.

* `vir::test::SKIP() << "details"`
  Ends the test, marking it as skipped but not failed.

* `vir::test::ADD_PASS() << "details"`
  Counts and prints an additional passed test

* `vir::test::expect_failure()`
  The next `COMPARE`, `VERIFY`, etc. is expected to fail. The failure will 
  still end the test, but it will print `XFAIL` instead of `FAIL` and will not 
  count as a failed test in the summary.

* `T vir::test::make_value_unknown(const T& x)`
  The value returned from this function will be unknown to the compiler, 
  inhibiting constant propagation optimization passes. This can be important to 
  fully test whether an operation works correctly under all circumstances. Most 
  importantly, some unit tests may compile to nothing (identified as dead code, 
  i.e. code without side effects) if the compiler can infer the result from 
  constant inputs. In such cases it may be important to make test values 
  unknown to the compiler so that runtime behavior is actually tested.

* `NOINLINE(<testable expression>)`
  When a test fails and you want to identify the exact instruction sequence 
  that lead to the failure, then wrapping the expression inside the `COMPARE` 
  or `VERIFY` macro with `NOINLINE` can help you. It places the expression 
  inside a return statement int a lambda which is executed in a function that 
  is guaranteed to not get inlined. Consequently, the instruction pointer 
  printed on failure takes you right after the `vir::test::noinline<...>` call 
  where the failing test was evaluated.

Example:
```cpp
TEST(test_name) {
  VERIFY(1 > 0);
  COMPARE(1, 1);

  struct A { int x; };
  COMPARE(A(), A());     // implicitly does memcmp
  MEMCOMPARE(A(), A());  // explicitly does memcmp
}
```

### Creating a test function instantiated from a typelist
```cpp
TEST_TYPES(T, test_name, (int, float, short)) {
  COMPARE(T(), T(0));
}
```

### Creating a test function that expects an exception
```cpp
TEST_CATCH(test_name, std::exception) {
  throw std::exception();
}
```

This expects the code to throw an exception of the type specified as second 
macro argument. If the code does not throw (or throws a different exception),
the test is considered failed.

### Output additional information on failure
Every compare/verify macro acts as an output stream, usable for printing more 
information in the failure case. Alternatively, the `.on_failure(...)` function 
can be used. Internally, `on_failure` still uses ostream operators to format 
the output string.

Example:
```cpp
TEST(test_name) {
  int test = 3;
  COMPARE(test, 2) << "more " << "details";
  VERIFY(1 > 0).on_failure("or ", "like ", "this");
}
```
Prints:
```
 FAIL: ┍ at tests/testfile.cpp:5 (0x40451f):
 FAIL: │ test (3) == 2 (2) -> false more details
 FAIL: ┕ test_name

 Testing done. 0 tests passed. 1 tests failed.
```

### Default output on failure
Note in the output above it shows:

1. The name of the test function that failed is at the end (`test_name`).

2. The first line points to the source location of the `COMPARE` macro that 
   failed. In this case it's on line 5 of tests/testfile.cpp

3. If you want to inspect the disassembly of the test, the failure was located 
   around 0x40451f.

4. The `COMPARE` macro compared the expression `test`, which had value `3`, 
   against the expression `2`, which had value `2`. The result of `operator==` 
   is `false` (this can be useful information if `operator==` returns a 
   non-bool type).

5. At the end of test executable, a summary of the test results is shown.

### Testing assertions
If you have assertions using `<cassert>`'s `assert(cond)` macro in your code, 
you can `#include <vir/testassert.h>` to replace the standard `assert` macro 
with an implementation of the test framework. This enables two features:

1. If an assertion is triggered from one of the tests, the test framework 
   recognizes it as a test failure and does not call `abort` like the standard 
   definition of the macro does.

2. You can test whether pre-condition violations are recognized by the 
   assertions in your code. Wrap the code that violates the pre-condition with 
   `vir::test::expect_assert_failure([]() { violate_pre_condition(); })`. Now 
   the test fails if the assertion holds.
