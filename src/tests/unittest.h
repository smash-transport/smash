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

#ifndef SRC_TESTS_UNITTEST_H_
#define SRC_TESTS_UNITTEST_H_

#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <tuple>
#include <typeinfo>

#if defined(__GNUC__) && !defined(_WIN32) && defined(_GLIBCXX_OSTREAM)
#define HACK_OSTREAM_FOR_TTY 1
#endif

#ifdef HACK_OSTREAM_FOR_TTY
#include <unistd.h>
#include <ext/stdio_sync_filebuf.h>
#endif

#include "tests/assert.h"
#include "tests/macros.h"
#include "tests/ulp.h"

namespace UnitTest {

using std::vector;
using std::tuple;
using std::get;

namespace AnsiColor {
struct Type {
  const char *data;
};
static const Type green = {"\033[1;40;32m"};
static const Type yellow = {"\033[1;40;33m"};
static const Type blue = {"\033[1;40;34m"};
static const Type normal = {"\033[0m"};
}  // namespace AnsiColor

#ifdef HACK_OSTREAM_FOR_TTY
class hacked_ostream : public std::ostream {
 public:
  using std::ostream::_M_streambuf;
};
static __attribute__((__const__)) bool mayUseColor(const std::ostream &os) {
  std::basic_streambuf<char> *hack1 = const_cast<std::basic_streambuf<char> *>(
      os.*(&hacked_ostream::_M_streambuf));
  __gnu_cxx::stdio_sync_filebuf<char> *hack =
      dynamic_cast<__gnu_cxx::stdio_sync_filebuf<char> *>(hack1);
  if (!hack) {
    return false;
  }
  FILE *file = hack->file();
  return 1 == isatty(fileno(file));
}
#else
static bool mayUseColor(const std::ostream &) {
  return false;
}
#endif

static std::ostream &operator<<(std::ostream &s, const AnsiColor::Type &color) {
  if (mayUseColor(s)) {
    return s << color.data;
  } else {
    return s;
  }
}

static inline void printPass() {
  std::cout << AnsiColor::green << " PASS: " << AnsiColor::normal;
}

class UnitTestFailure {};

typedef void (*TestFunction)(void);

class UnitTester {  // {{{1
 public:
  UnitTester()
      : status(true),
        expect_failure(false),
        assert_failure(0),
        expect_assert_failure(false),
        float_fuzzyness(1.f),
        double_fuzzyness(1.),
        only_name(0),
        m_finalized(false),
        failedTests(0),
        passedTests(0) {
  }

  int finalize() {
    m_finalized = true;
    std::cout << "\n Testing done. " << passedTests << " tests passed. "
              << failedTests << " tests failed." << std::endl;
    return failedTests;
  }

  void runTestInt(TestFunction fun, const char *name);

  bool status;
  bool expect_failure;
  int assert_failure;
  bool expect_assert_failure;
  float float_fuzzyness;
  double double_fuzzyness;
  const char *only_name;

 private:
  bool m_finalized;
  int failedTests;

 public:
  int passedTests;
};

static UnitTester global_unit_test_object_;

static inline void EXPECT_FAILURE() {  // {{{1
  global_unit_test_object_.expect_failure = true;
}

static const char *_unittest_fail() {  // {{{1
  if (global_unit_test_object_.expect_failure) {
    return "XFAIL: ";
  }
  static const char *str = 0;
  if (str == 0) {
    if (mayUseColor(std::cout)) {
      static const char *fail = " \033[1;40;31mFAIL:\033[0m ";
      str = fail;
    } else {
      static const char *fail = " FAIL: ";
      str = fail;
    }
  }
  return str;
}
static void initTest(int argc, char **argv) {  // {{{1
  for (int i = 1; i < argc; ++i) {
    if (0 == std::strcmp(argv[i], "--help") ||
        0 == std::strcmp(argv[i], "-h")) {
      std::cout << "Usage: " << argv[0] << " [-h|--help] [--only <testname>]\n";
      exit(0);
    }
    if (0 == std::strcmp(argv[i], "--only") && i + 1 < argc) {
      global_unit_test_object_.only_name = argv[i + 1];
    }
  }
}
// setFuzzyness {{{1
template <typename T>
static inline void setFuzzyness(T);
template <>
inline void setFuzzyness<float>(float fuzz) {
  global_unit_test_object_.float_fuzzyness = fuzz;
}
template <>
inline void setFuzzyness<double>(double fuzz) {
  global_unit_test_object_.double_fuzzyness = fuzz;
}
void UnitTester::runTestInt(TestFunction fun, const char *name) {  // {{{1
  if (global_unit_test_object_.only_name &&
      0 != std::strcmp(name, global_unit_test_object_.only_name)) {
    return;
  }
  global_unit_test_object_.status = true;
  global_unit_test_object_.expect_failure = false;
  try {
    setFuzzyness<float>(1);
    setFuzzyness<double>(1);
    fun();
  }
  catch (UnitTestFailure) {
  }
  if (global_unit_test_object_.expect_failure) {
    if (!global_unit_test_object_.status) {
      std::cout << "XFAIL: " << name << std::endl;
    } else {
      std::cout
          << "unexpected PASS: " << name
          << "\n    This test should have failed but didn't. Check the code!"
          << std::endl;
      ++failedTests;
    }
  } else {
    if (!global_unit_test_object_.status) {
      std::cout << _unittest_fail() << "┕ " << name << std::endl;
      ++failedTests;
    } else {
      UnitTest::printPass();
      std::cout << name;
      std::cout << std::endl;
      ++passedTests;
    }
  }
}

// ulpDiffToReferenceWrapper {{{1
template <typename T>
T ulpDiffToReferenceWrapper(T a, T b) {
  const T diff = ulpDiffToReference(a, b);
  return diff;
}
// unittest_fuzzyCompareHelper {{{1
template <typename T>
static inline bool unittest_fuzzyCompareHelper(const T &a, const T &b) {
  return a == b;
}
template <>
inline bool unittest_fuzzyCompareHelper<float>(const float &a, const float &b) {
  return ulpDiffToReferenceWrapper(a, b) <=
         global_unit_test_object_.float_fuzzyness;
}
template <>
inline bool unittest_fuzzyCompareHelper<double>(const double &a,
                                                const double &b) {
  return ulpDiffToReferenceWrapper(a, b) <=
         global_unit_test_object_.double_fuzzyness;
}

// unitttest_comparePrintHelper {{{1
template <typename T1, typename T2, typename M>
inline void unitttest_comparePrintHelper(const T1 &a, const T2 &b, const M &m,
                                         const char *aa, const char *bb,
                                         const char *file, int line,
                                         double fuzzyness = 0.) {
  std::cout << "       " << aa << " (" << std::setprecision(10) << a
            << std::setprecision(6) << ") == " << bb << " ("
            << std::setprecision(10) << b << std::setprecision(6) << ") -> "
            << m;
  if (fuzzyness > 0.) {
    std::cout << " with fuzzyness " << fuzzyness;
  }
  std::cout << " at " << file << ":" << line << " failed.\n";
}

// unittest_fuzzynessHelper {{{1
template <typename T>
inline double unittest_fuzzynessHelper(const T &) {
  return 0.;
}
template <>
inline double unittest_fuzzynessHelper<float>(const float &) {
  return global_unit_test_object_.float_fuzzyness;
}
template <>
inline double unittest_fuzzynessHelper<double>(const double &) {
  return global_unit_test_object_.double_fuzzyness;
}

class _UnitTest_Compare {  // {{{1
 public:
  enum OptionFuzzy { Fuzzy };

  // Normal Compare ctor {{{2
  template <typename T1, typename T2>
  ALWAYS_INLINE _UnitTest_Compare(const T1 &a, const T2 &b, const char *_a,
                                  const char *_b, const char *_file, int _line)
      : m_ip(getIp()), m_failed(!(a == b)) {
    if (IS_UNLIKELY(m_failed)) {
      printFirst();
      printPosition(_file, _line);
      print(":\n");
      print(_a);
      print(" (");
      print(std::setprecision(10));
      print(a);
      print(") == ");
      print(_b);
      print(" (");
      print(std::setprecision(10));
      print(b);
      print(std::setprecision(6));
      print(") -> ");
      print(a == b);
    }
  }

  // Fuzzy Compare ctor {{{2
  template <typename T>
  ALWAYS_INLINE _UnitTest_Compare(const T &a, const T &b, const char *_a,
                                  const char *_b, const char *_file, int _line,
                                  OptionFuzzy)
      : m_ip(getIp()), m_failed(!unittest_fuzzyCompareHelper(a, b)) {
    if (IS_UNLIKELY(m_failed)) {
      printFirst();
      printPosition(_file, _line);
      print(":\n");
      print(_a);
      print(" (");
      print(std::setprecision(10));
      print(a);
      print(") ≈ ");
      print(_b);
      print(" (");
      print(std::setprecision(10));
      print(b);
      print(std::setprecision(6));
      print(") -> ");
      print(a == b);
      printFuzzyInfo(a, b);
    }
  }

  // VERIFY ctor {{{2
  ALWAYS_INLINE _UnitTest_Compare(bool good, const char *cond,
                                  const char *_file, int _line)
      : m_ip(getIp()), m_failed(!good) {
    if (IS_UNLIKELY(m_failed)) {
      printFirst();
      printPosition(_file, _line);
      print(": ");
      print(cond);
    }
  }

  // FAIL ctor {{{2
  ALWAYS_INLINE _UnitTest_Compare(const char *_file, int _line)
      : m_ip(getIp()), m_failed(true) {
    printFirst();
    printPosition(_file, _line);
    print(":\n");
  }

  // stream operators {{{2
  template <typename T>
  ALWAYS_INLINE const _UnitTest_Compare &operator<<(const T &x) const {
    if (IS_UNLIKELY(m_failed)) {
      print(x);
    }
    return *this;
  }

  ALWAYS_INLINE const _UnitTest_Compare &operator<<(const char *str) const {
    if (IS_UNLIKELY(m_failed)) {
      print(str);
    }
    return *this;
  }

  ALWAYS_INLINE const _UnitTest_Compare &operator<<(const char ch) const {
    if (IS_UNLIKELY(m_failed)) {
      print(ch);
    }
    return *this;
  }

  ALWAYS_INLINE const _UnitTest_Compare &operator<<(bool b) const {
    if (IS_UNLIKELY(m_failed)) {
      print(b);
    }
    return *this;
  }

  ALWAYS_INLINE ~_UnitTest_Compare()  // {{{2
#ifdef _NO_NOEXCEPT
      throw(UnitTestFailure)
#else
      noexcept(false)
#endif
  {
    if (IS_UNLIKELY(m_failed)) {
      printLast();
    }
  }

  // }}}2
 private:
  static ALWAYS_INLINE size_t getIp() {  // {{{2
    size_t _ip;
#ifdef __GNUC__
#ifdef __x86_64__
    asm volatile("lea 0(%%rip),%0" : "=r"(_ip));
#else
    // asm volatile("call 1f\n\t1: pop %0" : "=r"(_ip));
    asm volatile("1: movl $1b,%0" : "=r"(_ip));
#endif
#else
    _ip = 0;
#endif
    return _ip;
  }
  static void printFirst() {  // {{{2
    std::cout << _unittest_fail() << "┍ ";
  }
  // print overloads {{{2
  template <typename T>
  static inline void print(const T &x) {
    std::cout << x;
  }
  static void print(const std::type_info &x) {
    std::cout << x.name();
  }
  static void print(const char *str) {
    const char *pos = 0;
    if (0 != (pos = std::strchr(str, '\n'))) {
      if (pos == str) {
        std::cout << '\n' << _unittest_fail() << "│ " << &str[1];
      } else {
        char *left = strdup(str);
        left[pos - str] = '\0';
        std::cout << left << '\n' << _unittest_fail() << "│ " << &pos[1];
        free(left);
      }
    } else {
      std::cout << str;
    }
  }
  static void print(const char ch) {
    if (ch == '\n') {
      std::cout << '\n' << _unittest_fail() << "│ ";
    } else {
      std::cout << ch;
    }
  }
  static void print(bool b) {
    std::cout << (b ? "true" : "false");
  }
  // printLast {{{2
  static void printLast() {
    std::cout << std::endl;
    global_unit_test_object_.status = false;
    throw UnitTestFailure();
  }
  // printPosition {{{2
  void printPosition(const char *_file, int _line) {
    std::cout << "at " << _file << ':' << _line << " (0x" << std::hex << m_ip
              << std::dec << ')';
  }
  // printFuzzy... {{{2
  template <typename T>
  static inline void printFuzzyInfo(T, T) {
  }
  template <typename T>
  static inline void printFuzzyInfoImpl(T a, T b, double fuzzyness) {
    print("\ndistance: ");
    print(ulpDiffToReferenceSigned(a, b));
    print(", allowed distance: ");
    print(fuzzyness);
  }
  // member variables {{{2
  const size_t m_ip;
  const bool m_failed;
};
// printFuzzyInfo specializations for float and double {{{1
template <>
inline void _UnitTest_Compare::printFuzzyInfo(float a, float b) {
  printFuzzyInfoImpl(a, b, global_unit_test_object_.float_fuzzyness);
}
template <>
inline void _UnitTest_Compare::printFuzzyInfo(double a, double b) {
  printFuzzyInfoImpl(a, b, global_unit_test_object_.double_fuzzyness);
}

// FUZZY_COMPARE {{{1
// Workaround for clang: The "<< ' '" is only added to silence the warnings
// about unused return values.
#define FUZZY_COMPARE(a, b)                                       \
  UnitTest::_UnitTest_Compare(a, b, #a, #b, __FILE__, __LINE__,   \
                              UnitTest::_UnitTest_Compare::Fuzzy) \
      << ' '
// COMPARE {{{1
#define COMPARE(a, b) \
  UnitTest::_UnitTest_Compare(a, b, #a, #b, __FILE__, __LINE__) << ' '
// VERIFY {{{1
#define VERIFY(cond) \
  UnitTest::_UnitTest_Compare(cond, #cond, __FILE__, __LINE__) << ' '
// FAIL {{{1
#define FAIL() UnitTest::_UnitTest_Compare(__FILE__, __LINE__) << ' '

// ADD_PASS() << "text" {{{1
class ADD_PASS {
 public:
  ADD_PASS() {
    ++global_unit_test_object_.passedTests;
    printPass();
  }
  ~ADD_PASS() {
    std::cout << std::endl;
  }
  template <typename T>
  ADD_PASS &operator<<(const T &x) {
    std::cout << x;
    return *this;
  }
};
// unittest_assert (called from assert macro) {{{1
void unittest_assert(bool cond, const char *code, const char *file, int line) {
  if (!cond) {
    if (global_unit_test_object_.expect_assert_failure) {
      ++global_unit_test_object_.assert_failure;
    } else {
      _UnitTest_Compare(file, line) << "assert(" << code << ") failed.";
    }
  }
}
// EXPECT_ASSERT_FAILURE {{{1
#define EXPECT_ASSERT_FAILURE(code)                                          \
  global_unit_test_object_.expect_assert_failure = true;                     \
  global_unit_test_object_.assert_failure = 0;                               \
  code;                                                                      \
  if (global_unit_test_object_.assert_failure == 0) {                        \
    /* failure expected but it didn't fail */                                \
    std::cout << "       " << #code << " at " << __FILE__ << ":" << __LINE__ \
              << " did not fail as was expected.\n";                         \
    global_unit_test_object_.status = false;                                 \
    throw UnitTestFailure();                                                 \
    return;                                                                  \
  }                                                                          \
  global_unit_test_object_.expect_assert_failure = false
}  // namespace UnitTest
namespace UnitTest {  // {{{1
// typeToString {{{2
template <typename T>
inline std::string typeToString();

template <typename T>
inline std::string typeToString_impl(T) {
  return typeid(T).name();
}

template <typename T>
inline std::string typeToString() {
  return typeToString_impl(T());
}
template <>
inline std::string typeToString<void>() {
  return "";
}

template <>
inline std::string typeToString<long double>() {
  return "long double";
}
template <>
inline std::string typeToString<double>() {
  return "double";
}
template <>
inline std::string typeToString<float>() {
  return " float";
}
template <>
inline std::string typeToString<long long>() {
  return " long long";
}
template <>
inline std::string typeToString<unsigned long long>() {
  return "ulong long";
}
template <>
inline std::string typeToString<long>() {
  return "  long";
}
template <>
inline std::string typeToString<unsigned long>() {
  return " ulong";
}
template <>
inline std::string typeToString<int>() {
  return "   int";
}
template <>
inline std::string typeToString<unsigned int>() {
  return "  uint";
}
template <>
inline std::string typeToString<short>() {
  return " short";
}
template <>
inline std::string typeToString<unsigned short>() {
  return "ushort";
}
template <>
inline std::string typeToString<char>() {
  return "  char";
}
template <>
inline std::string typeToString<unsigned char>() {
  return " uchar";
}
template <>
inline std::string typeToString<signed char>() {
  return " schar";
}
// runAll and TestData {{{2
typedef tuple<TestFunction, std::string> TestData;
vector<TestData> g_allTests;

static void runAll() {
  for (const auto &data : g_allTests) {
    global_unit_test_object_.runTestInt(get<0>(data), get<1>(data).c_str());
  }
}
// class Test {{{2
template <typename T>
class Test {
 public:
  Test(TestFunction fun, std::string name) {
    name += '<' + typeToString<T>() + '>';
    g_allTests.emplace_back(fun, name);
  }
};

// class Test2 {{{2
template <template <typename V> class TestFunctor, typename... TestTypes>
class Test2;

template <template <typename V> class TestFunctor>
class Test2<TestFunctor> {
 protected:
  explicit Test2(const std::string &) {
  }
};

template <template <typename V> class TestFunctor, typename TestType0,
          typename... TestTypes>
class Test2<TestFunctor, TestType0, TestTypes...> : public Test2<TestFunctor,
                                                                 TestTypes...> {
  typedef Test2<TestFunctor, TestTypes...> Base;

 public:
  static void call() {
    TestFunctor<TestType0>()();
  }

  explicit Test2(std::string name) : Base(name) {
    name += '<' + typeToString<TestType0>() + '>';
    g_allTests.emplace_back(&call, name);
  }
};
// hackTypelist {{{2
template <template <typename V> class F, typename... Typelist>
UnitTest::Test2<F, Typelist...> hackTypelist(void (*)(Typelist...));
}  // namespace UnitTest
// TEST_BEGIN / TEST_END / TEST macros {{{1
#define TEST_BEGIN(V__, fun__, typelist__)                                     \
  template <typename V__>                                                      \
  struct fun__;                                                                \
  static auto test_##fun__##__ = decltype(                                     \
      UnitTest::hackTypelist<fun__>(std::declval<void typelist__>()))(#fun__); \
  template <typename V__>                                                      \
  struct fun__ {                                                               \
    void operator()() {
#define TEST_END \
  }              \
  };

#define TEST(fun__)                                             \
  void fun__();                                                 \
  static UnitTest::Test<void> test_##fun__##__(&fun__, #fun__); \
  void fun__()
// main {{{1
int main(int argc, char **argv) {
  UnitTest::initTest(argc, argv);
  UnitTest::runAll();
  return UnitTest::global_unit_test_object_.finalize();
}
// }}}1
#endif  // SRC_TESTS_UNITTEST_H_

// vim: foldmethod=marker sw=2
