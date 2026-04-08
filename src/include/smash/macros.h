/*
 *    Copyright (c) 2013-2015,2017-2018,2020,2022-2023
 *      SMASH Team
 */
#ifndef SRC_INCLUDE_SMASH_MACROS_H_
#define SRC_INCLUDE_SMASH_MACROS_H_

/* support for gcc branch prediction */
#ifdef __GNUC__
#define likely(x) __builtin_expect((x), 1)
#define unlikely(x) __builtin_expect((x), 0)
#else
/// Tell the branch predictor that this expression is likely true.
#define likely(x) (x)
/// Tell the branch predictor that this expression is likely false.
#define unlikely(x) (x)
#endif

/// Mark as deprecated, generating compiler warnings when used.
#define SMASH_DEPRECATED(msg) __attribute__((deprecated(msg)))

/* Handy macros to suppress GNU diagnostics */
#if defined(__GNUC__) && !defined(__clang__)

/// Helper macros to turn macro arguments into _Pragma strings
#define DO_PRAGMA_(x) _Pragma(#x)

/// Extra level of indirection to allow expansion of the argument
#define DO_PRAGMA(x) DO_PRAGMA_(x)

/**
 * Safe suppression of GCC warnings for specific headers. This is needed for
 * some headers that trigger warnings that we cannot fix (e.g. from third-party
 * libraries, like Eigen3) and that we do not want to ignore globally. Example
 * usage:
 *
 * \code
 * GCC_IGNORE_BEGIN("-Wnull-dereference")
 * GCC_IGNORE_ALSO("-Wcast-align")
 *   #include <Eigen/Dense>
 * GCC_IGNORE_END()
 * \endcode
 *
 * Note that at the moment this is restricted to the GNU compiler and other
 * compilers are ignored. If future needs arise, this can be extended to support
 * other compilers as well.
 */
#define GCC_IGNORE_BEGIN(warning) \
  DO_PRAGMA(GCC diagnostic push)  \
  DO_PRAGMA(GCC diagnostic ignored warning)

/// Add more warnings while still inside the block
#define GCC_IGNORE_ALSO(warning) DO_PRAGMA(GCC diagnostic ignored warning)

/// End ignoring warnings
#define GCC_IGNORE_END() DO_PRAGMA(GCC diagnostic pop)

#else

/// Non-GCC compilers: fallback empty macro
#define GCC_IGNORE_BEGIN(warning)

/// Non-GCC compilers: fallback empty macro
#define GCC_IGNORE_ALSO(warning)

/// Non-GCC compilers: fallback empty macro
#define GCC_IGNORE_END()

#endif

#endif  // SRC_INCLUDE_SMASH_MACROS_H_
