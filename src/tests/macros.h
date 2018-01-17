/*
 *
 *    Copyright (c) 2014-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_TESTS_MACROS_H_
#define SRC_TESTS_MACROS_H_

#ifdef __GNUC__
#define ALWAYS_INLINE_L inline
#define ALWAYS_INLINE_R __attribute__((__always_inline__))
#define ALWAYS_INLINE ALWAYS_INLINE_L ALWAYS_INLINE_R
#define IS_UNLIKELY(x__) __builtin_expect((x__), 0)
#elif defined _MSC_VER
#define ALWAYS_INLINE inline __forceinline
#define ALWAYS_INLINE_L ALWAYS_INLINE
#define ALWAYS_INLINE_R
#define IS_UNLIKELY(x__) x
#else
#define ALWAYS_INLINE
#define ALWAYS_INLINE_L
#define ALWAYS_INLINE_R
#define IS_UNLIKELY(x__) x
#endif

#endif  // SRC_TESTS_MACROS_H_
