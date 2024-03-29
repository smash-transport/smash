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

#endif  // SRC_INCLUDE_SMASH_MACROS_H_
