/*
 *    Copyright (c) 2012-2014
 *              none
 */
#ifndef SRC_INCLUDE_MACROS_H_
#define SRC_INCLUDE_MACROS_H_

/* support for gcc branch prediction */
#ifdef __GNUC__
#define likely(x)       __builtin_expect((x), 1)
#define unlikely(x)     __builtin_expect((x), 0)
#else
#define likely(x)       (x)
#define unlikely(x)     (x)
#endif

#define SMASH_UNUSED(x) (void)(x)

#endif  // SRC_INCLUDE_MACROS_H_
