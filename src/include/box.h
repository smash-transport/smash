/*
 *    Copyright (c) 2012 maximilian attems <attems@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_BOX_H_
#define SRC_INCLUDE_BOX_H_

/* Cube edge length */
extern float A;

/* temporal time step */
extern float EPS;

/* number of steps */
extern int STEPS;

/* Debug runs generate more output */
extern bool verbose;

/* Compile time debug info */
#define DEBUG 1
#ifdef DEBUG
# define printd printf
#else
# define printd(...) ((void)0)
#endif

#endif  // SRC_INCLUDE_BOX_H_
