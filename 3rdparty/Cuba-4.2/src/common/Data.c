/*
	Data.c
		initialized data for Cuba
		by Thomas Hahn
		last modified 21 Jul 14 th
*/


#include "stddecl.h"

int cubaverb_ = uninitialized;

#ifdef HAVE_FORK
corespec cubaworkers_ = {
  uninitialized, uninitialized, 
  uninitialized, uninitialized };
#endif

