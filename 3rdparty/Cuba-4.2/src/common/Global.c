/*
	Global.c
		set global vars
		by Thomas Hahn
		last modified 21 Jul 14 th
*/


#include "stddecl.h"


coreinit cubafun_;
extern int cubaverb_;

#ifdef HAVE_FORK
extern corespec cubaworkers_;
#endif


Extern void SUFFIX(cubaverbose)(cint verb)
{
  cubaverb_ = verb;
}

/*********************************************************************/

Extern void SUFFIX(cubacores)(cint n, cint p)
{
#ifdef HAVE_FORK
  cubaworkers_.ncores = n;
  cubaworkers_.pcores = p;
#endif
}

Extern void SUFFIX(cubaaccel)(cint n, cint p)
{
#ifdef HAVE_FORK
  cubaworkers_.naccel = n;
  cubaworkers_.paccel = p;
#endif
}

/*********************************************************************/

Extern void SUFFIX(cubainit)(subroutine f, void *arg)
{
  cubafun_.initfun = f;
  cubafun_.initarg = arg;
}

/*********************************************************************/

Extern void SUFFIX(cubaexit)(subroutine f, void *arg)
{
  cubafun_.exitfun = f;
  cubafun_.exitarg = arg;
}

