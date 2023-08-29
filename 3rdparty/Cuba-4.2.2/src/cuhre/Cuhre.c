/*
	Cuhre.c
		Adaptive integration using cubature rules
		by Thomas Hahn
		last modified 22 Jul 14 th
*/


#define CUHRE
#define ROUTINE "Cuhre"

#include "decl.h"
#include "CSample.c"

/*********************************************************************/

Extern void EXPORT(Cuhre)(ccount ndim, ccount ncomp,
  Integrand integrand, void *userdata, cnumber nvec,
  creal epsrel, creal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  ccount key, cchar *statefile, Spin **pspin,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  This t;

  VerboseInit();

  t.ndim = ndim;
  t.ncomp = ncomp;
  t.integrand = integrand;
  t.userdata = userdata;
  t.nvec = nvec;
  t.epsrel = epsrel;
  t.epsabs = epsabs;
  t.flags = MaxVerbose(flags);
  t.mineval = mineval;
  t.maxeval = maxeval;
  t.key = key;
  t.statefile = statefile;
  FORK_ONLY(t.spin = Invalid(pspin) ? NULL : *pspin;)
 
  *pfail = Integrate(&t, integral, error, prob);
  *pnregions = t.nregions;
  *pneval = t.neval;

  WaitCores(&t, pspin);
}

/*********************************************************************/

Extern void EXPORT(cuhre)(ccount *pndim, ccount *pncomp,
  Integrand integrand, void *userdata, cnumber *pnvec,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, cnumber *pmineval, cnumber *pmaxeval,
  ccount *pkey, cchar *statefile, Spin **pspin,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob, cint statefilelen)
{
  This t;

  VerboseInit();

  t.ndim = *pndim;
  t.ncomp = *pncomp;
  t.integrand = integrand;
  t.userdata = userdata;
  t.nvec = *pnvec;
  t.epsrel = *pepsrel;
  t.epsabs = *pepsabs;
  t.flags = MaxVerbose(*pflags);
  t.mineval = *pmineval;
  t.maxeval = *pmaxeval;
  t.key = *pkey;
  CString(t.statefile, statefile, statefilelen);
  FORK_ONLY(t.spin = Invalid(pspin) ? NULL : *pspin;)

  *pfail = Integrate(&t, integral, error, prob);
  *pnregions = t.nregions;
  *pneval = t.neval;

  WaitCores(&t, pspin);
}

