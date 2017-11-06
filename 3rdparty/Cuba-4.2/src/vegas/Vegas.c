/*
	Vegas.c
		Vegas Monte Carlo integration
		by Thomas Hahn
		last modified 25 Nov 14 th
*/


#define VEGAS
#define ROUTINE "Vegas"

#include "decl.h"
#include "CSample.c"

/*********************************************************************/

Extern void EXPORT(Vegas)(ccount ndim, ccount ncomp,
  Integrand integrand, void *userdata, cnumber nvec,
  creal epsrel, creal epsabs, cint flags, cint seed,
  cnumber mineval, cnumber maxeval,
  cnumber nstart, cnumber nincrease,
  cnumber nbatch, cint gridno,
  cchar *statefile, Spin **pspin,
  number *pneval, int *pfail,
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
  t.seed = seed;
  t.mineval = mineval;
  t.maxeval = maxeval;
  t.nstart = nstart;
  t.nincrease = nincrease;
  t.nbatch = nbatch;
  t.gridno = gridno;
  t.statefile = statefile;
  FORK_ONLY(t.spin = Invalid(pspin) ? NULL : *pspin;)

  *pfail = Integrate(&t, integral, error, prob);
  *pneval = t.neval;

  WaitCores(&t, pspin);
}

/*********************************************************************/

Extern void EXPORT(vegas)(ccount *pndim, ccount *pncomp,
  Integrand integrand, void *userdata, cnumber *pnvec,
  creal *pepsrel, creal *pepsabs, cint *pflags, cint *pseed,
  cnumber *pmineval, cnumber *pmaxeval,
  cnumber *pnstart, cnumber *pnincrease, 
  cnumber *pnbatch, cint *pgridno,
  cchar *statefile, Spin **pspin,
  number *pneval, int *pfail,
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
  t.seed = *pseed;
  t.mineval = *pmineval;
  t.maxeval = *pmaxeval;
  t.nstart = *pnstart;
  t.nincrease = *pnincrease;
  t.nbatch = *pnbatch;
  t.gridno = *pgridno;
  CString(t.statefile, statefile, statefilelen);
  FORK_ONLY(t.spin = Invalid(pspin) ? NULL : *pspin;)

  *pfail = Integrate(&t, integral, error, prob);
  *pneval = t.neval;

  WaitCores(&t, pspin);
}

