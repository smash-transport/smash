/*
	Divonne.c
		Multidimensional integration by partitioning
		originally by J.H. Friedman and M.H. Wright
		(CERNLIB subroutine D151)
		this version by Thomas Hahn
		last modified 22 Jul 14 th
*/

#define DIVONNE
#define ROUTINE "Divonne"

#include "decl.h"
#include "CSample.c"

/*********************************************************************/

Extern void EXPORT(Divonne)(ccount ndim, ccount ncomp,
  Integrand integrand, void *userdata, cnumber nvec,
  creal epsrel, creal epsabs,
  cint flags, cint seed,
  cnumber mineval, cnumber maxeval,
  cint key1, cint key2, cint key3, ccount maxpass,
  creal border, creal maxchisq, creal mindeviation,
  cnumber ngiven, ccount ldxgiven, real *xgiven,
  cnumber nextra, PeakFinder peakfinder,
  cchar *statefile, Spin **pspin,
  int *pnregions, number *pneval, int *pfail,
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
  t.key1 = key1;
  t.key2 = key2;
  t.key3 = key3;
  t.maxpass = maxpass;
  t.border.upper = 1 - (t.border.lower = border);
  t.maxchisq = maxchisq;
  t.mindeviation = mindeviation;
  t.ngiven = ngiven;
  t.xgiven = xgiven;
  t.ldxgiven = ldxgiven;
  t.nextra = nextra;
  t.peakfinder = peakfinder;
  t.statefile = statefile;
  FORK_ONLY(t.spin = Invalid(pspin) ? NULL : *pspin;)

  *pfail = Integrate(&t, integral, error, prob);
  *pnregions = t.nregions;
  *pneval = t.neval;

  WaitCores(&t, pspin);
}

/*********************************************************************/

Extern void EXPORT(divonne)(ccount *pndim, ccount *pncomp,
  Integrand integrand, void *userdata, cnumber *pnvec,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, cint *pseed,
  cnumber *pmineval, cnumber *pmaxeval,
  cint *pkey1, cint *pkey2, cint *pkey3, ccount *pmaxpass,
  creal *pborder, creal *pmaxchisq, creal *pmindeviation,
  cnumber *pngiven, ccount *pldxgiven, real *xgiven,
  cnumber *pnextra, PeakFinder peakfinder,
  cchar *statefile, Spin **pspin,
  int *pnregions, number *pneval, int *pfail,
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
  t.key1 = *pkey1;
  t.key2 = *pkey2;
  t.key3 = *pkey3;
  t.maxpass = *pmaxpass;
  t.border.upper = 1 - (t.border.lower = *pborder);
  t.maxchisq = *pmaxchisq;
  t.mindeviation = *pmindeviation;
  t.ngiven = *pngiven;
  t.xgiven = xgiven;
  t.ldxgiven = *pldxgiven;
  t.nextra = *pnextra;
  t.peakfinder = peakfinder;
  CString(t.statefile, statefile, statefilelen);
  FORK_ONLY(t.spin = Invalid(pspin) ? NULL : *pspin;)

  *pfail = Integrate(&t, integral, error, prob);
  *pnregions = t.nregions;
  *pneval = t.neval;

  WaitCores(&t, pspin);
}

