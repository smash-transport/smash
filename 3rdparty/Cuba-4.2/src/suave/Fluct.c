/*
	Fluct.c
		compute the fluctuation in the left and right half
		this file is part of Suave
		last modified 14 Mar 15 th
*/


#if defined(HAVE_LONG_DOUBLE) && defined(HAVE_POWL) && REALSIZE <= 10

typedef long double realL;
#define REALL_MAX_EXP LDBL_MAX_EXP
#define REALL_MAX LDBL_MAX
#define powL powl
#define ldexpL ldexpl

#else

typedef real realL;
#define REALL_MAX_EXP REAL_MAX_EXP
#define REALL_MAX REAL_MAX
#define powL powx
#define ldexpL ldexpx

#endif

typedef const realL crealL;

typedef struct {
  realL fluct;
  number n;
} Var;

static inline realL MinL(crealL a, crealL b) {
  return (a < b) ? a : b;
}

static inline realL MaxL(crealL a, crealL b) {
  return (a > b) ? a : b;
}

/*********************************************************************/

static void Fluct(cThis *t, Var *var,
  cBounds *b, creal *w, number n, ccount comp, creal avg, creal err)
{
  creal *x = w + n;
  creal *f = x + n*t->ndim + comp;
  count nvar = 2*t->ndim;
  creal norm = 1/(err*Max(fabs(avg), err));
  creal flat = 2/3./t->flatness;
  crealL max = ldexpL(1., (int)((REALL_MAX_EXP - 2)/t->flatness));

  Clear(var, nvar);

  while( n-- ) {
    count dim;
    crealL arg = 1 + fabs(*w++)*Sq(*f - avg)*norm;
    crealL ft = powL(MinL(arg, max), t->flatness);

    f += t->ncomp;

    for( dim = 0; dim < t->ndim; ++dim ) {
      Var *v = &var[2*dim + (*x++ >= .5*(b[dim].lower + b[dim].upper))];
      crealL f = v->fluct + ft;
      v->fluct = MaxL(f, REALL_MAX/2);
      ++v->n;
    }
  }

  while( nvar-- ) {
    var->fluct = powL(var->fluct, flat);
    ++var;
  }
}

