/*
	Parallel.c
		the parallel sampling routine
		for the C versions of the Cuba routines
		by Thomas Hahn
		last modified 23 Apr 15 th
*/

#include "sock.h"

#define MINSLICE 10
#define MINCORES 1
/*#define MINCORES 2*/

typedef struct {
  number n, m, i;
  VES_ONLY(count iter;)
  DIV_ONLY(int phase SHM_ONLY(, shmid);)
} Slice;

#if defined HAVE_SHMGET && (defined SUAVE || defined DIVONNE)
#define FRAMECOPY
#endif

Extern void SUFFIX(cubafork)(Spin **);
Extern void SUFFIX(cubawait)(Spin **);

/*********************************************************************/

static inline void DoSampleParallel(This *t, number n, creal *x, real *f
  VES_ONLY(, creal *w, ccount iter))
{
  char out[128];
  Slice slice, rslice;
  fd_set ready;
  int core, abort, running = 0;
  const fdpid *pfp;
  Spin *spin = t->spin;
  cint paccel = spin->spec.paccel;
  cint naccel = IMin(spin->spec.naccel, (n + paccel - 1)/IMax(paccel, 1));
  cnumber nrest = IDim(n - naccel*paccel);
  cint ncores = IMin(spin->spec.ncores, nrest/MINSLICE);
  number pcores = IMin(spin->spec.pcores, nrest/IMax(ncores, 1));
  number nx = nrest - ncores*pcores;
  if( nx >= ncores ) nx = 0;

  t->neval += n;

  if( VERBOSE > 2 ) {
    sprintf(out, "sampling " NUMBER " points each on %d cores",
      pcores, ncores);
    Print(out);
  }

  slice.n = paccel;
  slice.m = IMax(slice.n, pcores);
  slice.i = 0;
  VES_ONLY(slice.iter = iter;)
  DIV_ONLY(slice.phase = t->phase;)

#ifdef DIVONNE
  if( n > t->nframe ) {
    FrameFree(t, Master);
    t->nframe = n;
    FrameAlloc(t, Master);
  }
  SHM_ONLY(slice.shmid = t->shmid;)
#endif

  SHM_ONLY(if( t->shmid != -1 ) {
    slice.m = n;
#ifdef FRAMECOPY
    VES_ONLY(Copy(t->frame, w, n);)
    Copy(t->frame + n*NW, x, n*t->ndim);
#endif
  })

#define PutSamples(fd) do { \
  slice.n = IMin(slice.n, n); \
  MASTER("sending samples (sli:%lu[+" VES_ONLY(NUMBER "w:%lu+") \
    NUMBER "x:%lu]) to fd %d", \
    sizeof slice, VES_ONLY(slice.n, sizeof *w,) \
    slice.n, t->ndim*sizeof *x, fd); \
  writesock(fd, &slice, sizeof slice); \
  SHM_ONLY(if( t->shmid == -1 )) { \
    VES_ONLY(writesock(fd, w, slice.n*sizeof *w); \
             w += slice.n;) \
    writesock(fd, x, slice.n*t->ndim*sizeof *x); \
    x += slice.n*t->ndim; \
  } \
  slice.i += slice.n; \
  n -= slice.n; \
  ++running; \
} while( 0 )

#define GetSamples(fd) do { \
  readsock(fd, &rslice, sizeof rslice); \
  MASTER("reading samples (sli:%lu[+" NUMBER "f:%lu]) from fd %d", \
    sizeof rslice, rslice.n, t->ncomp*sizeof *f, fd); \
  if( rslice.n == -1 ) abort = 1; \
  else SHM_ONLY(if( t->shmid == -1 )) \
    readsock(fd, f + rslice.i*t->ncomp, rslice.n*t->ncomp*sizeof *f); \
  --running; \
} while( 0 )

  ++pcores;
  pfp = spin->fp;
  for( core = -naccel; n && core < ncores; ++core ) {
    cint fd = pfp++->fd;
    pcores -= (core == nx);
    slice.n = (core < 0) ? paccel : pcores;
    PutSamples(fd);
  }

  abort = 0;

  while( running ) {
    int fdmax = 0;

    FD_ZERO(&ready);
    pfp = spin->fp;
    for( core = -naccel; core < ncores; ++core ) {
      cint fd = pfp++->fd;
      FD_SET(fd, &ready);
      fdmax = IMax(fdmax, fd);
    }
    fdmax = select(fdmax + 1, &ready, NULL, NULL, NULL);

    pfp = spin->fp;
    for( core = -naccel; core < ncores; ++core ) {
      cint fd = pfp++->fd;
      if( FD_ISSET(fd, &ready) ) {
        GetSamples(fd);
        if( abort ) break;
        if( n ) PutSamples(fd);
        if( --fdmax == 0 ) break;
      }
    }
  }

  if( abort ) longjmp(t->abort, -99);

#ifdef FRAMECOPY
  if( t->shmid != -1 )
    Copy(f, t->frame + slice.m*(NW + t->ndim), slice.m*t->ncomp);
#endif
}

/*********************************************************************/

static void DoSample(This *t, number n, creal *x, real *f
  VES_ONLY(, creal *w, ccount iter))
{
  if( t->spin == NULL ||
      t->spin->spec.ncores + t->spin->spec.naccel < MINCORES ||
      n < MINCORES*MINSLICE )
    DoSampleSerial(t, n, x, f VES_ONLY(, w, iter));
  else
    DoSampleParallel(t, n, x, f VES_ONLY(, w, iter));
}

/*********************************************************************/

#ifdef DIVONNE

typedef struct {
  number neval, neval_opt, neval_cut;
  count nregions, iregion, retval;
} ExploreResult;

static inline int ExploreParallel(This *t, cint iregion)
{
  Vector(Totals, totals, NCOMP);
  csize_t regionsize = RegionSize;
  Region *region;
  Spin *spin = t->spin;
  cint cores = spin->spec.naccel + spin->spec.ncores;
  int core = t->running;
  int ireg = iregion;

  if( core >= ((iregion < 0) ? 1 : cores) ) {
    fd_set ready;
    int fd = 0, fdmax = 0;
    ExploreResult res;
    count comp, succ;

    FD_ZERO(&ready);
    for( core = 0; core < cores; ++core ) {
      fd = spin->fp[core].fd;
      FD_SET(fd, &ready);
      fdmax = IMax(fd, fdmax);
    }
    select(fdmax + 1, &ready, NULL, NULL, NULL);

    for( core = 0; core < cores; ++core ) {
      fd = spin->fp[core].fd;
      if( FD_ISSET(fd, &ready) ) break;
    }

    --t->running;
    MASTER("reading res + region (res:%lu+reg:%lu) from fd %d",
      sizeof res, regionsize, fd);
    readsock(fd, &res, sizeof res);
    ireg = res.iregion;
    region = RegionPtr(ireg);
    succ = ireg + region->next;
    readsock(fd, region, regionsize);
    if( --res.nregions > 0 ) {
      region->next = t->nregions - ireg;
      EnlargeRegions(t, res.nregions);
      MASTER("reading regions (%dreg:%lu) from fd %d",
        res.nregions, regionsize, fd);
      readsock(fd, RegionPtr(t->nregions), res.nregions*regionsize);
      t->nregions += res.nregions;

      RegionPtr(t->nregions-1)->next = succ - t->nregions + 1;
    }

    MASTER("reading totals (tot:%lu) from fd %d",
      t->ncomp*sizeof(Totals), fd);
    readsock(fd, totals, t->ncomp*sizeof(Totals));
    for( comp = 0; comp < t->ncomp; ++comp )
      t->totals[comp].secondspread =
        Max(t->totals[comp].secondspread, totals[comp].secondspread);

    t->neval += res.neval;
    t->neval_opt += res.neval_opt;
    t->neval_cut += res.neval_cut;

    if( res.retval == -1 ) return -1;
  }

  if( iregion >= 0 ) {
    Slice slice;
    cint fd = spin->fp[core].fd;
    slice.n = 0;
    slice.i = iregion;
    slice.phase = t->phase;
    region = RegionPtr(iregion);
    MASTER("writing region (sli:%lu+sam:%lu+reg:%lu+tot:%lu) to fd %d",
      sizeof slice, sizeof(Samples), regionsize,
      t->ncomp*sizeof(Totals), fd);
    writesock(fd, &slice, sizeof slice);
    writesock(fd, &t->samples[region->isamples], sizeof(Samples));
    writesock(fd, region, regionsize);
    writesock(fd, t->totals, t->ncomp*sizeof(Totals));
    region->depth = 0;
    ++t->running;
  }

  return ireg;
}

/*********************************************************************/

static int Explore(This *t, cint iregion)
{
  if( t->spin == NULL ||
      t->spin->spec.ncores + t->spin->spec.naccel < MINCORES )
    return ExploreSerial(t, iregion);
  else
    return ExploreParallel(t, iregion);
}

#endif

/*********************************************************************/

static void Worker(This *t, const size_t alloc, cint core, cint fd)
{
  Slice slice;

  if( readsock(fd, &slice, sizeof slice) == sizeof slice &&
      slice.n != -1 ) {
#ifdef DIVONNE
    csize_t regionsize = RegionSize;
    Vector(Totals, totals, NCOMP);
    Spin spin = {{0, 0, 0, 0}};		/* no recursive forks */

    t->totals = totals;
    t->spin = &spin;
    t->size = 2*t->ndim + 2;
    AllocRegions(t);
#endif

    if( alloc ) {
#ifndef DIVONNE
      FrameAlloc(t, Worker);
#endif
#if defined DIVONNE || defined CUHRE
      RuleAlloc(t);
#endif
    }
#ifdef SUAVE
    else SHM_ONLY(if( t->shmid == -1 ))
      MemAlloc(t->frame, t->nframe*SAMPLESIZE);
#endif

    if( cubafun_.initfun ) cubafun_.initfun(cubafun_.initarg, &core);

    do {
      number n = slice.n;
      WORKER("received slice.n = " NUMBER, n);
      DIV_ONLY(t->phase = slice.phase;)

      if( n > 0 ) {
        real VES_ONLY(*w,) *x, *f;
        WORKER("reading samples (sli:%lu[+" VES_ONLY(NUMBER "w:%lu+")
          NUMBER "x:%lu]) from fd %d",
          sizeof slice, VES_ONLY(n, sizeof *w,) n, t->ndim*sizeof *x, fd);

#ifdef DIVONNE
        if( slice.m > t->nframe ) {
          FrameFree(t, Worker);
          t->nframe = slice.m;
          SHM_ONLY(t->shmid = slice.shmid;)
          FrameAlloc(t, Worker);
        }
#endif

        VES_ONLY(w = t->frame;)
        x = t->frame + slice.m*NW;
        f = x + slice.m*t->ndim;

        SHM_ONLY(if( t->shmid != -1 ) {
          VES_ONLY(w += slice.i;)
          x += slice.i*t->ndim;
          f += slice.i*t->ncomp;
        }
        else) {
          VES_ONLY(readsock(fd, w, n*sizeof *w);)
          readsock(fd, x, n*t->ndim*sizeof *x);
        }

        slice.n |= SampleRaw(t, n, x, f, core VES_ONLY(, w, slice.iter));
        WORKER("writing samples (sli:%lu[+" NUMBER "f:%lu]) to fd %d",
          sizeof slice, slice.n, t->ncomp*sizeof *f, fd);
        writesock(fd, &slice, sizeof slice);
        if( SHM_ONLY(t->shmid == -1 &&) slice.n != -1 )
          writesock(fd, f, slice.n*t->ncomp*sizeof *f);
      }
#ifdef DIVONNE
      else {
        Samples *samples, psamples;
        ExploreResult res;

        WORKER("reading region (sli:%lu+sam:%lu+reg:%lu+tot:%lu) from fd %d",
          sizeof slice, sizeof psamples, regionsize,
          t->ncomp*sizeof(Totals), fd);
        readsock(fd, &psamples, sizeof psamples);
        readsock(fd, t->region, regionsize);
        readsock(fd, totals, t->ncomp*sizeof(Totals));
        t->nregions = 1;
        t->neval = t->neval_opt = t->neval_cut = 0;

        samples = &t->samples[RegionPtr(0)->isamples];
        if( psamples.n != samples->n ) {
          SamplesFree(samples);
          *samples = psamples;
          SamplesAlloc(t, samples);
        }

        res.retval = ExploreSerial(t, 0);
        res.neval = t->neval;
        res.neval_opt = t->neval_opt;
        res.neval_cut = t->neval_cut;
        res.nregions = t->nregions;
        res.iregion = slice.i;
        WORKER("writing regions (res:%lu+%dreg:%lu+tot:%lu) to fd %d",
          sizeof res, t->nregions, regionsize,
          t->ncomp*sizeof(Totals), fd);
        writesock(fd, &res, sizeof res);
        writesock(fd, t->region, t->nregions*regionsize);
        writesock(fd, totals, t->ncomp*sizeof(Totals));
      }
#endif
    } while( readsock(fd, &slice, sizeof slice) == sizeof slice &&
             slice.n != -1 );

    if( cubafun_.exitfun ) cubafun_.exitfun(cubafun_.exitarg, &core);

#if defined DIVONNE || defined CUHRE
    RuleFree(t);
#endif

    FrameFree(t, Worker);

#ifdef DIVONNE
    free(t->region);
#endif
  }

  WORKER("worker wrapping up");
}

/*********************************************************************/

static inline void ForkCores(This *t)
{
  dispatch d;
  const fdpid *pfp;
  int ncores, core;

  DIV_ONLY(t->running = 0;)

  d.worker = Worker;
  d.thisptr = t;
  d.thissize = sizeof *t;

  if( t->spin == NULL ) {
    SUFFIX(cubafork)(&t->spin);
    if( t->spin == NULL ) return;
    d.thissize = 0;
  }

  pfp = t->spin->fp;
  ncores = t->spin->spec.ncores;
  for( core = -t->spin->spec.naccel; core < ncores; ++core ) {
    cint fd = pfp++->fd;
    writesock(fd, &d, sizeof d);
    if( d.thissize ) writesock(fd, t, d.thissize);
  }
}

/*********************************************************************/

static inline void WaitCores(This *t, Spin **pspin)
{
  if( Invalid(pspin) ) SUFFIX(cubawait)(&t->spin);
  else {
    Slice slice = { .n = -1 };
    cint cores = t->spin->spec.naccel + t->spin->spec.ncores;
    const fdpid *pfp = t->spin->fp;
    int core;
    for( core = 0; core < cores; ++core )
      writesock(pfp[core].fd, &slice, sizeof slice);
    *pspin = t->spin;
    MasterExit();
  }
}

