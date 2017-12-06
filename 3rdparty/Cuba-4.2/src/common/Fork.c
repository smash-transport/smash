/*
	Fork.c
		fork the cores for parallel sampling
		(C version only)
		by Thomas Hahn
		last modified 23 Apr 15 th
*/


#define ROUTINE "cubafork"
#include "stddecl.h"

#ifdef HAVE_FORK

#include "sock.h"

#define MINCORES 1

coreinit cubafun_;
extern int cubaverb_;
extern corespec cubaworkers_;

/*********************************************************************/

static inline void Child(cint fd, cint core)
{
  dispatch d;

  while( readsock(fd, &d, sizeof d) == sizeof d ) {
    if( d.thissize ) {
      MemAlloc(d.thisptr, d.thissize);
      WORKER("reading This (%lu)", d.thissize);
      readsock(fd, d.thisptr, d.thissize);
    }
    WORKER("running %p on fd %d", d.thisptr, fd);
    d.worker(d.thisptr, d.thissize, core, fd);
    if( d.thissize ) free(d.thisptr);
  }
}

/*********************************************************************/

Extern void SUFFIX(cubafork)(Spin **pspin)
{
  char out[128];
  int cores, core;
  fdpid *pfp;
  Spin *spin;

  VerboseInit();

  EnvInit(cubaworkers_.paccel, "CUBAACCELMAX", 1000);
  EnvInit(cubaworkers_.pcores, "CUBACORESMAX", 10000);
  EnvInit(cubaworkers_.naccel, "CUBAACCEL", 0);
  EnvInit(cubaworkers_.ncores, "CUBACORES", -sysconf(_SC_NPROCESSORS_ONLN));

#ifdef HAVE_GETLOADAVG
  if( cubaworkers_.ncores < 0 ) {
    static int load = uninitialized;
    if( load == uninitialized ) {
      double loadavg;
      getloadavg(&loadavg, 1);
      load = floor(loadavg);
    }
    cubaworkers_.ncores = IMax(-cubaworkers_.ncores - load, 0);
  }
#else
  cubaworkers_.ncores = abs(cubaworkers_.ncores);
#endif

  cores = cubaworkers_.naccel + cubaworkers_.ncores;
  if( cores < MINCORES ) {
    *pspin = NULL;
    return;
  }

  if( cubaverb_ ) {
    sprintf(out, "using %d cores %d accelerators via "
#ifdef HAVE_SHMGET
      "shared memory",
#else
      "pipes",
#endif
      cubaworkers_.ncores, cubaworkers_.naccel);
    Print(out);
  }

  fflush(NULL);		/* make sure all buffers are flushed,
			   or else buffered content will be written
			   out multiply, at each child's exit(0) */

  MemAlloc(spin, sizeof *spin + cores*sizeof *spin->fp);
  spin->spec = cubaworkers_;
  pfp = spin->fp;
  for( core = -spin->spec.naccel; core < spin->spec.ncores; ++core ) {
    int fd[2];
    pid_t pid;
    assert(
      socketpair(AF_LOCAL, SOCK_STREAM, 0, fd) != -1 &&
      (pid = fork()) != -1 );
    if( pid == 0 ) {
      close(fd[0]);
      free(spin);
      Child(fd[1], core);
      exit(0);
    }
    MASTER("forked pid %d pipe %d(master) -> %d(worker)",
      pid, fd[0], fd[1]);
    close(fd[1]);
    pfp->fd = fd[0];
    pfp->pid = pid;
    ++pfp;
  }

  *pspin = spin;
}

/*********************************************************************/

Extern void SUFFIX(cubawait)(Spin **pspin)
{
  int cores, core, status;
  Spin *spin;

  MasterExit();

  if( Invalid(pspin) || (spin = *pspin) == NULL ) return;

  cores = spin->spec.naccel + spin->spec.ncores;

  for( core = 0; core < cores; ++core ) {
    MASTER("closing fd %d", spin->fp[core].fd);
    close(spin->fp[core].fd);
  }

#ifdef KILL_WORKERS
  for( core = 0; core < cores; ++core ) {
    MASTER("killing pid %d", spin->fp[core].pid);
    kill(spin->fp[core].pid, SIGKILL);
  }
#endif

  for( core = 0; core < cores; ++core ) {
    DEB_ONLY(pid_t pid;)
    MASTER("waiting for child");
    DEB_ONLY(pid =) wait(&status);
    MASTER("pid %d terminated with exit code %d", pid, status);
  }

  free(spin);
  *pspin = NULL;
}

#else

Extern void SUFFIX(cubafork)(Spin **pspin) {}

Extern void SUFFIX(cubawait)(Spin **pspin)
{
  MasterExit();
}

#endif

