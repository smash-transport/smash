/*
	cuba.h
		Prototypes for the Cuba library
		this file is part of Cuba
		last modified 13 Mar 15 th
*/

typedef double cubareal;

	/* integrand_t is intentionally a minimalistic integrand type.
	   It includes neither the nvec and core arguments nor the
	   extra arguments passed by Vegas/Suave (weight, iter) and
	   Divonne (phase).
	   In most cases, integrand_t is just what you want, otherwise
	   simply use an explicit typecast to integrand_t in the Cuba
	   invocation. */
typedef int (*integrand_t)(const int *ndim, const cubareal x[],
  const int *ncomp, cubareal f[], void *userdata);

typedef void (*peakfinder_t)(const int *ndim, const cubareal b[],
  int *n, cubareal x[], void *userdata);

#ifdef __cplusplus
extern "C" {
#endif

void Vegas(const int ndim, const int ncomp,
  integrand_t integrand, void *userdata, const int nvec,
  const cubareal epsrel, const cubareal epsabs,
  const int flags, const int seed,
  const int mineval, const int maxeval,
  const int nstart, const int nincrease, const int nbatch,
  const int gridno, const char *statefile, void *spin,
  int *neval, int *fail,
  cubareal integral[], cubareal error[], cubareal prob[]);

void llVegas(const int ndim, const int ncomp,
  integrand_t integrand, void *userdata, const long long int nvec,
  const cubareal epsrel, const cubareal epsabs,
  const int flags, const int seed,
  const long long int mineval, const long long int maxeval,
  const long long int nstart, const long long int nincrease,
  const long long int nbatch,
  const int gridno, const char *statefile, void *spin,
  long long int *neval, int *fail,
  cubareal integral[], cubareal error[], cubareal prob[]);

void Suave(const int ndim, const int ncomp,
  integrand_t integrand, void *userdata, const int nvec,
  const cubareal epsrel, const cubareal epsabs,
  const int flags, const int seed,
  const int mineval, const int maxeval,
  const int nnew, const int nmin,
  const cubareal flatness, const char *statefile, void *spin,
  int *nregions, int *neval, int *fail,
  cubareal integral[], cubareal error[], cubareal prob[]);

void llSuave(const int ndim, const int ncomp,
  integrand_t integrand, void *userdata, const long long int nvec,
  const cubareal epsrel, const cubareal epsabs,
  const int flags, const int seed,
  const long long int mineval, const long long int maxeval,
  const long long int nnew, const long long int nmin,
  const cubareal flatness, const char *statefile, void *spin,
  int *nregions, long long int *neval, int *fail,
  cubareal integral[], cubareal error[], cubareal prob[]);

void Divonne(const int ndim, const int ncomp,
  integrand_t integrand, void *userdata, const int nvec,
  const cubareal epsrel, const cubareal epsabs,
  const int flags, const int seed,
  const int mineval, const int maxeval,
  const int key1, const int key2, const int key3, const int maxpass,
  const cubareal border, const cubareal maxchisq, const cubareal mindeviation,
  const int ngiven, const int ldxgiven, cubareal xgiven[],
  const int nextra, peakfinder_t peakfinder,
  const char *statefile, void *spin,
  int *nregions, int *neval, int *fail,
  cubareal integral[], cubareal error[], cubareal prob[]);

void llDivonne(const int ndim, const int ncomp,
  integrand_t integrand, void *userdata, const long long int nvec,
  const cubareal epsrel, const cubareal epsabs,
  const int flags, const int seed,
  const long long int mineval, const long long int maxeval,
  const int key1, const int key2, const int key3, const int maxpass,
  const cubareal border, const cubareal maxchisq, const cubareal mindeviation,
  const long long int ngiven, const int ldxgiven, cubareal xgiven[],
  const long long int nextra, peakfinder_t peakfinder,
  const char *statefile, void *spin,
  int *nregions, long long int *neval, int *fail,
  cubareal integral[], cubareal error[], cubareal prob[]);

void Cuhre(const int ndim, const int ncomp,
  integrand_t integrand, void *userdata, const int nvec,
  const cubareal epsrel, const cubareal epsabs,
  const int flags, const int mineval, const int maxeval,
  const int key,
  const char *statefile, void *spin,
  int *nregions, int *neval, int *fail,
  cubareal integral[], cubareal error[], cubareal prob[]);

void llCuhre(const int ndim, const int ncomp,
  integrand_t integrand, void *userdata, const long long int nvec,
  const cubareal epsrel, const cubareal epsabs,
  const int flags,
  const long long int mineval, const long long int maxeval,
  const int key,
  const char *statefile, void *spin,
  int *nregions, long long int *neval, int *fail,
  cubareal integral[], cubareal error[], cubareal prob[]);

void cubafork(void *pspin);
void cubawait(void *pspin);

void cubacores(const int n, const int p);
void cubaaccel(const int n, const int p);

void cubainit(void (*f)(), void *arg);
void cubaexit(void (*f)(), void *arg);

#ifdef __cplusplus
}
#endif

