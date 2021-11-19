/*
	sock.h
		socket read/write
		by Thomas Hahn
		last modified 27 May 14 th
*/

#include <sys/socket.h>

#ifdef DEBUG
#define TERM_RED "\e[31m"
#define TERM_BLUE "\e[34m"
#define TERM_RESET "\e[0m\n"
#define MASTER(s, ...) \
fprintf(stderr, TERM_RED ROUTINE " master %d(%d): " s TERM_RESET, core, getpid(), ##__VA_ARGS__)
#define WORKER(s, ...) \
fprintf(stderr, TERM_BLUE ROUTINE " worker %d(%d): " s TERM_RESET, core, getpid(), ##__VA_ARGS__)
#define DEB_ONLY(...) __VA_ARGS__
#else
#define MASTER(s, ...)
#define WORKER(s, ...)
#define DEB_ONLY(...)
#endif

#ifdef LOW_LEVEL_DEBUG
#define TERM_GREEN "\e[32m"
#define TERM_MAGENTA "\e[35m"
#define READ(s, ...) \
fprintf(stderr, TERM_GREEN ROUTINE " pid %d: read " s TERM_RESET, getpid(), ##__VA_ARGS__)
#define WRITE(s, ...) \
fprintf(stderr, TERM_MAGENTA ROUTINE " pid %d: write " s TERM_RESET, getpid(), ##__VA_ARGS__)
#else
#define READ(s, ...)
#define WRITE(s, ...)
#endif

/*********************************************************************/

#ifndef MSG_WAITALL
/* Windows */
#define MSG_WAITALL 0
#endif

static inline int readsock(cint fd, void *data, csize_t n)
{
  ssize_t got;
  size_t remain = n;
  do got = recv(fd, data, remain, MSG_WAITALL);
  while( got > 0 && (data += got, remain -= got) > 0 );
  READ("%lu bytes at %p from fd %d", n, data, fd);
  return got;
}

/*********************************************************************/

static inline int writesock(cint fd, const void *data, csize_t n)
{
  ssize_t got;
  size_t remain = n;
  do got = send(fd, data, remain, MSG_WAITALL);
  while( got > 0 && (data += got, remain -= got) > 0 );
  WRITE("%lu bytes at %p to fd %d", n, data, fd);
  return got;
}

