/*
	mkdist.c

	Usage: mkdist cvfz packagename.tar.gz packagedir files

	Creates packagename.tar.gz for distribution which unpacks
	into the directory "packagedir".
	1) Sets up a temporary tree in which
	   - symlinks are preserved if they point to files in the tree,
	   - all other files are hardlinked.
	2) Tars that tree.
	3) Removes the temporary tree.

	last modified 4 Dec 14 th
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <limits.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

char path[PATH_MAX], *path_end;


static inline char *copydir(char *d, const char *s) {
  const char *dir = strrchr(s, '/');
  const ssize_t n = (dir) ? dir - s + 1 : 0;
  return memcpy(d, s, n) + n;
}


static int depth(const char *path)
{
  int n = 0;
  for( ; ; ) {
    if( strncmp(path, "./", 2) == 0 ) path += 2;
    else if( strncmp(path, "../", 3) == 0 ) {
      path += 3;
      if( --n < 0 ) return n;
    }
    else if( (path = strchr(path, '/')) ) ++n;
    else return n;
    path = path + strspn(path, "/");
  }
}


static void mkdirhier(char *file)
{
  struct stat st;
  char *s = file;
  const char *e = strrchr(file, '/');

  while( s < e ) {
    s = strchr(s, '/');
    *s = 0;
    if( stat(file, &st) ) mkdir(file, 0777);
    else if( !S_ISDIR(st.st_mode) ) {
      fprintf(stderr, "Cannot create %s\n", file);
      exit(1);
    }
    *s++ = '/';
  }
}


static void inspect(const char *file, const int phase)
{
  struct stat st;
  if( lstat(file, &st) ) {
    fprintf(stderr, "Cannot stat %s\n", file);
    return;
  }

  strcpy(path_end, file);
  if( phase == 0 ) mkdirhier(path);

  if( S_ISREG(st.st_mode) ) {
    if( phase == 0 ) link(file, path);
    return;
  }

  if( S_ISDIR(st.st_mode) ) {
    char sub[PATH_MAX], *sub_end;
    struct dirent *entry;

    DIR *cwd = opendir(file);
    if( cwd == NULL ) {
      fprintf(stderr, "Cannot read directory %s\n", file);
      exit(1);
    }

    strcpy(sub, file);
    sub_end = sub + strlen(sub);
    *sub_end++ = '/';

    while( (entry = readdir(cwd)) )
      if( *entry->d_name != '.' ) {
        strcpy(sub_end, entry->d_name);
        inspect(sub, phase);
      }

    closedir(cwd);
    return;
  }

  if( S_ISLNK(st.st_mode) ) {
    char src[PATH_MAX], tmp[PATH_MAX];
    char *lnrel = copydir(src, file), *lnabs;
    ssize_t n = readlink(file, lnrel, PATH_MAX);
    if( n < 0 ) {
      fprintf(stderr, "Cannot read link %s\n", file);
      exit(1);
    }
    lnrel[n++] = 0;

    if( *(lnabs = lnrel) != '/' &&
        depth(lnabs = src) >= 0 ) {	/* not out-of-tree */
      if( phase == 0 ) return;
      memcpy(copydir(tmp, path), lnrel, n);
      if( stat(tmp, &st) == 0 ) {
        symlink(lnrel, path);
        return;
      }
    } else if( phase == 1 ) return;

    if( realpath(lnabs, tmp) == NULL ||
        link(tmp, path) == -1 )
      fprintf(stderr, "Dangling link %s\n", file);
  }
}


int main(int argc, char **argv)
{
  char cmd[PATH_MAX];
  struct stat st;
  int c;

  if( argc < 5 ) {
    fprintf(stderr,
      "Usage:\t%s tarflags packagename.tar[.gz] packagedir files\n\n"
      "Creates packagename.tar[.gz] for distribution which contains\n"
      "\"files\" and unpacks into the directory \"packagedir\".\n"
      "Symlinks are preserved if they point to files in the package.\n\n",
      argv[0]);
    exit(1);
  }

  if( strstr(argv[2], ".tar") == NULL ) {
    fprintf(stderr, "%s is not a tar file\n", argv[2]);
    exit(1);
  }

  if( stat(argv[3], &st) == 0 ) {
    fprintf(stderr, "%s exists already\n", argv[3]);
    exit(1);
  }
  strcpy(path, argv[3]);
  path_end = path + strlen(path);
  *path_end++ = '/';

  sprintf(cmd, "rm -fr %s %s", path, argv[2]);
  system(cmd);

  for( c = 4; c < argc; ++c ) inspect(argv[c], 0);
  for( c = 4; c < argc; ++c ) inspect(argv[c], 1);

  sprintf(cmd, "tar %s %s --owner=root --group=root %s",
    argv[1], argv[2], argv[3]);
  system(cmd);

  sprintf(cmd, "rm -fr %s", argv[3]);
  system(cmd);

  return 0;
}

