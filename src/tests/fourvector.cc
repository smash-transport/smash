/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#include "../include/FourVector.h"

int main() {
  FourVector a(0.12, 0.06, 0.003, -0.15), b(0.06, 0.03, 0.0015, -0.075);
  FourVector c(0.1, 0.6, 0.3, -0.15), d(0.01, 0.06, 0.0015, -0.75);

  /* check equality */
  if (a == b)
    return -1;

  /* check smaller equal */
  if (c <= d)
    return -1;

  /* check assignment */
  FourVector f = a;
  if (!(f == a))
    return -1;

  /* check addition! */
  FourVector g = b + b;
  if (a != g)
    return -1;

  /* check division */
  a /= 2;
  if (a != b)
    return -1;

  return 0;
}
