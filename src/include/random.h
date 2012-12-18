/* Implements the Polar form of the Box-Muller Transformation
   to generate gaussian random numbers.  Modified to work with
   rand()

   (c) Copyright 1994, Everett F. Carter Jr.
       Permission is granted by the author to use
       this software for any application provided this
       copyright notice is preserved.
*/
#ifndef SRC_INCLUDE_RANDOM_H_
#define SRC_INCLUDE_RANDOM_H_

#include <math.h>
#include <stdlib.h>

extern unsigned int seedp;

/* normal random variate generator
 * mean m, standard deviation s
 */
template <class T> T box_muller(T m, T s) {
  T x1, x2, w1, y1;
  static T y2;
  static int use_last = 0;

  /* use value from previous call */
  if (use_last) {
    y1 = y2;
    use_last = 0;
  } else {
     do {
       x1 = (T) 2.0 * rand_r(&seedp) / RAND_MAX - (T) 1.0;
       x2 = (T) 2.0 * rand_r(&seedp) / RAND_MAX - (T) 1.0;
       w1 = x1 * x1 + x2 * x2;
     } while (w1 >= (T) 1.0);

       w1 = sqrt(((T) - 2.0 * log(w1)) / w1);
       y1 = x1 * w1;
       y2 = x2 * w1;
       use_last = 1;
  }

  return m + y1 * s;
}

// generates gaussian random number with standard deviation s and mean 0
template <class T> T randGauss(T s) {
  return box_muller((T) 0, s);
}
#endif  // SRC_INCLUDE_RANDOM_H_
