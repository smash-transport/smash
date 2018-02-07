// PythiaStdlib.cc is a part of the PYTHIA event generator.
// Copyright (C) 2017 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the Gamma function.

#include "Pythia8/PythiaStdlib.h"

namespace Pythia8 {

//==========================================================================

// Convert string to lowercase for case-insensitive comparisons.
// By default remove any initial and trailing blanks or escape characters.

string toLower(const string& name, bool trim) {

  // Copy string without initial and trailing blanks or escape characters.
  string temp = name;
  if (trim) {
    if (name.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos) return "";
    int firstChar = name.find_first_not_of(" \n\t\v\b\r\f\a");
    int lastChar  = name.find_last_not_of(" \n\t\v\b\r\f\a");
    temp          = name.substr( firstChar, lastChar + 1 - firstChar);
  }

  // Convert to lowercase letter by letter.
  for (int i = 0; i < int(temp.length()); ++i) temp[i] = tolower(temp[i]);
  return temp;

}

//--------------------------------------------------------------------------

// The Gamma function for real arguments, using the Lanczos approximation.
// Code based on http://en.wikipedia.org/wiki/Lanczos_approximation

double GammaCoef[9] = {
     0.99999999999980993,     676.5203681218851,   -1259.1392167224028,
      771.32342877765313,   -176.61502916214059,    12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};

double GammaReal(double x) {

  // Reflection formula (recursive!) for x < 0.5.
  if (x < 0.5) return M_PI / (sin(M_PI * x) * GammaReal(1 - x));

  // Iterate through terms.
  double z = x - 1.;
  double gamma = GammaCoef[0];
  for (int i = 1; i < 9; ++i) gamma += GammaCoef[i] / (z + i);

  // Answer.
  double t = z + 7.5;
  gamma *= sqrt(2. * M_PI) * pow(t, z + 0.5) * exp(-t);
  return gamma;

}

//--------------------------------------------------------------------------

// Polynomial approximation for modified Bessel function of the first kind.
// Based on Abramowitz & Stegun, Handbook of mathematical functions (1964).

double besselI0(double x){

  // Parametrize in terms of t.
  double result = 0.;
  double t  = x / 3.75;
  double t2 = pow2(t);

  // Only positive values relevant.
  if      ( t < 0.) ;
  else if ( t < 1.) {
    result = 1.0 + 3.5156229 * t2 + 3.0899424 * pow2(t2)
           + 1.2067492 * pow3(t2) + 0.2659732 * pow4(t2)
           + 0.0360768 * pow5(t2) * 0.0045813 * pow6(t2);
  } else {
    double u = 1. / t;
    result = exp(x) / sqrt(x) * ( 0.39894228 + 0.01328592 * u
           + 0.00225319 * pow2(u) - 0.00157565 * pow3(u)
           + 0.00916281 * pow4(u) - 0.02057706 * pow5(u)
           + 0.02635537 * pow6(u) - 0.01647633 * pow7(u)
           + 0.00392377 * pow8(u) );
  }

  return result;
}

//--------------------------------------------------------------------------

// Polynomial approximation for modified Bessel function of the first kind.
// Based on Abramowitz & Stegun, Handbook of mathematical functions (1964).

double besselI1(double x){

  // Parametrize in terms of t.
  double result = 0.;
  double t = x / 3.75;
  double t2 = pow2(t);

  // Only positive values relevant.
  if      ( t < 0.) ;
  else if ( t < 1.) {
    result = x * ( 0.5 + 0.87890594 * t2 + 0.51498869 * pow2(t2)
           + 0.15084934 * pow3(t2) + 0.02658733 * pow4(t2)
           + 0.00301532 * pow5(t2) * 0.00032411 * pow6(t2) );
  } else {
    double u = 1. / t;
    result = exp(x) / sqrt(x) * ( 0.39894228 - 0.03988024 * u
           - 0.00368018 * pow2(u) + 0.00163801 * pow3(u)
           - 0.01031555 * pow4(u) + 0.02282967 * pow5(u)
           - 0.02895312 * pow6(u) + 0.01787654 * pow7(u)
           - 0.00420059 * pow8(u) );
  }

  return result;
}

//--------------------------------------------------------------------------

// Polynomial approximation for modified Bessel function of a second kind.
// Based on Abramowitz & Stegun, Handbook of mathematical functions (1964).

double besselK0(double x){

  double result = 0.;

  // Polynomial approximation valid ony for x > 0.
  if      ( x < 0.) ;
  else if ( x < 2.) {
    double x2 = pow2(0.5 * x);
    result = -log(0.5 * x) * besselI0(x) - 0.57721566
           + 0.42278420 * x2 + 0.23069756 * pow2(x2)
           + 0.03488590 * pow3(x2) + 0.00262698 * pow4(x2)
           + 0.00010750 * pow5(x2) + 0.00000740 * pow6(x2);
  } else {
    double z = 2. / x;
    result = exp(-x) / sqrt(x) * ( 1.25331414 - 0.07832358 * z
           + 0.02189568 * pow2(z) - 0.01062446 * pow3(z)
           + 0.00587872 * pow4(z) - 0.00251540 * pow5(z)
           + 0.00053208 * pow6(z) );
  }

  return result;
}

//--------------------------------------------------------------------------

// Polynomial approximation for modified Bessel function of a second kind.
// Based on Abramowitz & Stegun, Handbook of mathematical functions (1964).

double besselK1(double x){

  double result = 0.;

  // Polynomial approximation valid ony for x > 0.
  if      ( x < 0.) ;
  else if ( x < 2.) {
    double x2 = pow2(0.5 * x);
    result = log(0.5 * x) * besselI1(x) + 1./x * ( 1. + 0.15443144 * x2
           - 0.67278579 * pow2(x2) - 0.18156897 * pow3(x2)
           - 0.01919402 * pow4(x2) - 0.00110404 * pow5(x2)
           - 0.00004686 * pow6(x2) );
  } else {
    double z = 2. / x;
    result = exp(-x) / sqrt(x) * ( 1.25331414 + 0.23498619 * z
           - 0.03655620 * pow2(z) + 0.01504268 * pow3(z)
           - 0.00780353 * pow4(z) + 0.00325614 * pow5(z)
           - 0.00068245 * pow6(z) );
  }

  return result;
}

//==========================================================================

} // end namespace Pythia8
