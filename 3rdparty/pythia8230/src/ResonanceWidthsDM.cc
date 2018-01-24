// ResonanceWidthsDM.cc is a part of the PYTHIA event generator.
// Copyright (C) 2017 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for DM resonance properties
// in the ResonanceS and ResonanceZp classes.

#include "Pythia8/ParticleData.h"
#include "Pythia8/ResonanceWidthsDM.h"
#include "Pythia8/PythiaComplex.h"

namespace Pythia8 {

//==========================================================================

// The ResonanceS class.
// Derived class for S properties. (DMmed(s=0), PDG id 54)

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceS::initConstants() {

  // Locally stored properties and couplings.
  double vq = settingsPtr->parm("Sdm:vf");
  double vX = settingsPtr->parm("Sdm:vX");
  double aq = settingsPtr->parm("Sdm:af");
  double aX = settingsPtr->parm("Sdm:aX");

  gq = abs(aq) > 0 ? aq : vq;
  gX = abs(aX) > 0 ? aX : vX;

  if (abs(aX) > 0) pScalar = true;
  else pScalar = false;

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceS::calcPreFac(bool) {

  // Common coupling factors.
  preFac      = 1.0 / (12.0 * M_PI * mRes);
  alpS = couplingsPtr->alphaS(mHat * mHat);

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceS::calcWidth(bool) {

  // Check that above threshold.
  if (ps == 0.) return;

  double mRat2 = pow2(mf1 / mRes);
  double kinfac = (1 - 4 * mRat2) * (1. + 2 * mRat2);

  widNow = 0.;

  if (id1Abs == 21)
    widNow = pow2(gq) * preFac * pow2(alpS / M_PI) * eta2gg();

  if(id1Abs < 7)
    widNow = 3. * pow2(gq * mf1) * preFac * kinfac;

  if(id1Abs == 52)
    widNow = pow2(gX * mf1) * preFac * kinfac;

}

//--------------------------------------------------------------------------

double ResonanceS::eta2gg() {

  // Initial values.
  complex eta = complex(0., 0.);
  double  mLoop, epsilon, root, rootLog;
  complex phi, etaNow;

  // Loop over s, c, b, t quark flavours.
  for (int idNow = 3; idNow < 7; ++idNow) {
    mLoop   = particleDataPtr->m0(idNow);
    epsilon = pow2(2. * mLoop / mHat);

    // Value of loop integral.
    if (epsilon <= 1.) {
      root    = sqrt(1. - epsilon);
      rootLog = (epsilon < 1e-4) ? log(4. / epsilon - 2.)
                : log( (1. + root) / (1. - root) );
      phi     = complex( -0.25 * (pow2(rootLog) - pow2(M_PI)),
                0.5 * M_PI * rootLog );
    }
    else phi  = complex( pow2( asin(1. / sqrt(epsilon)) ), 0.);

    // Factors that depend on Higgs and flavour type.
    if (!pScalar) etaNow = -0.5 * epsilon
       * (complex(1., 0.) + (1. - epsilon) * phi);
    else etaNow = -0.5 * epsilon * phi;


    // Sum up contribution and return square of absolute value.
    eta += etaNow;
  }
  return (pow2(eta.real()) + pow2(eta.imag()));

}

//==========================================================================

// The ResonanceZp class.
// Derived class for Z'^0 properties. (DMmed(s=1), PDG id 55)

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceZp::initConstants() {

  // Locally stored properties and couplings.
  gZp = settingsPtr->parm("Zp:gZp");
  vu = settingsPtr->parm("Zp:vu");
  vd = settingsPtr->parm("Zp:vd");
  vl = settingsPtr->parm("Zp:vl");
  vv = settingsPtr->parm("Zp:vv");
  vX = settingsPtr->parm("Zp:vX");
  au = settingsPtr->parm("Zp:au");
  ad = settingsPtr->parm("Zp:ad");
  al = settingsPtr->parm("Zp:al");
  av = settingsPtr->parm("Zp:av");
  aX = settingsPtr->parm("Zp:aX");

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceZp::calcPreFac(bool) {

  // Common coupling factors.
  preFac      = pow2(gZp) * mRes / 12.0 / M_PI;

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceZp::calcWidth(bool) {

  // Check that above threshold.
  if (ps == 0. || id1 * id2 > 0) return;

  double kinFacA  = pow3(ps);
  double kinFacV  = ps * (1. + 2. * mr1);

  double fac = 0;
  widNow = 0.;

  if (id1Abs < 7) {
    if (id1Abs %2 ) fac = vu * vu * kinFacV + au * au * kinFacA;
    else fac = vd * vd * kinFacV + ad * ad * kinFacA;
  }
  else if (id1Abs > 7 && id1Abs < 17){
    if (id1Abs %2 ) fac = vl * vl * kinFacV + al * al * kinFacA;
    else fac = vv * vv * kinFacV + av * av * kinFacA;
  }

  if (id1Abs == 52) fac = vX * vX * kinFacV + aX * aX * kinFacA;
  widNow = fac * preFac;
  if (id1Abs < 7) widNow *= 3;

}

//==========================================================================

} // end namespace Pythia8
