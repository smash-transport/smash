// SigmaDM.cc is a part of the PYTHIA event generator.
// Copyright (C) 2017 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// Dark Matter simulation classes.

#include "Pythia8/SigmaDM.h"

namespace Pythia8 {

//==========================================================================

// Sigma2ffbar2Zp2XX.
// Cross section for f fbar' -> Zprime -> XX. (Zprime a.k.a. DMmed(s=1).)

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2ffbar2Zp2XX::initProc() {

  // Store mass and width for propagator.
  mRes      = particleDataPtr->m0(55);
  GammaRes  = particleDataPtr->mWidth(55);
  m2Res     = mRes*mRes;

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(55);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma2ffbar2Zp2XX::sigmaKin() {

  sigma0 = sH / ( pow2(sH - m2Res) + pow2(mRes * GammaRes) );

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma2ffbar2Zp2XX::sigmaHat() {

  // Check for allowed flavour combinations
  if (id1 + id2 != 0 || abs(id1) > 6 ) return 0.;

  // double trace  = 8.0 * pow2(s3 - tH) + pow2(s3 - uH) + 2.0 * s3 * sH ;
  // double sigma = sigma0 * trace;

  // // Colour factors.
  // if (abs(id1) < 7) sigma /= 3.;

  double widthIn = particlePtr->resWidthChan(mRes, abs(id1), -abs(id1));
  double sigma = widthIn * sigma0 * particlePtr->resWidthChan(mRes, 52, -52);

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ffbar2Zp2XX::setIdColAcol() {

  setId(id1, id2, 52, -52);
  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma3ffbar2Zp2XXj.
// Cross section for f fbar' -> Zprime -> XX + jet.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma3ffbar2Zp2XXj::initProc() {

  // Store mass and width for propagator.
  mRes      = particleDataPtr->m0(55);
  GammaRes  = particleDataPtr->mWidth(55);
  m2Res     = mRes*mRes;

  preFac = 8.0 * alpS / 9.0  ;

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(55);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma3ffbar2Zp2XXj::sigmaKin() {

  double propZp = sH / ( pow2(sH - m2Res) + pow2(mRes * GammaRes) );
  sigma0        = preFac * propZp;

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma3ffbar2Zp2XXj::sigmaHat() {

  // Check for allowed flavour combinations
  if (id1 + id2 != 0 || abs(id1) > 6 ) return 0.;

    // Answer.
  return sigma0;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma3ffbar2Zp2XXj::setIdColAcol() {

  setId(id1, id2, 52, -52, 21);

  // Colour flow topologies.
  if (id1 > 0) setColAcol( 1, 0, 0, 2, 0, 0, 0, 0, 1, 2);
  else         setColAcol( 0, 2, 1, 0, 0, 0, 0, 0, 1, 2);

}

//==========================================================================

// Sigma2ffbar2S2XX
// Cross section for f fbar' -> S -> XX

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2gg2S2XX::initProc() {

  // Store mass and width for propagator.
  mRes      = particleDataPtr->m0(54);
  GammaRes  = particleDataPtr->mWidth(54);
  m2Res     = mRes*mRes;

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(54);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma2gg2S2XX::sigmaKin() {

  double propS = pow2(sH) / ( pow2(sH - m2Res) + pow2(mRes * GammaRes) );
  sigma0        = 8. * M_PI * propS;

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma2gg2S2XX::sigmaHat() {

  // Check for allowed flavour combinations
  if (id1 != id2 || abs(id1) != 21 ) return 0.;

  double widthIn  = particlePtr->resWidthChan( mH, 21, 21) / 64.;

  // Width out only includes open channels.
  double widthOut = particlePtr->resWidthChan( mH, 52, -52);

  // Done.
  double sigma = widthIn * sigma0 * widthOut;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2S2XX::setIdColAcol() {

  setId(id1, id2, 52, -52);
  setColAcol( 1, 2, 2, 1, 0, 0);

}

//==========================================================================

// Sigma3gg2S2XXj
// Cross section for f fbar' -> S -> XX + jet

//--------------------------------------------------------------------------

// Initialize process.

void Sigma3gg2S2XXj::initProc() {

  // Store mass and width for propagator.
  mRes      = particleDataPtr->m0(54);
  GammaRes  = particleDataPtr->mWidth(54);
  m2Res     = mRes*mRes;

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(54);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma3gg2S2XXj::sigmaKin() {

  double wid = particlePtr->resWidthChan(mRes, 21, 21);

  Vec4 resonance = p3cm + p4cm;
  double sHCalc = resonance.m2Calc();
  double sHCalc2 = sHCalc * sHCalc;

  Vec4 p1cm = sH * x1Save * Vec4(0., 0., 1., 1.);

  Vec4 tChan = resonance - p1cm;
  double tHCalc = tChan.m2Calc();
  double tHCalc2 = tHCalc * tHCalc;
  Vec4 uChan = p5cm - p1cm;
  double uHCalc = uChan.m2Calc();
  double uHCalc2 = uHCalc * uHCalc;

  propS = sHCalc / ( pow2(sHCalc - m2Res) + pow2(mRes * GammaRes) );

  // Expression adapted from g g -> H g
  sigma0  = (M_PI / sH2) * (3. / 16.) * alpS * (wid / mRes)
    * (sH2 * sH2 + tHCalc2 * tHCalc2 + uHCalc2 * uHCalc2 + pow2(sHCalc2))
    / (sH * tHCalc * uHCalc * sHCalc);

  sigma0 *= propS * sHCalc2 * particlePtr->resWidthChan(mRes, 52, -52);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma3gg2S2XXj::sigmaHat() {

  return sigma0;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma3gg2S2XXj::setIdColAcol() {

  setId(id1, id2, 52, -52, 21);

  if( rndmPtr->flat() < 0.5)
    setColAcol( 1, 2, 3, 1, 0, 0, 0, 0, 3, 2);
  else
    setColAcol( 1, 2, 2, 3, 0, 0, 0, 0, 1, 3);

}

//==========================================================================

} // end namespace Pythia8
