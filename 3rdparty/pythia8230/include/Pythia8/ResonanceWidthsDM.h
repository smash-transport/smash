// ResonanceWidthsDM.h is a part of the PYTHIA event generator.
// Copyright (C) 2017 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for DM resonance properties: dynamical widths etc.
// ResonanceS, ResonanceZp: derived classes for individual resonances.

#ifndef Pythia8_ResonanceWidthsDM_H
#define Pythia8_ResonanceWidthsDM_H

#include "Pythia8/Settings.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/ResonanceWidths.h"

namespace Pythia8 {

//==========================================================================

// The ResonanceS class. (S a.k.a. DMmed(s=0), PDG id 54.)

class ResonanceS : public ResonanceWidths {

public:

  // Constructor.
  ResonanceS(int idResIn) {initBasic(idResIn);}

private:

  // Couplings.
  double gq, gX;
  double preFac, alpS;
  bool pScalar;

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool calledFromInit = false);

  virtual double eta2gg();

};

//==========================================================================

// The ResonanceZp class. (Zp a.k.a. DMmed(s=1), PDG id 55.)

class ResonanceZp : public ResonanceWidths {

public:

  // Constructor.
  ResonanceZp(int idResIn) {initBasic(idResIn);}

private:

  // Couplings.
  double vu, vd, vl, vv, vX, au, ad, al, av, aX;
  double preFac, gZp;

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool calledFromInit = false);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_ResonanceWidthsDM_H
