// SigmaDM.h is a part of the PYTHIA event generator.
// Copyright (C) 2017 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for Dark Matter process differential cross sections.
// Contains classes derived from SigmaProcess.

#ifndef Pythia8_SigmaDM_H
#define Pythia8_SigmaDM_H

#include "Pythia8/PythiaComplex.h"
#include "Pythia8/SigmaProcess.h"

namespace Pythia8 {

//==========================================================================

// A derived class for f fbar' -> Zprime -> X X. (Zprime a.k.a. DMmed(s=1).)

class Sigma2ffbar2Zp2XX : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2Zp2XX() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()       const {return "f fbar -> Zp -> XX";}
  virtual int    code()       const {return 6001;}
  virtual string inFlux()     const {return "qqbar";}
  virtual int    resonanceA() const {return 55;} // Zprime
  virtual bool   convertM2()  const {return true;} // Use |M|^2

private:

  // Parameters set at initialization.
  double mRes, GammaRes, m2Res, sigma0;

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* particlePtr;

};

//==========================================================================

class Sigma3ffbar2Zp2XXj : public Sigma3Process {

public:

  // Constructor.
  Sigma3ffbar2Zp2XXj() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()       const {return "f fbar -> Zp -> XX + jet";}
  virtual int    code()       const {return 6002;}
  virtual string inFlux()     const {return "qqbar";}
  virtual int    resonanceA() const {return 55;} // Zprime

protected:

  // Parameters set at initialization.
  double mRes, GammaRes, m2Res, preFac, sigma0;

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* particlePtr;

};

//==========================================================================

class Sigma3qg2Zp2XXj : public Sigma3ffbar2Zp2XXj {

  // Constructor.
  Sigma3qg2Zp2XXj() {}

  // Info on the subprocess.
  virtual string name()       const {return "q g -> Zp -> XX + jet";}
  virtual int    code()       const {return 6003;}
  virtual string inFlux()     const {return "qg";}
  virtual int    resonanceA() const {return 55;} // Zprime

};

//==========================================================================

class Sigma2gg2S2XX : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2S2XX() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()       const {return "g g -> S -> XX";}
  virtual int    code()       const {return 6011;}
  virtual string inFlux()     const {return "gg";}
  virtual int    resonanceA() const {return 54;} // scalar mediator
  virtual bool   convertM2()  const {return true;} // Use |M|^2

private:

  // Parameters set at initialization.
  double mRes, GammaRes, m2Res, sigma0;

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* particlePtr;

};

//==========================================================================

class Sigma3gg2S2XXj : public Sigma3Process {

public:

  // Constructor.
  Sigma3gg2S2XXj() {}

  // Initialize process.
  virtual void initProc();

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat).
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()       const {return "g g -> S -> XX + jet";}
  virtual int    code()       const {return 6012;}
  virtual string inFlux()     const {return "gg";}
  virtual int    resonanceA() const {return 54;} // scalar mediator

protected:

  // Parameters set at initialization.
  double mRes, GammaRes, m2Res, propS, sigma0;

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* particlePtr;

};

//==========================================================================

class Sigma3qg2S2XXj : public Sigma3gg2S2XXj {

public:

  // Constructor.
  Sigma3qg2S2XXj() {}

  // Info on the subprocess.
  virtual string name()       const {return "q g -> S -> XX + jet";}
  virtual int    code()       const {return 6013;}
  virtual string inFlux()     const {return "qg";}

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia_SigmaDM_H
