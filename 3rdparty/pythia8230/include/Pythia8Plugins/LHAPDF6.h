// LHAPDF6.h is a part of the PYTHIA event generator.
// Copyright (C) 2017 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the LHAPDF6 PDF plugin class.

#ifndef Pythia8_LHAPDF6_H
#define Pythia8_LHAPDF6_H

#include "Pythia8/PartonDistributions.h"
#include "LHAPDF/LHAPDF.h"

namespace Pythia8 {

//==========================================================================

// Global tracking of opened PDF sets.

//--------------------------------------------------------------------------

namespace LHAPDF6Interface {

  // Structure to hold all of the data needed to reproduce a set.
  struct pdf_Info {
    ::LHAPDF::PDF *pdf;
    vector< ::LHAPDF::PDF*> pdfs;
    ::LHAPDF::PDFSet pdfSet;
  };
  // Define opened PDF sets global variable.
  map< int, pair< pdf_Info, int> > initializedSets;

}

//==========================================================================

// Provide interface to the LHAPDF6 library of parton densities.

class LHAPDF6 : public PDF {

public:

  // Constructor.
  LHAPDF6(int idBeamIn, string setName, int member,  int,
          Info* infoPtr) : PDF(idBeamIn), id (-1), pdf(0), extrapol(false)
    {init(setName, member, infoPtr);}

  // Destructor.
  ~LHAPDF6();

  // Allow extrapolation beyond boundaries (not implemented).
  void setExtrapolate(bool extrapolIn) {extrapol = extrapolIn;}

private:

  // The LHAPDF objects.
  int id;
  LHAPDF6Interface::pdf_Info pdfInfo;
  ::LHAPDF::PDF *pdf;
  ::LHAPDF::Extrapolator *ext;
  bool extrapol;

  // Initialization of PDF set.
  void init(string setName, int member, Info* infoPtr);

  // Update parton densities.
  void xfUpdate(int id, double x, double Q2);

  // Check whether x and Q2 values fall inside the fit bounds.
  bool insideBounds(double x, double Q2) {
   return ( x > pdf->xMin()  &&  x < pdf->xMax()
        && Q2 > pdf->q2Min() && Q2 < pdf->q2Max() ); }

  // Return the running alpha_s shipped with the LHAPDF set.
  double alphaS(double Q2) { return pdf->alphasQ2(Q2); }

  // Return quark masses used in the PDF fit.
  double muPDFSave, mdPDFSave, mcPDFSave, msPDFSave, mbPDFSave;
  double mQuarkPDF(int id) {
    if (abs(id) == 1) return mdPDFSave;
    if (abs(id) == 2) return muPDFSave;
    if (abs(id) == 3) return msPDFSave;
    if (abs(id) == 4) return mcPDFSave;
    if (abs(id) == 5) return mbPDFSave;
    return -1.;
 }

  // Calculate uncertainties using the LHAPDF prescription.
  void calcPDFEnvelope(int, double, double, int);
  void calcPDFEnvelope(pair<int,int>, pair<double,double>, double, int);
  PDFEnvelope pdfEnvelope;
  PDFEnvelope getPDFEnvelope() { return pdfEnvelope; }
  static const double PDFMINVALUE;

};

//--------------------------------------------------------------------------

// Destructor.

LHAPDF6::~LHAPDF6() {
  map< int, pair<LHAPDF6Interface::pdf_Info, int> >::iterator set =
    LHAPDF6Interface::initializedSets.find(id);
  if (set == LHAPDF6Interface::initializedSets.end()) return;
  set->second.second--;
  if (set->second.second == 0) {
    LHAPDF6Interface::initializedSets.erase(set);
  }

}

//--------------------------------------------------------------------------

const double LHAPDF6::PDFMINVALUE = 1e-10;

//--------------------------------------------------------------------------

// Initialize a parton density function from LHAPDF6.

void LHAPDF6::init(string setName, int member, Info *info) {
  isSet = false;

  // Initialize the set if not initialized.
  id = ::LHAPDF::lookupLHAPDFID(setName, member);
  if (id < 0) {
    info->errorMsg("Error in LHAPDF6::init: unknown PDF " + setName);
    return;
  } else if (LHAPDF6Interface::initializedSets.find(id) ==
             LHAPDF6Interface::initializedSets.end()) {
    pdfInfo.pdfSet = ::LHAPDF::PDFSet(setName);
    pdfInfo.pdfs = pdfInfo.pdfSet.mkPDFs();
    pdfInfo.pdf = ::LHAPDF::mkPDF(id);
    pdf = pdfInfo.pdf;
    if (!pdf) {
      info->errorMsg("Error in LHAPDF6::init: could not initialize PDF "
                     + setName);
      return;
    } else LHAPDF6Interface::initializedSets[id] = make_pair(pdfInfo,0);
  } else {
    pair< LHAPDF6Interface::pdf_Info, int > &set
      = LHAPDF6Interface::initializedSets[id];
    pdfInfo.pdf = set.first.pdf;
    pdfInfo.pdfs = set.first.pdfs;
    pdfInfo.pdfSet = set.first.pdfSet;
    pdf = pdfInfo.pdf;
    set.second++;
  }
  isSet = true;

  // Store quark masses used in PDF fit.
  muPDFSave = pdf->info().get_entry_as<double>("MUp");
  mdPDFSave = pdf->info().get_entry_as<double>("MDown");
  mcPDFSave = pdf->info().get_entry_as<double>("MCharm");
  msPDFSave = pdf->info().get_entry_as<double>("MStrange");
  mbPDFSave = pdf->info().get_entry_as<double>("MBottom");

}

//--------------------------------------------------------------------------

// Give the parton distribution function set from LHAPDF6.

void LHAPDF6::xfUpdate(int, double x, double Q2) {

  // Freeze at boundary value if PDF is evaluated outside the fit region.
  if (x < pdf->xMin() && !extrapol) x = pdf->xMin();
  if (x > pdf->xMax() )    x = pdf->xMax();
  if (Q2 < pdf->q2Min() ) Q2 = pdf->q2Min();
  if (Q2 > pdf->q2Max() ) Q2 = pdf->q2Max();

  // Update values.
  xg     = pdf->xfxQ2(21, x, Q2);
  xu     = pdf->xfxQ2(2,  x, Q2);
  xd     = pdf->xfxQ2(1,  x, Q2);
  xs     = pdf->xfxQ2(3,  x, Q2);
  xubar  = pdf->xfxQ2(-2, x, Q2);
  xdbar  = pdf->xfxQ2(-1, x, Q2);
  xsbar  = pdf->xfxQ2(-3, x, Q2);
  xc     = pdf->xfxQ2(4,  x, Q2);
  xb     = pdf->xfxQ2(5,  x, Q2);
  xgamma = pdf->xfxQ2(22, x, Q2);

  // Subdivision of valence and sea.
  xuVal  = xu - xubar;
  xuSea  = xubar;
  xdVal  = xd - xdbar;
  xdSea  = xdbar;

  // idSav = 9 to indicate that all flavours reset.
  idSav = 9;

}

//--------------------------------------------------------------------------

// Call the PDF uncertainty calculation in the LHAPDF library

void LHAPDF6::calcPDFEnvelope(int idNow, double xNow, double Q2NowIn,
  int valSea) {

  // Freeze at boundary value if PDF is evaluated outside the fit region.
  double x1 = (xNow < pdfInfo.pdf->xMin() && !extrapol)
            ? pdfInfo.pdf->xMin() : xNow;
  if (x1 > pdfInfo.pdf->xMax() ) x1 = pdfInfo.pdf->xMax();
  double Q2Now = (Q2NowIn < pdfInfo.pdf->q2Min() )
               ? pdfInfo.pdf->q2Min() : Q2NowIn;
  if (Q2Now > pdfInfo.pdf->q2Max() ) Q2Now = pdfInfo.pdf->q2Max();

  vector<double> xfCalc(pdfInfo.pdfs.size());
  for(int imem=0; imem < pdfInfo.pdfs.size(); ++imem) {
    if( valSea==0 || (idNow != 1 && idNow != 2 ) ) {
      xfCalc[imem] = pdfInfo.pdfs[imem]->xfxQ2(idNow, x1, Q2Now);
    } else if ( valSea==1 && (idNow == 1 || idNow == 2 ) ) {
      xfCalc[imem] = pdfInfo.pdfs[imem]->xfxQ2(idNow, x1, Q2Now) -
        pdfInfo.pdfs[imem]->xfxQ2(-idNow, x1, Q2Now);
    } else if ( valSea==2 && (idNow == 1 || idNow == 2 ) ) {
      xfCalc[imem] = pdfInfo.pdfs[imem]->xfxQ2(-idNow, x1, Q2Now);
    }
  }
  ::LHAPDF::PDFUncertainty xfErr = pdfInfo.pdfSet.uncertainty(xfCalc);
  pdfEnvelope.centralPDF = xfErr.central;
  pdfEnvelope.errplusPDF = xfErr.errplus;
  pdfEnvelope.errminusPDF = xfErr.errminus;
  pdfEnvelope.errsymmPDF = xfErr.errsymm;
  pdfEnvelope.scalePDF = xfErr.scale;
}

//--------------------------------------------------------------------------

void LHAPDF6::calcPDFEnvelope(pair<int,int> idNows, pair<double,double> xNows,
  double Q2NowIn, int valSea) {

  // Freeze at boundary value if PDF is evaluated outside the fit region.
  double x1 = (xNows.first < pdfInfo.pdf->xMin() && !extrapol)
            ? pdfInfo.pdf->xMin() : xNows.first;
  if (x1 > pdfInfo.pdf->xMax() ) x1 = pdfInfo.pdf->xMax();
  double x2 = (xNows.second < pdfInfo.pdf->xMin() && !extrapol)
            ? pdfInfo.pdf->xMin() : xNows.second;
  if (x2 > pdfInfo.pdf->xMax() ) x2 = pdfInfo.pdf->xMax();
  double Q2Now = (Q2NowIn < pdfInfo.pdf->q2Min() )
               ? pdfInfo.pdf->q2Min() : Q2NowIn;
  if (Q2Now > pdfInfo.pdf->q2Max() ) Q2Now = pdfInfo.pdf->q2Max();

  vector<double> xfCalc(pdfInfo.pdfs.size());
  for(int imem=0; imem < pdfInfo.pdfs.size(); ++imem) {
    if ( valSea==0 || (idNows.first != 1 && idNows.first != 2 ) ) {
      xfCalc[imem] = pdfInfo.pdfs[imem]->xfxQ2(idNows.first, x1, Q2Now);
    } else if ( valSea==1 && (idNows.first == 1 || idNows.first == 2 ) ) {
      xfCalc[imem] = pdfInfo.pdfs[imem]->xfxQ2(idNows.first, x1, Q2Now)
        - pdfInfo.pdfs[imem]->xfxQ2(-idNows.first, x1, Q2Now);
    } else if ( valSea==2 && (idNows.first == 1 || idNows.first == 2 ) ) {
      xfCalc[imem] = pdfInfo.pdfs[imem]->xfxQ2(-idNows.first, x1, Q2Now);
    }
    if( valSea==0 || (idNows.second != 1 && idNows.second != 2 ) ) {
      xfCalc[imem] /= max(PDFMINVALUE,pdfInfo.pdfs[imem]->xfxQ2(idNows.second,
        x2, Q2Now));
    } else if ( valSea==1 && (idNows.second == 1 || idNows.second == 2 ) ) {
      xfCalc[imem] /= max(pdfInfo.pdfs[imem]->xfxQ2(idNows.second, x2, Q2Now)
        - pdfInfo.pdfs[imem]->xfxQ2(-idNows.second, x2, Q2Now),PDFMINVALUE);
    } else if ( valSea==2 && (idNows.second == 1 || idNows.second == 2 ) ) {
      xfCalc[imem] /= max( pdfInfo.pdfs[imem]->xfxQ2(-idNows.second, x2,
        Q2Now), PDFMINVALUE);
    }
  }

  ::LHAPDF::PDFUncertainty xfErr = pdfInfo.pdfSet.uncertainty(xfCalc);
  pdfEnvelope.centralPDF = xfErr.central;
  pdfEnvelope.errplusPDF = xfErr.errplus;
  pdfEnvelope.errminusPDF = xfErr.errminus;
  pdfEnvelope.errsymmPDF = xfErr.errsymm;
  pdfEnvelope.scalePDF = xfErr.scale;

}

//--------------------------------------------------------------------------

// Define external handles to the plugin for dynamic loading.

extern "C" LHAPDF6* newLHAPDF(int idBeamIn, string setName, int member,
                               Info* infoPtr) {
  return new LHAPDF6(idBeamIn, setName, member, 1, infoPtr);

}

extern "C" void deleteLHAPDF(LHAPDF6* pdf) {
  delete pdf;

}

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_LHAPDF6_H
