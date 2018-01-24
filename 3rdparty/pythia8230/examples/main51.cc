// main51.cc is a part of the PYTHIA event generator.
// Copyright (C) 2017 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Test of LHAPDF interface and whether PDF's behave sensibly.
// October 2017: updated to test external LHAPDF6 vs internal LHAGrid1.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

// Integration to check momentum sum rule.

double integrate(PDF* nowPDF, double Q2) {

  // Number of points, x ranges and initial values.
  int    nLin  = 980;
  int    nLog  = 1000;
  double xLin  = 0.02;
  double xLog  = 1e-8;
  double dxLin = (1. - xLin) / nLin;
  double dxLog = log(xLin / xLog) / nLog;
  double sum   = 0.;
  double x, sumNow;

  // Integration at large x in linear steps.
  for (int iLin = 0; iLin < nLin; ++iLin) {
    x      = xLin + (iLin + 0.5) * dxLin;
    sumNow = nowPDF->xf( 21, x, Q2);
    for (int i = 1; i < 6; ++i)
      sumNow += nowPDF->xf( i, x, Q2) + nowPDF->xf( -i, x, Q2);
    sum   += dxLin * sumNow;
  }

  // Integration at small x in logarithmic steps.
  for (int iLog = 0; iLog < nLog; ++iLog) {
    x      = xLog * pow( xLin / xLog, (iLog + 0.5) / nLog );
    sumNow = nowPDF->xf( 21, x, Q2);
    for (int i = 1; i < 6; ++i)
      sumNow += nowPDF->xf( i, x, Q2) + nowPDF->xf( -i, x, Q2);
    sum   += dxLog * x * sumNow;
  }

  // Done.
  return sum;

}

//==========================================================================

int main() {

  // Chosen new PDF set; LHAPDF5 file name conventions. Usage deprecated.
  //string pdfSet = "LHAPDF5:cteq5l.LHgrid";
  //string pdfSet = "LHAPDF5:cteq61.LHpdf";
  //string pdfSet = "LHAPDF5:cteq61.LHgrid";
  //string pdfSet = "LHAPDF5:MRST2004nlo.LHgrid";
  //string pdfSet = "LHAPDF5:MRST2001lo.LHgrid";
  //string pdfSet = "LHAPDF5:MRST2007lomod.LHgrid";
  // Initialize some external PDF and compare it with internal CTEQ 5L.
  //Info info;
  //PDF* extPDF = new LHAPDF(2212, pdfSet, &info);
  //PDF* intPDF = new CTEQ5L(2212);

  // Chosen new PDF set; LHAPDF6 file name conventions.
  // By default the internal .dat file is in share/Pythia/xmldoc
  // and the external in share/LHAPDF of the respective library.
  //string pdfSet = "NNPDF31_lo_as_0130";
  // The next one does not yet have official LHAPDF numbers,
  // so it has to be added by hand to share/LHAPDF/pdfsets.
  string pdfSet = "mcpdf_test_replicas";

  // Pointers to external LHAPDF6 and internal LHAGrid1 PDF packages,
  // for the same PDF set, specifically the central member.
  Info info;
  PDF* extPDF = new LHAPDF(2212, "LHAPDF6:" + pdfSet, &info);
  PDF* intPDF = new LHAGrid1(2212, pdfSet + "_0000.dat");

  // Alternative: compare two Pomeron PDF's. Boost second by factor 2.
  //PDF* extPDF = new PomFix( 990, -0.2, 2.5, 0., 3., 0.4, 0.5);
  //PDF* intPDF = new PomH1Jets( 990, 2.);
  //PDF* extPDF = new PomH1FitAB( 990, 2);
  //PDF* intPDF = new PomH1FitAB( 990, 3);

  // Allow extrapolation of PDF's beyond x and Q2 boundaries, at own risk.
  // Default behaviour is to freeze PDF's at boundaries.
  intPDF->setExtrapolate(true);
  extPDF->setExtrapolate(true);

  // Histogram F(x, Q2) = (9/4) x*g(x, Q2) + sum_{i = q, qbar} x*f_i(x, Q2)
  // for range 10^{-8} < x < 1 logarithmic in x and for Q2 = 4 and 100.
  Hist extF4("F( x, Q2 = 4) external", 80 , -8., 0.);
  Hist intF4("F( x, Q2 = 4) internal", 80 , -8., 0.);
  Hist ratF4("F( x, Q2 = 4) internal/external", 80 , -8., 0.);
  Hist extF100("F( x, Q2 = 100) external", 80 , -8., 0.);
  Hist intF100("F( x, Q2 = 100) internal", 80 , -8., 0.);
  Hist ratF100("F( x, Q2 = 100) internal/external", 80 , -8., 0.);

  // Loop over the two Q2 values.
  for (int iQ = 0; iQ < 2; ++iQ) {
    double Q2 = (iQ == 0) ? 4. : 100;

    // Loop over x values, in a logarithmic scale.
    for (int iX = 0; iX < 80; ++iX) {
      double xLog = -(0.1 * iX + 0.05);
      double x = pow( 10., xLog);

      // Evaluate external summed PDF.
      double extSum = 2.25 * extPDF->xf( 21, x, Q2);
      for (int i = 1; i < 6; ++i)
        extSum += extPDF->xf( i, x, Q2) + extPDF->xf( -i, x, Q2);
      if (iQ == 0) extF4.fill ( xLog, extSum );
      else       extF100.fill ( xLog, extSum );

      // Evaluate internal summed PDF.
      double intSum = 2.25 * intPDF->xf( 21, x, Q2);
      for (int i = 1; i < 6; ++i)
        intSum += intPDF->xf( i, x, Q2) + intPDF->xf( -i, x, Q2);
      if (iQ == 0) intF4.fill ( xLog, intSum );
      else       intF100.fill ( xLog, intSum );

    // End loops over x and Q2 values.
    }
  }

  // Show F(x, Q2) and their ratio internal/external.
  ratF4 = intF4 / extF4;
  ratF100 = intF100 / extF100;
  cout << extF4 << intF4 << ratF4 << extF100 << intF100 << ratF100;

  // Histogram momentum sum as a function of Q2 (or rather log10(Q2)).
  Hist extXSum("momentum sum(log10(Q2)) external", 100, -2., 8.);
  Hist intXSum("momentum sum(log10(Q2)) internal", 100, -2., 8.);
  Hist ratXSum("momentum sum(log10(Q2)) internal/external", 100, -2., 8.);

  // Loop over Q2 values.
  for (int iQ = 0; iQ < 100; ++iQ) {
    double log10Q2 = -2.0 + 0.1 * iQ + 0.05;
    double Q2 = pow( 10., log10Q2);

    // Evaluate external and internal momentum sums.
    double extSum = integrate( extPDF, Q2);
    extXSum.fill( log10Q2, extSum);
    double intSum = integrate( intPDF, Q2);
    intXSum.fill( log10Q2, intSum);
  }

  // Show momentum sum as a function of Q2.
  ratXSum = intXSum / extXSum;
  cout << extXSum << intXSum << ratXSum;

  // Done.
  delete extPDF;
  delete intPDF;
  return 0;
}
