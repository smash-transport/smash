#include "include/stringprocess.h"
#include "include/angles.h"
#include "include/kinematics.h"
#include "include/random.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

namespace Smash {

// constructor
StringProcess::StringProcess() {
  //for (int ic = 0; ic < 4; ic++) {
  //  idqsetA[ic] = 0;
  //  idqsetB[ic] = 0;
  //}

  for (int imu = 0; imu < 4; imu++) {
    evecBasisAB[imu] = ThreeVector(0., 0., 0.);
  }

  PDGidA = PDGidB = 0;
  baryonA = baryonB = 0;
  //chargeA = chargeB = 0;
  PPosA = 0.;
  PNegA = 0.;
  PPosB = 0.;
  PNegB = 0.;
  massA = massB = 0.;
  sqrtsAB = 0.;
  pabscomAB = 0.;

  pLightConeMin = 0.001;
  betapowS = 0.5;

  alphapowV = 1.;
  betapowV = 2.5;

  sigmaQperp = 1.2;
  kappaString = 1.;

  reset_finalArray();
}

// destructor
StringProcess::~StringProcess() {
}

void StringProcess::set_pythia(Pythia8::Pythia *pythiaIn) {
  pythia = pythiaIn;
}

void StringProcess::reset_finalArray() {
  NpartFinal = 0;
  NpartString1 = 0;
  NpartString2 = 0;

  final_PDGid[0].resize(0);
  final_PDGid[1].resize(0);

  final_pvec[0].resize(0);
  final_pvec[1].resize(0);
  final_pvec[2].resize(0);
  final_pvec[3].resize(0);
  final_pvec[4].resize(0);

  final_tform[0].resize(0);
  final_tform[1].resize(0);
}

// compute the formation time and fill the arrays with final-state particles
int StringProcess::append_finalArray(FourVector &uString, ThreeVector &evecLong) {
  int ret;

  int nfrag;
  int ipyth;
  int ipstr, jpstr, kpstr;
  int id, islead, baryon, bstring;
  double E, px, py, pz, mass;
  double pPosTot, pNegTot;
  double tProd, xtotfac;

  bool foundFW;
  bool foundBW;

  std::vector<int> idfrag;
  std::vector<double> Efrag;
  std::vector<double> pxfrag;
  std::vector<double> pyfrag;
  std::vector<double> pzfrag;
  std::vector<double> mfrag;

  std::vector<int> indY;
  std::vector<double> pparallel;
  std::vector<double> Yparallel;
  std::vector<double> XVertexPos;
  std::vector<double> XVertexNeg;

  ThreeVector vstring;
  FourVector pvRS;
  FourVector pvCM;

  FourVector xfRS;
  FourVector xfCM;

  nfrag = 0;
  bstring = 0;
  for (ipyth = 0; ipyth < pythia->event.size(); ipyth++) {
    if (pythia->event[ipyth].isFinal()) {
      nfrag = nfrag + 1;
      id = pythia->event[ipyth].id();
      bstring = bstring + pythia->particleData.baryonNumberType(id);
    }
  }

  vstring = uString.velocity();

  idfrag.resize(nfrag);
  Efrag.resize(nfrag);
  pxfrag.resize(nfrag);
  pyfrag.resize(nfrag);
  pzfrag.resize(nfrag);
  mfrag.resize(nfrag);

  indY.resize(nfrag);
  pparallel.resize(nfrag);
  Yparallel.resize(nfrag);
  XVertexPos.resize(nfrag + 1);
  XVertexNeg.resize(nfrag + 1);

  pPosTot = 0.;
  pNegTot = 0.;
  ipstr = 0;
  for (ipyth = 0; ipyth < pythia->event.size(); ipyth++) {
    if (pythia->event[ipyth].isFinal()) {
      idfrag[ipstr] = pythia->event[ipyth].id();
      Efrag[ipstr] = pythia->event[ipyth].e();
      pxfrag[ipstr] = pythia->event[ipyth].px();
      pyfrag[ipstr] = pythia->event[ipyth].py();
      pzfrag[ipstr] = pythia->event[ipyth].pz();
      mfrag[ipstr] = pythia->event[ipyth].m();

      pparallel[ipstr] = pxfrag[ipstr] * evecLong.x1() +
                         pyfrag[ipstr] * evecLong.x2() +
                         pzfrag[ipstr] * evecLong.x3();
      Yparallel[ipstr] = 0.5 * std::log((Efrag[ipstr] + pparallel[ipstr]) /
                                   (Efrag[ipstr] - pparallel[ipstr]));

      pPosTot = pPosTot + (Efrag[ipstr] + pparallel[ipstr]) / std::sqrt(2.);
      pNegTot = pNegTot + (Efrag[ipstr] - pparallel[ipstr]) / std::sqrt(2.);

      ipstr = ipstr + 1;
    }
  }

  for (ipstr = 0; ipstr < nfrag; ipstr++) {
    kpstr = 0;
    for (jpstr = 0; jpstr < nfrag; jpstr++) {
      if ((ipstr != jpstr) && (Yparallel[ipstr] < Yparallel[jpstr])) {
        kpstr = kpstr + 1;
      }
    }
    indY[kpstr] = ipstr;
  }

  XVertexPos[0] = pPosTot / kappaString;
  for (kpstr = 0; kpstr < nfrag; kpstr++) {
    ipstr = indY[kpstr];

    XVertexPos[kpstr + 1] =
        XVertexPos[kpstr] -
        (Efrag[ipstr] + pparallel[ipstr]) / (kappaString * std::sqrt(2.));
  }

  XVertexNeg[nfrag] = pNegTot / kappaString;
  for (kpstr = nfrag - 1; kpstr >= 0; kpstr--) {
    ipstr = indY[kpstr];

    XVertexNeg[kpstr] =
        XVertexNeg[kpstr + 1] -
        (Efrag[ipstr] - pparallel[ipstr]) / (kappaString * std::sqrt(2.));
  }

  ret = 0;
  foundFW = false;
  foundBW = false;
  for (kpstr = 0; kpstr < nfrag; kpstr++) {
    ipstr = indY[kpstr];

    id = idfrag[ipstr];
    baryon = pythia->particleData.baryonNumberType(id);
    pvRS.set_x0( Efrag[ipstr] );
    pvRS.set_x1( pxfrag[ipstr] );
    pvRS.set_x2( pyfrag[ipstr] );
    pvRS.set_x3( pzfrag[ipstr] );
    mass = mfrag[ipstr];

    xfRS.set_x0( (XVertexPos[kpstr] + XVertexNeg[kpstr + 1]) / std::sqrt(2.) );
    xfRS.set_x1(
        evecLong.x1() * (XVertexPos[kpstr] - XVertexNeg[kpstr + 1]) / std::sqrt(2.) );
    xfRS.set_x2(
        evecLong.x2() * (XVertexPos[kpstr] - XVertexNeg[kpstr + 1]) / std::sqrt(2.) );
    xfRS.set_x3(
        evecLong.x3() * (XVertexPos[kpstr] - XVertexNeg[kpstr + 1]) / std::sqrt(2.) );

    tProd = (XVertexPos[kpstr] + XVertexNeg[kpstr + 1]) / std::sqrt(2.);

    if (abs(bstring) == 0) {  // mesonic string
      if ((kpstr == 0) && (foundFW == false)) {
        islead = 1;
        if (abs(baryon) == 3) {
          xtotfac = 1. / 3.;
        } else if (baryon == 0) {
          xtotfac = 0.5;
        } else {
          fprintf(stderr,
                  "  StringProcess::append_finalArray warning : particle is not "
                  "meson or baryon.\n");
        }
        foundFW = true;
      } else if ((kpstr == (nfrag - 1)) && (foundBW == false)) {
        islead = 1;
        if (abs(baryon) == 3) {
          xtotfac = 1. / 3.;
        } else if (baryon == 0) {
          xtotfac = 0.5;
        } else {
          fprintf(stderr,
                  "  StringProcess::append_finalArray warning : particle is not "
                  "meson or baryon.\n");
        }
        foundBW = true;
      } else {
        islead = 0;
        xtotfac = 0.;
      }
    } else if (abs(bstring) == 3) {  // baryonic string
      if ((baryon == bstring) && (foundFW == false)) {
        islead = 1;
        xtotfac = 2. / 3.;
        foundFW = true;
      } else if ((kpstr == (nfrag - 2)) && (baryon == 0) &&
                 (foundFW == false) && (foundBW == false)) {
        islead = 1;
        xtotfac = 0.5;
        foundBW = true;
      } else if ((kpstr == (nfrag - 1)) && (baryon == 0) && (foundFW == true) &&
                 (foundBW == false)) {
        islead = 1;
        xtotfac = 0.5;
        foundBW = true;
      } else {
        islead = 0;
        xtotfac = 0.;
      }
    } else {  // otherwise
      fprintf(stderr,
              "  StringProcess::append_finalArray warning : string is neither "
              "mesonic nor baryonic.\n");
      islead = 0;
      xtotfac = 0.;
    }

    final_PDGid[0].push_back(id);
    final_PDGid[1].push_back(islead);

    pvCM = pvRS.LorentzBoost( -vstring );
    E = pvCM.x0();
    px = pvCM.x1();
    py = pvCM.x2();
    pz = pvCM.x3();
    xfCM = xfRS.LorentzBoost( -vstring );
    tProd = xfCM.x0();

    final_pvec[0].push_back(E);
    final_pvec[1].push_back(px);
    final_pvec[2].push_back(py);
    final_pvec[3].push_back(pz);
    final_pvec[4].push_back(mass);

    final_tform[0].push_back(tProd);
    final_tform[1].push_back(xtotfac);

    ret = ret + 1;
  }

  return ret;
}

bool StringProcess::init(const ParticleList &incomingList){
  bool ret;
  //std::array<int, 3> qcontent;

  PDGcodeA = incomingList[0].pdgcode();
  PDGcodeB = incomingList[1].pdgcode();
  PDGidA = PDGcodeA.get_decimal();
  PDGidB = PDGcodeB.get_decimal();
  massA = incomingList[0].effective_mass();
  massB = incomingList[1].effective_mass();

  plabA = incomingList[0].momentum();
  plabB = incomingList[1].momentum();

  sqrtsAB = ( plabA + plabB ).abs();
  pabscomAB = pCM(sqrtsAB, massA, massB);

  make_incoming_com_momenta();
  make_orthonormal_basis();
  compute_incoming_lightcone_momenta();

  //qcontent = incomingList[0].pdgcode().quark_content();
  //idqsetA[0] = static_cast<int>( incomingList[0].pdgcode().spin_degeneracy() );
  //idqsetA[3] = qcontent[0];
  //idqsetA[2] = qcontent[1];
  //idqsetA[1] = qcontent[2];
  //qcontent = incomingList[1].pdgcode().quark_content();
  //idqsetB[0] = static_cast<int>( incomingList[1].pdgcode().spin_degeneracy() );
  //idqsetB[3] = qcontent[0];
  //idqsetB[2] = qcontent[1];
  //idqsetB[1] = qcontent[2];

  xfracMin = pLightConeMin / sqrtsAB;

  // quantum numbers of hadron A
  baryonA = 3*PDGcodeA.baryon_number();
  //chargeA = 3*PDGcodeA.charge();
  /*
  baryonA = pythia->particleData.baryonNumberType(idqsetA[1]) +
            pythia->particleData.baryonNumberType(idqsetA[2]);
  //chargeA = pythia->particleData.chargeType(idqsetA[1]) +
  //          pythia->particleData.chargeType(idqsetA[2]);
  if (idqsetA[3] != 0) {
    baryonA = baryonA + pythia->particleData.baryonNumberType(idqsetA[3]);
    //chargeA = chargeA + pythia->particleData.chargeType(idqsetA[3]);
  }
  */
  // quantum numbers of hadron B
  baryonB = 3*PDGcodeB.baryon_number();
  //chargeB = 3*PDGcodeB.charge();
  /*
  baryonB = pythia->particleData.baryonNumberType(idqsetB[1]) +
            pythia->particleData.baryonNumberType(idqsetB[2]);
  //chargeB = pythia->particleData.chargeType(idqsetB[1]) +
  //          pythia->particleData.chargeType(idqsetB[2]);
  if (idqsetB[3] != 0) {
    baryonB = baryonB + pythia->particleData.baryonNumberType(idqsetB[3]);
    //chargeB = chargeB + pythia->particleData.chargeType(idqsetB[3]);
  }
  */

  ret = true;
  return ret;
}

bool StringProcess::init_lab(PdgCode &idAIn, PdgCode &idBIn, double massAIn, double massBIn,
                       Pythia8::Vec4 plabAIn, Pythia8::Vec4 plabBIn) {
  bool ret;
  //std::array<int, 3> qcontent;

  PDGcodeA = idAIn;
  PDGcodeB = idBIn;
  PDGidA = PDGcodeA.get_decimal();
  PDGidB = PDGcodeB.get_decimal();
  massA = massAIn;
  massB = massBIn;

  plabA = FourVector(plabAIn.e(), plabAIn.px(), plabAIn.py(), plabAIn.pz());
  plabB = FourVector(plabBIn.e(), plabBIn.px(), plabBIn.py(), plabBIn.pz());

  sqrtsAB = ( plabA + plabB ).abs();
  pabscomAB = pCM(sqrtsAB, massA, massB);

  make_incoming_com_momenta();
  make_orthonormal_basis();
  compute_incoming_lightcone_momenta();

  //qcontent = idAIn.quark_content();
  //idqsetA[0] = static_cast<int>( idAIn.spin_degeneracy() );
  //idqsetA[3] = qcontent[0];
  //idqsetA[2] = qcontent[1];
  //idqsetA[1] = qcontent[2];
  //qcontent = idBIn.quark_content();
  //idqsetB[0] = static_cast<int>( idBIn.spin_degeneracy() );
  //idqsetB[3] = qcontent[0];
  //idqsetB[2] = qcontent[1];
  //idqsetB[1] = qcontent[2];

  xfracMin = pLightConeMin / sqrtsAB;

  // quantum numbers of hadron A
  baryonA = 3*PDGcodeA.baryon_number();
  //chargeA = 3*PDGcodeA.charge();
  /*
  baryonA = pythia->particleData.baryonNumberType(idqsetA[1]) +
            pythia->particleData.baryonNumberType(idqsetA[2]);
  //chargeA = pythia->particleData.chargeType(idqsetA[1]) +
  //          pythia->particleData.chargeType(idqsetA[2]);
  if (idqsetA[3] != 0) {
    baryonA = baryonA + pythia->particleData.baryonNumberType(idqsetA[3]);
    //chargeA = chargeA + pythia->particleData.chargeType(idqsetA[3]);
  }
  */
  // quantum numbers of hadron B
  baryonB = 3*PDGcodeB.baryon_number();
  //chargeB = 3*PDGcodeB.charge();
  /*
  baryonB = pythia->particleData.baryonNumberType(idqsetB[1]) +
            pythia->particleData.baryonNumberType(idqsetB[2]);
  //chargeB = pythia->particleData.chargeType(idqsetB[1]) +
  //          pythia->particleData.chargeType(idqsetB[2]);
  if (idqsetB[3] != 0) {
    baryonB = baryonB + pythia->particleData.baryonNumberType(idqsetB[3]);
    //chargeB = chargeB + pythia->particleData.chargeType(idqsetB[3]);
  }
  */

  ret = true;
  return ret;
}

bool StringProcess::init_com(PdgCode &idAIn, PdgCode &idBIn, double massAIn, double massBIn,
                       double sqrtsABIn) {
  bool ret;
  //std::array<int, 3> qcontent;

  PDGcodeA = idAIn;
  PDGcodeB = idBIn;
  PDGidA = PDGcodeA.get_decimal();
  PDGidB = PDGcodeB.get_decimal();
  massA = massAIn;
  massB = massBIn;
  sqrtsAB = sqrtsABIn;
  pabscomAB = pCM(sqrtsAB, massA, massB);

  plabA = FourVector( std::sqrt(massA*massA + pabscomAB*pabscomAB), 0., 0., pabscomAB );
  plabB = FourVector( std::sqrt(massB*massB + pabscomAB*pabscomAB), 0., 0., -pabscomAB );

  make_incoming_com_momenta();
  make_orthonormal_basis();
  compute_incoming_lightcone_momenta();

  //qcontent = idAIn.quark_content();
  //idqsetA[0] = static_cast<int>( idAIn.spin_degeneracy() );
  //idqsetA[3] = qcontent[0];
  //idqsetA[2] = qcontent[1];
  //idqsetA[1] = qcontent[2];
  //qcontent = idBIn.quark_content();
  //idqsetB[0] = static_cast<int>( idBIn.spin_degeneracy() );
  //idqsetB[3] = qcontent[0];
  //idqsetB[2] = qcontent[1];
  //idqsetB[1] = qcontent[2];

  xfracMin = pLightConeMin / sqrtsAB;

  // quantum numbers of hadron A
  baryonA = 3*PDGcodeA.baryon_number();
  //chargeA = 3*PDGcodeA.charge();
  /*
  baryonA = pythia->particleData.baryonNumberType(idqsetA[1]) +
            pythia->particleData.baryonNumberType(idqsetA[2]);
  //chargeA = pythia->particleData.chargeType(idqsetA[1]) +
  //          pythia->particleData.chargeType(idqsetA[2]);
  if (idqsetA[3] != 0) {
    baryonA = baryonA + pythia->particleData.baryonNumberType(idqsetA[3]);
    //chargeA = chargeA + pythia->particleData.chargeType(idqsetA[3]);
  }
  */
  // quantum numbers of hadron B
  baryonB = 3*PDGcodeB.baryon_number();
  //chargeB = 3*PDGcodeB.charge();
  /*
  baryonB = pythia->particleData.baryonNumberType(idqsetB[1]) +
            pythia->particleData.baryonNumberType(idqsetB[2]);
  //chargeB = pythia->particleData.chargeType(idqsetB[1]) +
  //          pythia->particleData.chargeType(idqsetB[2]);
  if (idqsetB[3] != 0) {
    baryonB = baryonB + pythia->particleData.baryonNumberType(idqsetB[3]);
    //chargeB = chargeB + pythia->particleData.chargeType(idqsetB[3]);
  }
  */

  ret = true;
  return ret;
}

// single diffractive AB > AX or XB
bool StringProcess::next_SDiff(int channel) {
  bool ret;

  int ntry;
  bool foundPabsX, foundMassX;

  int pdgidH;
  double mstrMin = 0.;
  double mstrMax = 0.;
  double pabscomHX, massH, massX, rmass;
  double QTrn, QTrx, QTry;

  int nfrag;
  int idqX1, idqX2;

  FourVector pstrHcom;
  FourVector pstrHlab;
  FourVector pstrXcom;
  FourVector pstrXlab;
  ThreeVector threeMomentum;

  FourVector ustrHcom;
  FourVector ustrHlab;
  FourVector ustrXcom;
  FourVector ustrXlab;

  double pabs;
  FourVector pnull;
  FourVector prs;
  ThreeVector evec;

  reset_finalArray();

  ntry = 0;
  foundPabsX = false;
  foundMassX = false;

  if( channel == 1 ) { // AB > AX
    mstrMin = massB;
    mstrMax = sqrtsAB - massA;
    pdgidH = PDGidA;
    massH = massA;
  } else if( channel == 2 ) { // AB > XB
    mstrMin = massA;
    mstrMax = sqrtsAB - massB;
    pdgidH = PDGidB;
    massH = massB;
  } else {
    throw std::runtime_error("invalid argument for StringProcess::next_SDiff");
  }

  while (((foundPabsX == false) || (foundMassX == false)) && (ntry < 100)) {
    ntry = ntry + 1;

    if( channel == 1 ) { // AB > AX
      makeStringEnds(PDGcodeB, idqX1, idqX2);
      //makeStringEnds(idqsetB, idqX1, idqX2);
    } else if( channel == 2 ) { // AB > XB
      makeStringEnds(PDGcodeA, idqX1, idqX2);
      //makeStringEnds(idqsetA, idqX1, idqX2);
    }

    QTrx = Random::normal(0., sigmaQperp/std::sqrt(2.) );
    QTry = Random::normal(0., sigmaQperp/std::sqrt(2.) );
    QTrn = std::sqrt(QTrx*QTrx + QTry*QTry);

    rmass = std::log(mstrMax / mstrMin) * Random::uniform(0., 1.);
    massX = mstrMin * exp(rmass);
    pabscomHX = pCM( sqrtsAB, massH, massX );

    foundPabsX = pabscomHX > QTrn;
    foundMassX = massX > (pythia->particleData.m0(idqX1) +
                          pythia->particleData.m0(idqX2));
  }

  ret = false;
  if ((foundPabsX == true) && (foundMassX == true)) {
    double sign_direction = 0.;
    if( channel == 1 ) { // AB > AX
      sign_direction = 1.;
    } else if( channel == 2 ) { // AB > XB
      sign_direction = -1.;
    }
    threeMomentum = sign_direction * (
                        evecBasisAB[3] * std::sqrt(pabscomHX*pabscomHX - QTrn*QTrn) +
                        evecBasisAB[1] * QTrx +
                        evecBasisAB[2] * QTry );
    pstrHcom = FourVector( std::sqrt(pabscomHX*pabscomHX + massH*massH), threeMomentum );
    threeMomentum = -sign_direction * (
                        evecBasisAB[3] * std::sqrt(pabscomHX*pabscomHX - QTrn*QTrn) +
                        evecBasisAB[1] * QTrx +
                        evecBasisAB[2] * QTry );
    pstrXcom = FourVector( std::sqrt(pabscomHX*pabscomHX + massX*massX), threeMomentum );

    pstrHlab = pstrHcom.LorentzBoost( -vcomAB );
    pstrXlab = pstrXcom.LorentzBoost( -vcomAB );

    ustrHcom = pstrHcom / massH;
    ustrXcom = pstrXcom / massX;
    ustrHlab = pstrHlab / massH;
    ustrXlab = pstrXlab / massX;

    threeMomentum = pstrXcom.threevec();
    pnull = FourVector( threeMomentum.abs(), threeMomentum );
    prs = pnull.LorentzBoost( ustrXcom.velocity() );
    pabs = prs.threevec().abs();
    evec = prs.threevec() / pabs;

    nfrag = fragmentString(idqX1, idqX2, massX, evec, false);
    if (nfrag > 0) {
      NpartString1 = append_finalArray(ustrXlab, evec);
    } else {
      nfrag = 0;
      NpartString1 = 0;
      ret = false;
    }

    NpartString2 = 1;
    final_PDGid[0].push_back(pdgidH);
    final_PDGid[1].push_back(1);
    final_pvec[0].push_back(pstrHlab.x0());
    final_pvec[1].push_back(pstrHlab.x1());
    final_pvec[2].push_back(pstrHlab.x2());
    final_pvec[3].push_back(pstrHlab.x3());
    final_pvec[4].push_back(massH);
    final_tform[0].push_back(0.);
    final_tform[1].push_back(1.);

    if ((NpartString1 > 0) && (NpartString2 > 0) && (nfrag == NpartString1)) {
      NpartFinal = NpartString1 + NpartString2;
      ret = true;
    }
  }

  return ret;
}

/*
// single diffractive AB > AX
bool StringProcess::next_SDiff_AX() {
  bool ret;

  int ntry;
  bool foundPabsX, foundMassX;

  double mstrMin, mstrMax;
  double pabscomAX, massX, rmass;
  double QTrn, QTrx, QTry;

  int nfrag;
  int idqX1, idqX2;

  FourVector pstrHcom;
  FourVector pstrHlab;
  FourVector pstrXcom;
  FourVector pstrXlab;
  ThreeVector threeMomentum;

  FourVector ustrHcom;
  FourVector ustrHlab;
  FourVector ustrXcom;
  FourVector ustrXlab;

  double pabs;
  FourVector pnull;
  FourVector prs;
  ThreeVector evec;

  reset_finalArray();

  ntry = 0;
  foundPabsX = false;
  foundMassX = false;
  mstrMin = massB;
  mstrMax = sqrtsAB - massA;
  if (mstrMin > mstrMax) {
    fprintf(stderr, "  StringProcess::next_SDiff_AX : mstrMin > mstrMax\n");
    fprintf(stderr,
            "  StringProcess::next_SDiff_AX : mstrMin = %e GeV, mstrMax = %e GeV\n",
            mstrMin, mstrMax);
    ntry = 100;
  }
  while (((foundPabsX == false) || (foundMassX == false)) && (ntry < 100)) {
    ntry = ntry + 1;

    makeStringEnds(PDGcodeB, idqX1, idqX2);
    //makeStringEnds(idqsetB, idqX1, idqX2);

    QTrx = Random::normal(0., sigmaQperp/std::sqrt(2.) );
    QTry = Random::normal(0., sigmaQperp/std::sqrt(2.) );
    QTrn = std::sqrt(QTrx*QTrx + QTry*QTry);

    rmass = std::log(mstrMax / mstrMin) * Random::uniform(0., 1.);
    massX = mstrMin * exp(rmass);
    pabscomAX = pCM( sqrtsAB, massA, massX );

    foundPabsX = pabscomAX > QTrn;
    foundMassX = massX > (pythia->particleData.m0(idqX1) +
                          pythia->particleData.m0(idqX2));
  }

  ret = false;
  if ((foundPabsX == true) && (foundMassX == true)) {
    threeMomentum = evecBasisAB[3] * std::sqrt(pabscomAX*pabscomAX - QTrn*QTrn) +
                        evecBasisAB[1] * QTrx +
                        evecBasisAB[2] * QTry;
    pstrHcom = FourVector( std::sqrt(pabscomAX*pabscomAX + massA*massA), threeMomentum );
    threeMomentum = -evecBasisAB[3] * std::sqrt(pabscomAX*pabscomAX - QTrn*QTrn) -
                        evecBasisAB[1] * QTrx -
                        evecBasisAB[2] * QTry;
    pstrXcom = FourVector( std::sqrt(pabscomAX*pabscomAX + massX*massX), threeMomentum );

    pstrHlab = pstrHcom.LorentzBoost( -vcomAB );
    pstrXlab = pstrXcom.LorentzBoost( -vcomAB );

    ustrHcom = pstrHcom / massA;
    ustrXcom = pstrXcom / massX;
    ustrHlab = pstrHlab / massA;
    ustrXlab = pstrXlab / massX;

    threeMomentum = pstrXcom.threevec();
    pnull = FourVector( threeMomentum.abs(), threeMomentum );
    prs = pnull.LorentzBoost( ustrXcom.velocity() );
    pabs = prs.threevec().abs();
    evec = prs.threevec() / pabs;

    nfrag = fragmentString(idqX1, idqX2, massX, evec, false);
    if (nfrag > 0) {
      NpartString1 = append_finalArray(ustrXlab, evec);
    } else {
      nfrag = 0;
      NpartString1 = 0;
      ret = false;
    }

    NpartString2 = 1;
    final_PDGid[0].push_back(PDGidA);
    final_PDGid[1].push_back(1);
    final_pvec[0].push_back(pstrHlab.x0());
    final_pvec[1].push_back(pstrHlab.x1());
    final_pvec[2].push_back(pstrHlab.x2());
    final_pvec[3].push_back(pstrHlab.x3());
    final_pvec[4].push_back(massA);
    final_tform[0].push_back(0.);
    final_tform[1].push_back(1.);

    if ((NpartString1 > 0) && (NpartString2 > 0) && (nfrag == NpartString1)) {
      NpartFinal = NpartString1 + NpartString2;
      ret = true;
    }
  }

  return ret;
}

// single diffractive AB > XB
bool StringProcess::next_SDiff_XB() {
  bool ret;

  int ntry;
  bool foundPabsX, foundMassX;

  double mstrMin, mstrMax;
  double pabscomXB, massX, rmass;
  double QTrn, QTrx, QTry;

  int nfrag;
  int idqX1, idqX2;

  FourVector pstrHcom;
  FourVector pstrHlab;
  FourVector pstrXcom;
  FourVector pstrXlab;
  ThreeVector threeMomentum;

  FourVector ustrHcom;
  FourVector ustrHlab;
  FourVector ustrXcom;
  FourVector ustrXlab;

  double pabs;
  FourVector pnull;
  FourVector prs;
  ThreeVector evec;

  reset_finalArray();

  ntry = 0;
  foundPabsX = false;
  foundMassX = false;
  mstrMin = massA;
  mstrMax = sqrtsAB - massB;
  if (mstrMin > mstrMax) {
    fprintf(stderr, "  StringProcess::next_SDiff_XB : mstrMin > mstrMax\n");
    fprintf(stderr,
            "  StringProcess::next_SDiff_XB : mstrMin = %e GeV, mstrMax = %e GeV\n",
            mstrMin, mstrMax);
    ntry = 100;
  }
  while (((foundPabsX == false) || (foundMassX == false)) && (ntry < 100)) {
    ntry = ntry + 1;

    makeStringEnds(PDGcodeA, idqX1, idqX2);
    //makeStringEnds(idqsetA, idqX1, idqX2);

    QTrx = Random::normal(0., sigmaQperp/std::sqrt(2.) );
    QTry = Random::normal(0., sigmaQperp/std::sqrt(2.) );
    QTrn = std::sqrt(QTrx*QTrx + QTry*QTry);

    rmass = log(mstrMax / mstrMin) * Random::uniform(0., 1.);
    massX = mstrMin * exp(rmass);
    pabscomXB = pCM( sqrtsAB, massX, massB );

    foundPabsX = pabscomXB > QTrn;
    foundMassX = massX > (pythia->particleData.m0(idqX1) +
                          pythia->particleData.m0(idqX2));
  }

  ret = false;
  if ((foundPabsX == true) && (foundMassX == true)) {
    threeMomentum = -evecBasisAB[3] * std::sqrt(pabscomXB*pabscomXB - QTrn*QTrn) -
                        evecBasisAB[1] * QTrx -
                        evecBasisAB[2] * QTry;
    pstrHcom = FourVector( std::sqrt(pabscomXB*pabscomXB + massB*massB), threeMomentum );
    threeMomentum = evecBasisAB[3] * std::sqrt(pabscomXB*pabscomXB - QTrn*QTrn) +
                        evecBasisAB[1] * QTrx +
                        evecBasisAB[2] * QTry;
    pstrXcom = FourVector( std::sqrt(pabscomXB*pabscomXB + massX*massX), threeMomentum );

    pstrHlab = pstrHcom.LorentzBoost( -vcomAB );
    pstrXlab = pstrXcom.LorentzBoost( -vcomAB );

    ustrHcom = pstrHcom / massB;
    ustrXcom = pstrXcom / massX;
    ustrHlab = pstrHlab / massB;
    ustrXlab = pstrXlab / massX;

    threeMomentum = pstrXcom.threevec();
    pnull = FourVector( threeMomentum.abs(), threeMomentum );
    prs = pnull.LorentzBoost( ustrXcom.velocity() );
    pabs = prs.threevec().abs();
    evec = prs.threevec() / pabs;

    nfrag = fragmentString(idqX1, idqX2, massX, evec, false);
    if (nfrag > 0) {
      NpartString1 = append_finalArray(ustrXlab, evec);
    } else {
      nfrag = 0;
      NpartString1 = 0;
      ret = false;
    }

    NpartString2 = 1;
    final_PDGid[0].push_back(PDGidB);
    final_PDGid[1].push_back(1);
    final_pvec[0].push_back(pstrHlab.x0());
    final_pvec[1].push_back(pstrHlab.x1());
    final_pvec[2].push_back(pstrHlab.x2());
    final_pvec[3].push_back(pstrHlab.x3());
    final_pvec[4].push_back(massB);
    final_tform[0].push_back(0.);
    final_tform[1].push_back(1.);

    if ((NpartString1 > 0) && (NpartString2 > 0) && (nfrag == NpartString1)) {
      NpartFinal = NpartString1 + NpartString2;
      ret = true;
    }
  }

  return ret;
}
*/

// double diffractive AB > XX
bool StringProcess::next_DDiff() {
  bool ret;

  int ntry;
  bool foundMass1, foundMass2;

  double xfracA, xfracB;
  double QPos, QNeg;
  double QTrn, QTrx, QTry;

  int nfrag1, nfrag2;
  int idq11, idq12;
  int idq21, idq22;
  double mstr1, mstr2;

  FourVector pstr1com;
  FourVector pstr1lab;
  FourVector pstr2com;
  FourVector pstr2lab;
  ThreeVector threeMomentum;

  FourVector ustr1com;
  FourVector ustr1lab;
  FourVector ustr2com;
  FourVector ustr2lab;

  double pabs;
  FourVector pnull;
  FourVector prs;
  ThreeVector evec;

  reset_finalArray();

  ntry = 0;
  foundMass1 = false;
  foundMass2 = false;
  while (((foundMass1 == false) || (foundMass2 == false)) && (ntry < 100)) {
    ntry = ntry + 1;

    makeStringEnds(PDGcodeA, idq11, idq12);
    makeStringEnds(PDGcodeB, idq21, idq22);
    //makeStringEnds(idqsetA, idq11, idq12);
    //makeStringEnds(idqsetB, idq21, idq22);

    xfracA = Random::beta_a0(xfracMin, betapowS + 1.);
    xfracB = Random::beta_a0(xfracMin, betapowS + 1.);
 
    QTrx = Random::normal(0., sigmaQperp/std::sqrt(2.) );
    QTry = Random::normal(0., sigmaQperp/std::sqrt(2.) );
    QTrn = std::sqrt(QTrx*QTrx + QTry*QTry);

    QPos = -QTrn*QTrn / (2. * xfracB * PNegB);
    QNeg = QTrn*QTrn / (2. * xfracA * PPosA);

    threeMomentum = evecBasisAB[3] * (PPosA + QPos - PNegA - QNeg) / std::sqrt(2.) +
                        evecBasisAB[1] * QTrx + evecBasisAB[2] * QTry;
    pstr1com = FourVector( (PPosA + QPos + PNegA + QNeg) / std::sqrt(2.), threeMomentum );
    mstr1 = pstr1com.sqr();

    threeMomentum = evecBasisAB[3] * (PPosB - QPos - PNegB + QNeg) / std::sqrt(2.) -
                        evecBasisAB[1] * QTrx - evecBasisAB[2] * QTry;
    pstr2com = FourVector( (PPosB - QPos + PNegB - QNeg) / std::sqrt(2.), threeMomentum );
    mstr2 = pstr2com.sqr();

    mstr1 = (mstr1 > 0.) ? std::sqrt(mstr1) : 0.;
    foundMass1 = mstr1 > (pythia->particleData.m0(idq11) +
                          pythia->particleData.m0(idq12));

    mstr2 = (mstr2 > 0.) ? std::sqrt(mstr2) : 0.;
    foundMass2 = mstr2 > (pythia->particleData.m0(idq21) +
                          pythia->particleData.m0(idq22));
  }

  ret = false;
  bool both_masses_above_pythia_threshold = foundMass1 && foundMass2;
  if ( both_masses_above_pythia_threshold ) {
    pstr1lab = pstr1com.LorentzBoost( -vcomAB );
    pstr2lab = pstr2com.LorentzBoost( -vcomAB );

    ustr1com = pstr1com / mstr1;
    ustr2com = pstr2com / mstr2;
    ustr1lab = pstr1lab / mstr1;
    ustr2lab = pstr2lab / mstr2;

    threeMomentum = pstr1com.threevec();
    pnull = FourVector( threeMomentum.abs(), threeMomentum );
    prs = pnull.LorentzBoost( ustr1com.velocity() );
    pabs = prs.threevec().abs();
    evec = prs.threevec() / pabs;

    nfrag1 = fragmentString(idq11, idq12, mstr1, evec, false);
    if (nfrag1 > 0) {
      NpartString1 = append_finalArray(ustr1lab, evec);
    } else {
      nfrag1 = 0;
      NpartString1 = 0;
      ret = false;
    }

    threeMomentum = pstr2com.threevec();
    pnull = FourVector( threeMomentum.abs(), threeMomentum );
    prs = pnull.LorentzBoost( ustr2com.velocity() );
    pabs = prs.threevec().abs();
    evec = prs.threevec() / pabs;

    nfrag2 = fragmentString(idq21, idq22, mstr2, evec, false);
    if (nfrag2 > 0) {
      NpartString2 = append_finalArray(ustr2lab, evec);
    } else {
      nfrag2 = 0;
      NpartString2 = 0;
      ret = false;
    }

    if ((NpartString1 > 0) && (NpartString2 > 0) && (nfrag1 == NpartString1) &&
        (nfrag2 == NpartString2)) {
      NpartFinal = NpartString1 + NpartString2;
      ret = true;
    }
  }

  return ret;
}

// non-diffractive
bool StringProcess::next_NDiff() {
  bool ret;

  int ntry;
  bool foundMass1, foundMass2;

  double xfracA, xfracB;
  double QPos, QNeg;
  double dPPos, dPNeg;
  double QTrn, QTrx, QTry;

  int nfrag1, nfrag2;
  int idqA1, idqA2;
  int idqB1, idqB2;
  int idq11, idq12;
  int idq21, idq22;
  double mstr1, mstr2;

  FourVector pstr1com;
  FourVector pstr1lab;
  FourVector pstr2com;
  FourVector pstr2lab;
  ThreeVector threeMomentum;

  FourVector ustr1com;
  FourVector ustr1lab;
  FourVector ustr2com;
  FourVector ustr2lab;

  double pabs;
  FourVector pnull;
  FourVector prs;
  ThreeVector evec;

  reset_finalArray();

  ntry = 0;
  foundMass1 = false;
  foundMass2 = false;
  while (((foundMass1 == false) || (foundMass2 == false)) && (ntry < 100)) {
    ntry = ntry + 1;

    makeStringEnds(PDGcodeA, idqA1, idqA2);
    makeStringEnds(PDGcodeB, idqB1, idqB2);
    //makeStringEnds(idqsetA, idqA1, idqA2);
    //makeStringEnds(idqsetB, idqB1, idqB2);

    if ((baryonA == 3) && (baryonB == 3)) {  // baryon-baryon
      idq11 = idqB1;
      idq12 = idqA2;
      idq21 = idqA1;
      idq22 = idqB2;
    } else if ((baryonA == 3) && (baryonB == 0)) {  // baryon-meson
      idq11 = idqB1;
      idq12 = idqA2;
      idq21 = idqA1;
      idq22 = idqB2;
    } else if ((baryonA == 3) && (baryonB == -3)) {  // baryon-antibaryon
      idq11 = idqB1;
      idq12 = idqA2;
      idq21 = idqA1;
      idq22 = idqB2;
    } else if ((baryonA == 0) && (baryonB == 3)) {  // meson-baryon
      idq11 = idqB1;
      idq12 = idqA2;
      idq21 = idqA1;
      idq22 = idqB2;
    } else if ((baryonA == 0) && (baryonB == 0)) {  // meson-meson
      idq11 = idqB1;
      idq12 = idqA2;
      idq21 = idqA1;
      idq22 = idqB2;
    } else if ((baryonA == 0) && (baryonB == -3)) {  // meson-antibaryon
      idq11 = idqA1;
      idq12 = idqB2;
      idq21 = idqB1;
      idq22 = idqA2;
    } else if ((baryonA == -3) && (baryonB == 3)) {  // antibaryon-baryon
      idq11 = idqA1;
      idq12 = idqB2;
      idq21 = idqB1;
      idq22 = idqA2;
    } else if ((baryonA == -3) && (baryonB == 0)) {  // antibaryon-meson
      idq11 = idqA1;
      idq12 = idqB2;
      idq21 = idqB1;
      idq22 = idqA2;
    } else if ((baryonA == -3) && (baryonB == -3)) {  // antibaryon-antibaryon
      idq11 = idqA1;
      idq12 = idqB2;
      idq21 = idqB1;
      idq22 = idqA2;
    } else {
      fprintf(stderr,
              "  StringProcess::next_NDiff : incorrect baryon number of incoming "
              "hadrons.\n");
      fprintf(stderr, "  StringProcess::next_NDiff : baryonA = %d, baryonB = %d\n",
              baryonA, baryonB);
      exit(1);
    }

    xfracA = Random::beta(alphapowV, betapowV);
    xfracB = Random::beta(alphapowV, betapowV);

    QTrx = Random::normal(0., sigmaQperp/std::sqrt(2.) );
    QTry = Random::normal(0., sigmaQperp/std::sqrt(2.) );
    QTrn = std::sqrt(QTrx*QTrx + QTry*QTry);

    QPos = -QTrn*QTrn / (2. * xfracB * PNegB);
    QNeg = QTrn*QTrn / (2. * xfracA * PPosA);

    dPPos = -xfracA * PPosA - QPos;
    dPNeg = xfracB * PNegB - QNeg;

    threeMomentum = evecBasisAB[3] * (PPosA + dPPos - PNegA - dPNeg) / std::sqrt(2.) +
                        evecBasisAB[1] * QTrx + evecBasisAB[2] * QTry;
    pstr1com = FourVector( (PPosA + dPPos + PNegA + dPNeg) / std::sqrt(2.), threeMomentum );
    mstr1 = pstr1com.sqr();

    threeMomentum = evecBasisAB[3] * (PPosB - dPPos - PNegB + dPNeg) / std::sqrt(2.) -
                        evecBasisAB[1] * QTrx - evecBasisAB[2] * QTry;
    pstr2com = FourVector( (PPosB - dPPos + PNegB - dPNeg) / std::sqrt(2.), threeMomentum );
    mstr2 = pstr2com.sqr();

    mstr1 = (mstr1 > 0.) ? std::sqrt(mstr1) : 0.;
    foundMass1 = mstr1 > (pythia->particleData.m0(idq11) +
                          pythia->particleData.m0(idq12));

    mstr2 = (mstr2 > 0.) ? std::sqrt(mstr2) : 0.;
    foundMass2 = mstr2 > (pythia->particleData.m0(idq21) +
                          pythia->particleData.m0(idq22));
  }

  ret = false;
  bool both_masses_above_pythia_threshold = foundMass1 && foundMass2;
  if ( both_masses_above_pythia_threshold ) {
    pstr1lab = pstr1com.LorentzBoost( -vcomAB );
    pstr2lab = pstr2com.LorentzBoost( -vcomAB );

    ustr1com = pstr1com / mstr1;
    ustr2com = pstr2com / mstr2;
    ustr1lab = pstr1lab / mstr1;
    ustr2lab = pstr2lab / mstr2;

    threeMomentum = pstr1com.threevec();
    pnull = FourVector( threeMomentum.abs(), threeMomentum );
    prs = pnull.LorentzBoost( ustr1com.velocity() );
    pabs = prs.threevec().abs();
    evec = prs.threevec() / pabs;

    nfrag1 = fragmentString(idq11, idq12, mstr1, evec, false);
    if (nfrag1 > 0) {
      NpartString1 = append_finalArray(ustr1lab, evec);
    } else {
      nfrag1 = 0;
      NpartString1 = 0;
      ret = false;
    }

    threeMomentum = pstr2com.threevec();
    pnull = FourVector( threeMomentum.abs(), threeMomentum );
    prs = pnull.LorentzBoost( ustr2com.velocity() );
    pabs = prs.threevec().abs();
    evec = prs.threevec() / pabs;

    nfrag2 = fragmentString(idq21, idq22, mstr2, evec, false);
    if (nfrag2 > 0) {
      NpartString2 = append_finalArray(ustr2lab, evec);
    } else {
      nfrag2 = 0;
      NpartString2 = 0;
      ret = false;
    }

    if ((NpartString1 > 0) && (NpartString2 > 0) && (nfrag1 == NpartString1) &&
        (nfrag2 == NpartString2)) {
      NpartFinal = NpartString1 + NpartString2;
      ret = true;
    }
  }

  return ret;
}

// baryon-antibaryon annihilation
bool StringProcess::next_BBbarAnn(){
	bool ret;

	int ntry;
	bool isBBbarpair, isAnnihilating;

	int ic, jc;
	int ijc, ipr, npr;
	std::array<int,3> quark_content_A;
	std::array<int,3> quark_content_B;
	std::vector<int> indexAnn;

	int idq11, idq12, idq12prev;
	int idq21, idq22, idq22prev;
	int nfrag1, nfrag2;
	double mstr1, mstr2;
	double mstr1Min, mstr2Min;

	FourVector ustr1lab;
	FourVector ustr2lab;

	double pabs;
	ThreeVector evec;

	ustr1lab = ucomAB;
	ustr2lab = ucomAB;

	indexAnn.resize(0);

	reset_finalArray();

	quark_content_A = PDGcodeA.quark_content();
	quark_content_B = PDGcodeB.quark_content();

	isBBbarpair = ( (baryonA == 3) && (baryonB == -3) ) || ( (baryonA == -3) && (baryonB == 3) );
	isAnnihilating = false;

	// if it is baryon-antibaryon pair
	if( isBBbarpair == true ){ // if it is
		mstr1 = 0.5*sqrtsAB;
		mstr2 = 0.5*sqrtsAB;

		for(ic = 0; ic < 3; ic++){
			for(jc = 0; jc < 3; jc++){
				if( quark_content_A[ic] == -quark_content_B[jc] ){
					ijc = ic*10 + jc;
					indexAnn.push_back( ijc );
				}
			}
		}
		//for(ic = 1; ic <= 3; ic++){
		//	for(jc = 1; jc <= 3; jc++){
		//		if( idqsetA[ic] == -idqsetB[jc] ){
		//			ijc = ic*10 + jc;
		//			indexAnn.push_back( ijc );
		//		}
		//	}
		//}
		npr = indexAnn.size();
		fprintf(stderr,"  StringProcess::next_BBarAnn : %d possible pairs for qqbar annihilation\n", npr);
		/* if it is a BBbar pair but there is no qqbar pair to annihilate,
		 * nothing happens */
		if( npr == 0 ){
			NpartString1 = 1;
			final_PDGid[0].push_back(PDGidA);
			final_PDGid[1].push_back(1);
			final_pvec[0].push_back(plabA.x0());
			final_pvec[1].push_back(plabA.x1());
			final_pvec[2].push_back(plabA.x2());
			final_pvec[3].push_back(plabA.x3());
			final_pvec[4].push_back(massA);
			final_tform[0].push_back(0.);
			final_tform[1].push_back(1.);

			NpartString2 = 1;
			final_PDGid[0].push_back(PDGidB);
			final_PDGid[1].push_back(1);
			final_pvec[0].push_back(plabB.x0());
			final_pvec[1].push_back(plabB.x1());
			final_pvec[2].push_back(plabB.x2());
			final_pvec[3].push_back(plabB.x3());
			final_pvec[4].push_back(massB);
			final_tform[0].push_back(0.);
			final_tform[1].push_back(1.);

			isAnnihilating = false;

			NpartFinal = NpartString1 + NpartString2;
			ret = true;
		}// endif no qqbar pair to annihilate

		ntry = 0;
		while( ( npr > 0 ) && ( isAnnihilating == false ) && ( ntry < 100 ) ){
			ntry = ntry + 1;

			// randomly choose a qqbar pair to annihilate
			ipr = Random::uniform_int(0, npr - 1);
			ijc = indexAnn.at(ipr);
			ic = ( ijc - (ijc%10) )/10;
			jc = ijc%10;
			fprintf(stderr,"  StringProcess::next_BBarAnn : ic = %d, jc = %d chosen\n", ic, jc);
			// make two qqbar pairs to excite strings
			if( (baryonA == 3) && (baryonB == -3) ){
				idq11 = quark_content_A[(ic + 1)%3];
				idq12 = quark_content_B[(jc + 1)%3];
				idq21 = quark_content_A[(ic + 2)%3];
				idq22 = quark_content_B[(jc + 2)%3];
				//idq11 = idqsetA[1 + ic%3];
				//idq12 = idqsetB[1 + jc%3];
				//idq21 = idqsetA[1 + (ic + 1)%3];
				//idq22 = idqsetB[1 + (jc + 1)%3];
			}
			else if( (baryonA == -3) && (baryonB == 3) ){
				idq11 = quark_content_B[(ic + 1)%3];
				idq12 = quark_content_A[(jc + 1)%3];
				idq21 = quark_content_B[(ic + 2)%3];
				idq22 = quark_content_A[(jc + 2)%3];
				//idq11 = idqsetB[1 + jc%3];
				//idq12 = idqsetA[1 + ic%3];
				//idq21 = idqsetB[1 + (jc + 1)%3];
				//idq22 = idqsetA[1 + (ic + 1)%3];
			}
			// randomly choose if we flip the antiquark contents
			if( Random::uniform_int(0, 1) == 0 ){
				idq12prev = idq12;
				idq22prev = idq22;
				idq12 = idq22prev;
				idq22 = idq12prev;
			}
			fprintf(stderr,"  StringProcess::next_BBarAnn : string 1 with %d, %d\n", idq11, idq12);
			fprintf(stderr,"  StringProcess::next_BBarAnn : string 2 with %d, %d\n", idq21, idq22);

			mstr1Min = pythia->particleData.m0(idq11) + pythia->particleData.m0(idq12);
			mstr2Min = pythia->particleData.m0(idq21) + pythia->particleData.m0(idq22);
			isAnnihilating = ( mstr1 > mstr1Min ) && ( mstr2 > mstr2Min );
		}
	}
	else{ // if it is not
		fprintf(stderr,"  StringProcess::next_BBarAnn failure : it is not BBbar pair.\n");
		isAnnihilating = false;
	}
	// endif baryon-antibaryon pair

	// implement collision in the case of annihilating BBbar pair
	if( isAnnihilating == true ){
		ret = false;

		// string 1
		pabs = pcomA.threevec().abs();
		evec = pcomA.threevec() / pabs;

		nfrag1 = fragmentString(idq11, idq12, mstr1, evec, false);
		if( nfrag1 > 0 ){
			NpartString1 = append_finalArray(ustr1lab, evec);
		}
		else{
			nfrag1 = 0;
			NpartString1 = 0;
			ret = false;
		}

		// string 2
		pabs = pcomB.threevec().abs();
		evec = pcomB.threevec() / pabs;

		nfrag2 = fragmentString(idq21, idq22, mstr2, evec, false);
		if( nfrag2 > 0 ){
			NpartString2 = append_finalArray(ustr2lab, evec);
		}
		else{
			nfrag2 = 0;
			NpartString2 = 0;
			ret = false;
		}

		if( ( NpartString1 > 0 ) && ( NpartString2 > 0 )
			&& ( nfrag1 == NpartString1 ) && ( nfrag2 == NpartString2 ) ){
			NpartFinal = NpartString1 + NpartString2;
			ret = true;
		}
	}

	return ret;
}

void StringProcess::make_incoming_com_momenta(){
  ucomAB = ( plabA + plabB )/sqrtsAB;
  vcomAB = ucomAB.velocity();

  pcomA = plabA.LorentzBoost(vcomAB);
  pcomB = plabB.LorentzBoost(vcomAB);
}

void StringProcess::make_orthonormal_basis(){
  if (std::abs(pcomA.x3()) < (1. - 1.0e-8) * pabscomAB) {
    double ex, ey, et;
    double theta, phi;

    evecBasisAB[3] = pcomA.threevec() / pabscomAB;

    theta = std::acos(evecBasisAB[3].x3());

    ex = evecBasisAB[3].x1();
    ey = evecBasisAB[3].x2();
    et = std::sqrt(ex*ex + ey*ey);
    if (ey > 0.) {
      phi = std::acos(ex / et);
    } else {
      phi = -std::acos(ex / et);
    }

    evecBasisAB[1].set_x1( cos(theta) * cos(phi) );
    evecBasisAB[1].set_x2( cos(theta) * sin(phi) );
    evecBasisAB[1].set_x3( -sin(theta) );

    evecBasisAB[2].set_x1( -sin(phi) );
    evecBasisAB[2].set_x2( cos(phi) );
    evecBasisAB[2].set_x3( 0. );
  } else {
    if (pcomA.x3() > 0.) {
      evecBasisAB[1] = ThreeVector(1., 0., 0.);
      evecBasisAB[2] = ThreeVector(0., 1., 0.);
      evecBasisAB[3] = ThreeVector(0., 0., 1.);
    } else {
      evecBasisAB[1] = ThreeVector(0., 1., 0.);
      evecBasisAB[2] = ThreeVector(1., 0., 0.);
      evecBasisAB[3] = ThreeVector(0., 0., -1.);
    }
  }
}

void StringProcess::compute_incoming_lightcone_momenta(){
  PPosA = ( pcomA.x0() + evecBasisAB[3] * pcomA.threevec() ) / std::sqrt(2.);
  PNegA = ( pcomA.x0() - evecBasisAB[3] * pcomA.threevec() ) / std::sqrt(2.);
  PPosB = ( pcomB.x0() + evecBasisAB[3] * pcomB.threevec() ) / std::sqrt(2.);
  PNegB = ( pcomB.x0() - evecBasisAB[3] * pcomB.threevec() ) / std::sqrt(2.);
}

void StringProcess::makeStringEnds(PdgCode &pdgcodeIn, int &idq1, int &idq2){
  int ir, ic, jc;
  int idq1tmp, idq2tmp;
  std::array<int,3> qcontent;
  std::array<int,2> idqtmp;

  qcontent = pdgcodeIn.quark_content();

  // if it is meson/baryon
  if (qcontent[0] == 0) {  // meson
    ir = 1 + Random::uniform_int(0, 1);

    idq1tmp = qcontent[ir];
    jc = 1 + ir % 2;
    idq2tmp = qcontent[jc];
  } else {  // baryon
    ir = Random::uniform_int(0, 2);

    idq1tmp = qcontent[ir];
    for (ic = 0; ic < 2; ic++) {
      jc = (ir + ic + 1) % 3;
      idqtmp[ic] = std::abs(qcontent[jc]);
    }

    if (idqtmp[0] == idqtmp[1]) {
      idq2tmp = idqtmp[0] * 1000 + idqtmp[1] * 100 + 3;
    } else {
      if (idqtmp[0] > idqtmp[1]) {
        idq2tmp = idqtmp[0] * 1000 + idqtmp[1] * 100;
      } else {
        idq2tmp = idqtmp[1] * 1000 + idqtmp[0] * 100;
      }

      double rspin = Random::uniform(0., 1.);
      if ( rspin < 0.25 ) {
        idq2tmp = idq2tmp + 1;
      } else {
        idq2tmp = idq2tmp + 3;
      }
    }

    if (idq1tmp < 0) {
      idq2tmp = -idq2tmp;
    }
  }  // endif meson/baryon

  if (idq1tmp > 0) {
    idq1 = idq1tmp;
    idq2 = idq2tmp;
  } else {
    idq1 = idq2tmp;
    idq2 = idq1tmp;
  }

  /* some mesons with PDG id 11X are actually mixed state of uubar and ddbar.
   * have a random selection whether we have uubar or ddbar in this case. */
  if ( (qcontent[0] == 0) && (idq1 == 1) && (idq2 == -1) ) {
    if ( Random::uniform_int(0, 1) == 0 ) {
      idq1 = 2;
      idq2 = -2;
    }
  }
}

/*
void StringProcess::makeStringEnds(std::array<int,4> &idqset, int &idq1, int &idq2) {
  int ir, ic, jc;
  int idq1tmp, idq2tmp;
  std::array<int,3> idqtmp;

  double rspin;

  // if it is meson/baryon
  if (idqset[3] == 0) {  // meson
    ir = 1 + Random::uniform_int(0, 1);

    idq1tmp = idqset[ir];
    jc = 1 + ir % 2;
    idq2tmp = idqset[jc];
  } else {  // baryon
    ir = 1 + Random::uniform_int(0, 2);
    idqtmp[0] = idqset[0];

    idq1tmp = idqset[ir];
    for (ic = 0; ic < 2; ic++) {
      jc = 1 + (ir + ic) % 3;
      idqtmp[ic + 1] = std::abs(idqset[jc]);
    }

    if (idqtmp[1] == idqtmp[2]) {
      idq2tmp = idqtmp[1] * 1000 + idqtmp[2] * 100 + 3;
    } else {
      if (idqtmp[1] > idqtmp[2]) {
        idq2tmp = idqtmp[1] * 1000 + idqtmp[2] * 100;
      } else {
        idq2tmp = idqtmp[2] * 1000 + idqtmp[1] * 100;
      }

      rspin = Random::uniform(0., 1.);
      if ( rspin < 0.25 ) {
        idq2tmp = idq2tmp + 1;
      } else {
        idq2tmp = idq2tmp + 3;
      }
    }

    if (idq1tmp < 0) {
      idq2tmp = -idq2tmp;
    }
  }  // endif meson/baryon

  if (idq1tmp > 0) {
    idq1 = idq1tmp;
    idq2 = idq2tmp;
  } else {
    idq1 = idq2tmp;
    idq2 = idq1tmp;
  }

  // some mesons with PDG id 11X are actually mixed state of uubar and ddbar.
  // have a random selection whether we have uubar or ddbar in this case.
  if ((idqset[3] == 0) && (idq1 == 1) && (idq2 == -1)) {
    if ( Random::uniform_int(0, 1) == 0 ) {
      idq1 = 2;
      idq2 = -2;
    }
  }
}
*/

int StringProcess::fragmentString(int idq1, int idq2, double mString,
                            ThreeVector &evecLong, bool random_rotation) {
  int number_of_fragments;
  bool successful_hadronization;

  int bstring;
  int ipart;
  int status;
  int col, acol;
  double sign_direction;
  double pCMquark;
  double m1, m2;
  Pythia8::Vec4 pquark;

  ThreeVector p3vec;
  FourVector pvRS;

  pythia->event.reset();

  bstring = pythia->particleData.baryonNumberType(idq1) +
            pythia->particleData.baryonNumberType(idq2);
  if( bstring == -3 ){
    sign_direction = -1;
  }
  else{
    sign_direction = 1.;
  }

  m1 = pythia->particleData.m0(idq1);
  m2 = pythia->particleData.m0(idq2);
  pCMquark = pCM( mString, m1, m2 );

  if (random_rotation == true) {
    Angles phitheta;
    phitheta.distribute_isotropically();

    p3vec.set_x1( pCMquark * phitheta.threevec().x1() );
    p3vec.set_x2( pCMquark * phitheta.threevec().x2() );
    p3vec.set_x3( pCMquark * phitheta.threevec().x3() );
  } else {
    if ( Random::uniform_int(0, 1) == 0 ) {
      evecLong.set_x1( -evecLong.x1() );
      evecLong.set_x2( -evecLong.x2() );
      evecLong.set_x3( -evecLong.x3() );
    }
  }

  if (m1 + m2 > mString) {
    throw std::runtime_error("String fragmentation: m1 + m2 > mString");
  }

  status = 1;
  col = 1;
  acol = 0;
  if (random_rotation == true) {
    pvRS = FourVector(0., -p3vec);
  } else {
    pvRS = FourVector(0., -sign_direction * pCMquark * evecLong);
  }
  pvRS.set_x0( std::sqrt(m1*m1 + pCMquark*pCMquark) );

  pquark.e( pvRS.x0() );
  pquark.px( pvRS.x1() );
  pquark.py( pvRS.x2() );
  pquark.pz( pvRS.x3() );

  pythia->event.append(idq1, status, col, acol, pquark, m1);

  status = 1;
  col = 0;
  acol = 1;
  if (random_rotation == true) {
    pvRS = FourVector(0., p3vec);
  } else {
    pvRS = FourVector(0., sign_direction * pCMquark * evecLong);
  }
  pvRS.set_x0( std::sqrt(m2*m2 + pCMquark*pCMquark) );

  pquark.e( pvRS.x0() );
  pquark.px( pvRS.x1() );
  pquark.py( pvRS.x2() );
  pquark.pz( pvRS.x3() );

  pythia->event.append(idq2, status, col, acol, pquark, m2);

  successful_hadronization = pythia->forceHadronLevel();
  number_of_fragments = 0;
  if (successful_hadronization == true) {
    for (ipart = 0; ipart < pythia->event.size(); ipart++) {
      if (pythia->event[ipart].isFinal()) {
        number_of_fragments++;
      }
    }
  }

  return number_of_fragments;
}

}  // namespace Smash
