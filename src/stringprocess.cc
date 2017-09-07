#include "include/stringprocess.h"
#include "include/random.h"
#include "include/kinematics.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

namespace Smash {

// constructor
StringProcess::StringProcess() {
  int ic, imu;
  //idqsetA = static_cast<int *>(malloc(4 * sizeof(int)));
  //idqsetB = static_cast<int *>(malloc(4 * sizeof(int)));
  for (ic = 0; ic < 4; ic++) {
    idqsetA[ic] = 0;
    idqsetB[ic] = 0;
  }

  //plabA = static_cast<double *>(malloc(4 * sizeof(double)));
  //plabB = static_cast<double *>(malloc(4 * sizeof(double)));
  //pcomA = static_cast<double *>(malloc(4 * sizeof(double)));
  //pcomB = static_cast<double *>(malloc(4 * sizeof(double)));
  //ucomAB = static_cast<double *>(malloc(4 * sizeof(double)));
  //evecBasisAB = new ThreeVector[4];
  //evecBasisAB = static_cast<double **>(malloc(4 * sizeof(double *)));
  for (imu = 0; imu < 4; imu++) {
    //plabA[imu] = 0.;
    //plabB[imu] = 0.;
    //pcomA[imu] = 0.;
    //pcomB[imu] = 0.;
    //ucomAB[imu] = 0.;
    evecBasisAB[imu] = ThreeVector(0., 0., 0.);
    //evecBasisAB[imu].set_x1( 0. );
    //evecBasisAB[imu].set_x2( 0. );
    //evecBasisAB[imu].set_x3( 0. );
    //evecBasisAB[imu] = static_cast<double *>(malloc(4 * sizeof(double)));
    //for (inu = 0; inu < 4; inu++) {
    //  evecBasisAB[imu][inu] = 0.;
    //}
  }

  //XSecSummed = static_cast<double *>(malloc(5 * sizeof(double)));
  //for (iproc = 0; iproc < 5; iproc++) {
  //  XSecSummed[iproc] = 0.;
  //}

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

  //lorentz = new Lorentz();
  //final_PDGid = new vector<int>[2];
  //final_pvec = new vector<double>[5];
  //final_tform = new vector<double>[2];
  reset_finalArray();
}

// destructor
StringProcess::~StringProcess() {
  //int imu;

  //free(idqsetA);
  //free(idqsetB);

  //free(plabA);
  //free(plabB);
  //free(pcomA);
  //free(pcomB);
  //free(ucomAB);
  //delete[] evecBasisAB;
  //for (imu = 0; imu < 4; imu++) {
  //  free(evecBasisAB[imu]);
  //}
  //free(evecBasisAB);

  //free(XSecSummed);

  //delete lorentz;

  //delete[] final_PDGid;
  //delete[] final_pvec;
  //delete[] final_tform;
}

void StringProcess::set_pythia(Pythia8::Pythia *pythiaIn) {
  pythia = pythiaIn;

  //sigmaTot.init(&pythia->info, pythia->settings, &pythia->particleData);
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
  //int *idfrag;
  //double *Efrag;
  //double *pxfrag;
  //double *pyfrag;
  //double *pzfrag;
  //double *mfrag;

  std::vector<int> indY;
  std::vector<double> pparallel;
  std::vector<double> Yparallel;
  std::vector<double> XVertexPos;
  std::vector<double> XVertexNeg;
  //int *indY;
  // double zparallel;
  //double *pparallel;
  //double *Yparallel;
  //double *XVertexPos;
  //double *XVertexNeg;

  ThreeVector vstring;
  FourVector pvRS;
  FourVector pvCM;
  //double *pvRS;
  //double *pvCM;

  FourVector xfRS;
  FourVector xfCM;
  //double *xfRS;
  //double *xfCM;

  //pvRS = static_cast<double *>(malloc(4 * sizeof(double)));
  //pvCM = static_cast<double *>(malloc(4 * sizeof(double)));

  //xfRS = static_cast<double *>(malloc(4 * sizeof(double)));
  //xfCM = static_cast<double *>(malloc(4 * sizeof(double)));

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

  // fprintf(stderr,"  StringProcess::append_finalArray nfrag = %d\n", nfrag);

  idfrag.resize(nfrag);
  Efrag.resize(nfrag);
  pxfrag.resize(nfrag);
  pyfrag.resize(nfrag);
  pzfrag.resize(nfrag);
  mfrag.resize(nfrag);
  //idfrag = static_cast<int *>(malloc(nfrag * sizeof(int)));
  //Efrag = static_cast<double *>(malloc(nfrag * sizeof(double)));
  //pxfrag = static_cast<double *>(malloc(nfrag * sizeof(double)));
  //pyfrag = static_cast<double *>(malloc(nfrag * sizeof(double)));
  //pzfrag = static_cast<double *>(malloc(nfrag * sizeof(double)));
  //mfrag = static_cast<double *>(malloc(nfrag * sizeof(double)));

  indY.resize(nfrag);
  pparallel.resize(nfrag);
  Yparallel.resize(nfrag);
  XVertexPos.resize(nfrag + 1);
  XVertexNeg.resize(nfrag + 1);
  //indY = static_cast<int *>(malloc(nfrag * sizeof(int)));
  //pparallel = static_cast<double *>(malloc(nfrag * sizeof(double)));
  //Yparallel = static_cast<double *>(malloc(nfrag * sizeof(double)));
  //XVertexPos = static_cast<double *>(malloc((nfrag + 1) * sizeof(double)));
  //XVertexNeg = static_cast<double *>(malloc((nfrag + 1) * sizeof(double)));

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
      //if (evecLong == NULL) {
      //  pparallel[ipstr] = pzfrag[ipstr];
      //  // fprintf(stderr,"  StringProcess::append_finalArray : NULL evecLong\n");
      //} else {
      //  pparallel[ipstr] = pxfrag[ipstr] * evecLong[1] +
      //                     pyfrag[ipstr] * evecLong[2] +
      //                     pzfrag[ipstr] * evecLong[3];
      //  // fprintf(stderr,"  StringProcess::append_finalArray : non-zero evecLong\n");
      //}
      Yparallel[ipstr] = 0.5 * std::log((Efrag[ipstr] + pparallel[ipstr]) /
                                   (Efrag[ipstr] - pparallel[ipstr]));

      pPosTot = pPosTot + (Efrag[ipstr] + pparallel[ipstr]) / std::sqrt(2.);
      pNegTot = pNegTot + (Efrag[ipstr] - pparallel[ipstr]) / std::sqrt(2.);

      ipstr = ipstr + 1;
    }
  }
  // fprintf(stderr,"  StringProcess::append_finalArray pPosTot = %e\n", pPosTot);
  // fprintf(stderr,"  StringProcess::append_finalArray pNegTot = %e\n", pNegTot);

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
  // fprintf(stderr,"  StringProcess::append_finalArray XVertexPos[nfrag] = %e\n",
  // XVertexPos[nfrag]);

  XVertexNeg[nfrag] = pNegTot / kappaString;
  for (kpstr = nfrag - 1; kpstr >= 0; kpstr--) {
    ipstr = indY[kpstr];

    XVertexNeg[kpstr] =
        XVertexNeg[kpstr + 1] -
        (Efrag[ipstr] - pparallel[ipstr]) / (kappaString * std::sqrt(2.));
  }
  // fprintf(stderr,"  StringProcess::append_finalArray XVertexNeg[0] = %e\n",
  // XVertexNeg[0]);

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
    //if (evecLong == NULL) {
    //  xfRS[1] = 0.;
    //  xfRS[2] = 0.;
    //  xfRS[3] = (XVertexPos[kpstr] - XVertexNeg[kpstr + 1]) / sqrt(2.);
    //} else {
    //  xfRS[1] =
    //      evecLong[1] * (XVertexPos[kpstr] - XVertexNeg[kpstr + 1]) / sqrt(2.);
    //  xfRS[2] =
    //      evecLong[2] * (XVertexPos[kpstr] - XVertexNeg[kpstr + 1]) / sqrt(2.);
    //  xfRS[3] =
    //      evecLong[3] * (XVertexPos[kpstr] - XVertexNeg[kpstr + 1]) / sqrt(2.);
    //}

    tProd = (XVertexPos[kpstr] + XVertexNeg[kpstr + 1]) / std::sqrt(2.);
    // zparallel = (XVertexPos[kpstr] - XVertexNeg[kpstr + 1])/sqrt(2.);
    // fprintf(stderr,"  StringProcess::append_finalArray tProd[kpstr = %d] = %e\n",
    // kpstr, tProd);

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

    // fprintf(stderr,"  StringProcess::append_finalArray adding to final arrays\n");

    final_PDGid[0].push_back(id);
    final_PDGid[1].push_back(islead);

    pvCM = pvRS.LorentzBoost( -vstring );
    E = pvCM.x0();
    px = pvCM.x1();
    py = pvCM.x2();
    pz = pvCM.x3();
    xfCM = xfRS.LorentzBoost( -vstring );
    tProd = xfCM.x0();
    //if (uString == NULL) {
    //  E = pvRS[0];
    //  px = pvRS[1];
    //  py = pvRS[2];
    //  pz = pvRS[3];

    //  tProd = xfRS[0];
    //  // fprintf(stderr,"  StringProcess::append_finalArray : NULL uString\n");
    //} else {
    //  lorentz->Boost1_RestToLab(3, uString, pvRS, pvCM);
    //  E = pvCM[0];
    //  px = pvCM[1];
    //  py = pvCM[2];
    //  pz = pvCM[3];

    //  lorentz->Boost1_RestToLab(3, uString, xfRS, xfCM);
    //  tProd = xfCM[0];
    //  // fprintf(stderr,"  StringProcess::append_finalArray : non-zero uString\n");
    //}

    final_pvec[0].push_back(E);
    final_pvec[1].push_back(px);
    final_pvec[2].push_back(py);
    final_pvec[3].push_back(pz);
    final_pvec[4].push_back(mass);

    final_tform[0].push_back(tProd);
    final_tform[1].push_back(xtotfac);

    ret = ret + 1;
  }

  //free(idfrag);
  //free(Efrag);
  //free(pxfrag);
  //free(pyfrag);
  //free(pzfrag);
  //free(mfrag);

  //free(indY);
  //free(pparallel);
  //free(Yparallel);
  //free(XVertexPos);
  //free(XVertexNeg);

  //free(pvRS);
  //free(pvCM);

  //free(xfRS);
  //free(xfCM);

  return ret;
}

// energy and charge conservation check
/*
bool StringProcess::check_conservation() {
  bool ret;

  int ipart, npart;

  CFIN = 0;
  BFIN = 0;
  EFIN = 0.;
  pxFIN = 0.;
  pyFIN = 0.;
  pzFIN = 0.;
  npart = final_PDGid[0].size();
  for (ipart = 0; ipart < npart; ipart++) {
    CFIN = CFIN + pythia->particleData.chargeType(final_PDGid[0].at(ipart));
    BFIN =
        BFIN + pythia->particleData.baryonNumberType(final_PDGid[0].at(ipart));
    EFIN = EFIN + final_pvec[0].at(ipart);
    pxFIN = pxFIN + final_pvec[1].at(ipart);
    pyFIN = pyFIN + final_pvec[2].at(ipart);
    pzFIN = pzFIN + final_pvec[3].at(ipart);
  }
  ret = true;

  if (CFIN != CINI) {
    fprintf(stderr,
            "  StringProcess::check_conservation : charge violation - CINI = %d, "
            "CFIN = %d\n",
            CINI, CFIN);
    ret = false;
  }

  if (BFIN != BINI) {
    fprintf(stderr,
            "  StringProcess::check_conservation : baryon number violation - BINI = "
            "%d, BFIN = %d\n",
            BINI, BFIN);
    ret = false;
  }

  if (fabs(EFIN - EINI) > 1.0e-8 * EINI) {
    fprintf(stderr,
            "  StringProcess::check_conservation : energy violation - EINI = %e GeV, "
            "EFIN = %e GeV\n",
            EINI, EFIN);
    ret = false;
  }

  if (fabs(pxFIN - pxINI) > 1.0e-8 * EINI) {
    fprintf(stderr,
            "  StringProcess::check_conservation : momentum x violation - pxINI = %e "
            "GeV, pxFIN = %e GeV\n",
            pxINI, pxFIN);
    ret = false;
  }

  if (fabs(pyFIN - pyINI) > 1.0e-8 * EINI) {
    fprintf(stderr,
            "  StringProcess::check_conservation : momentum y violation - pyINI = %e "
            "GeV, pyFIN = %e GeV\n",
            pyINI, pyFIN);
    ret = false;
  }

  if (fabs(pzFIN - pzINI) > 1.0e-8 * EINI) {
    fprintf(stderr,
            "  StringProcess::check_conservation : momentum z violation - pzINI = %e "
            "GeV, pzFIN = %e GeV\n",
            pzINI, pzFIN);
    ret = false;
  }

  return ret;
}
*/

bool StringProcess::init(const ParticleList &incomingList){
  bool ret;

  //int imu, inu;
  std::array<int, 3> qcontent;

  //double E, px, py, pz;
  //double ex, ey, et;
  //double theta, phi;

  PDGidA = incomingList[0].pdgcode().get_decimal();
  PDGidB = incomingList[1].pdgcode().get_decimal();
  massA = incomingList[0].effective_mass();
  massB = incomingList[1].effective_mass();

  plabA = incomingList[0].momentum();
  //plabA[0] = incomingList[0].momentum().x0();
  //plabA[1] = incomingList[0].momentum().x1();
  //plabA[2] = incomingList[0].momentum().x2();
  //plabA[3] = incomingList[0].momentum().x3();

  plabB = incomingList[1].momentum();
  //plabB[0] = incomingList[1].momentum().x0();
  //plabB[1] = incomingList[1].momentum().x1();
  //plabB[2] = incomingList[1].momentum().x2();
  //plabB[3] = incomingList[1].momentum().x3();

  //E = plabA[0] + plabB[0];
  //px = plabA[1] + plabB[1];
  //py = plabA[2] + plabB[2];
  //pz = plabA[3] + plabB[3];
  sqrtsAB = ( plabA + plabB ).abs();
  //sqrtsAB = sqrt(pow(fabs(E), 2.) - pow(fabs(px), 2.) - pow(fabs(py), 2.) -
  //               pow(fabs(pz), 2.));
  pabscomAB = pCM(sqrtsAB, massA, massB);
  //pabscomAB = sqrt((pow(fabs(sqrtsAB), 2.) - pow(fabs(massA + massB), 2.)) *
  //                 (pow(fabs(sqrtsAB), 2.) - pow(fabs(massA - massB), 2.))) /
  //            (2. * sqrtsAB);
  ucomAB = ( plabA + plabB )/sqrtsAB;
  //ucomAB[0] = E / sqrtsAB;
  //ucomAB[1] = px / sqrtsAB;
  //ucomAB[2] = py / sqrtsAB;
  //ucomAB[3] = pz / sqrtsAB;
  vcomAB = ucomAB.velocity();

  pcomA = plabA.LorentzBoost(vcomAB);
  //lorentz->Boost1_LabToRest(3, ucomAB, plabA, pcomA);
  pcomB = plabB.LorentzBoost(vcomAB);
  //lorentz->Boost1_LabToRest(3, ucomAB, plabB, pcomB);

  make_orthonormal_basis();
  /*
  if (fabs(pcomA.x3()) < (1. - 1.0e-8) * pabscomAB) {
    evecBasisAB[3] = pcomA.threevec() / pabscomAB;
    //for (imu = 1; imu <= 3; imu++) {
    //  evecBasisAB[3][imu] = pcomA[imu] / pabscomAB;
    //}
    theta = acos(evecBasisAB[3].x3());
    //theta = acos(evecBasisAB[3][3]);

    ex = evecBasisAB[3].x1();
    ey = evecBasisAB[3].x2();
    //ex = evecBasisAB[3][1];
    //ey = evecBasisAB[3][2];
    et = std::sqrt(ex*ex + ey*ey);
    if (ey > 0.) {
      phi = acos(ex / et);
    } else {
      phi = -acos(ex / et);
    }

    evecBasisAB[1].set_x1( cos(theta) * cos(phi) );
    evecBasisAB[1].set_x2( cos(theta) * sin(phi) );
    evecBasisAB[1].set_x3( -sin(theta) );
    //evecBasisAB[1][1] = cos(theta) * cos(phi);
    //evecBasisAB[1][2] = cos(theta) * sin(phi);
    //evecBasisAB[1][3] = -sin(theta);

    evecBasisAB[2].set_x1( -sin(phi) );
    evecBasisAB[2].set_x2( cos(phi) );
    evecBasisAB[2].set_x3( 0. );
    //evecBasisAB[2][1] = -sin(phi);
    //evecBasisAB[2][2] = cos(phi);
    //evecBasisAB[2][3] = 0.;
  } else {
    if (pcomA.x3() > 0.) {
      evecBasisAB[1] = ThreeVector(1., 0., 0.);
      //evecBasisAB[1].set_x1( 1. );
      //evecBasisAB[1].set_x2( 0. );
      //evecBasisAB[1].set_x3( 0. );

      evecBasisAB[2] = ThreeVector(0., 1., 0.);
      //evecBasisAB[2].set_x1( 0. );
      //evecBasisAB[2].set_x2( 1. );
      //evecBasisAB[2].set_x3( 0. );

      evecBasisAB[3] = ThreeVector(0., 0., 1.);
      //evecBasisAB[3].set_x1( 0. );
      //evecBasisAB[3].set_x2( 0. );
      //evecBasisAB[3].set_x3( 1. );
      //for (imu = 1; imu <= 3; imu++) {
      //  for (inu = 1; inu <= 3; inu++) {
      //    if (imu == inu) {
      //      evecBasisAB[imu][inu] = 1.;
      //    } else {
      //      evecBasisAB[imu][inu] = 0.;
      //    }
      //  }
      //}
    } else {
      evecBasisAB[1] = ThreeVector(0., 1., 0.);
      //evecBasisAB[1].set_x1( 0. );
      //evecBasisAB[1].set_x2( 1. );
      //evecBasisAB[1].set_x3( 0. );
      //evecBasisAB[1][1] = 0.;
      //evecBasisAB[1][2] = 1.;
      //evecBasisAB[1][3] = 0.;

      evecBasisAB[2] = ThreeVector(1., 0., 0.);
      //evecBasisAB[2].set_x1( 1. );
      //evecBasisAB[2].set_x2( 0. );
      //evecBasisAB[2].set_x3( 0. );
      //evecBasisAB[2][1] = 1.;
      //evecBasisAB[2][2] = 0.;
      //evecBasisAB[2][3] = 0.;

      evecBasisAB[3] = ThreeVector(0., 0., -1.);
      //evecBasisAB[3].set_x1( 0. );
      //evecBasisAB[3].set_x2( 0. );
      //evecBasisAB[3].set_x3( -1. );
      //evecBasisAB[3][1] = 0.;
      //evecBasisAB[3][2] = 0.;
      //evecBasisAB[3][3] = -1.;
    }
  }
  */
  compute_incoming_lightcone_momenta();

  // fprintf(stderr,"  StringProcess::init : evecBasisAB1 = (%e, %e, %e)\n",
  //	evecBasisAB[1][1], evecBasisAB[1][2], evecBasisAB[1][3]);
  // fprintf(stderr,"  StringProcess::init : evecBasisAB2 = (%e, %e, %e)\n",
  //	evecBasisAB[2][1], evecBasisAB[2][2], evecBasisAB[2][3]);
  // fprintf(stderr,"  StringProcess::init : evecBasisAB3 = (%e, %e, %e)\n",
  //	evecBasisAB[3][1], evecBasisAB[3][2], evecBasisAB[3][3]);

  qcontent = incomingList[0].pdgcode().quark_content();
  idqsetA[0] = static_cast<int>( incomingList[0].pdgcode().spin_degeneracy() );
  idqsetA[3] = qcontent[0];
  idqsetA[2] = qcontent[1];
  idqsetA[1] = qcontent[2];
  //PDGid2idqset(PDGidA, idqsetA);
  qcontent = incomingList[1].pdgcode().quark_content();
  idqsetB[0] = static_cast<int>( incomingList[1].pdgcode().spin_degeneracy() );
  idqsetB[3] = qcontent[0];
  idqsetB[2] = qcontent[1];
  idqsetB[1] = qcontent[2];
  //PDGid2idqset(PDGidB, idqsetB);

  xfracMin = pLightConeMin / sqrtsAB;

  // quantum numbers of hadron A
  baryonA = pythia->particleData.baryonNumberType(idqsetA[1]) +
            pythia->particleData.baryonNumberType(idqsetA[2]);
  //chargeA = pythia->particleData.chargeType(idqsetA[1]) +
  //          pythia->particleData.chargeType(idqsetA[2]);
  if (idqsetA[3] != 0) {
    baryonA = baryonA + pythia->particleData.baryonNumberType(idqsetA[3]);
    //chargeA = chargeA + pythia->particleData.chargeType(idqsetA[3]);
  }
  // quantum numbers of hadron B
  baryonB = pythia->particleData.baryonNumberType(idqsetB[1]) +
            pythia->particleData.baryonNumberType(idqsetB[2]);
  //chargeB = pythia->particleData.chargeType(idqsetB[1]) +
  //          pythia->particleData.chargeType(idqsetB[2]);
  if (idqsetB[3] != 0) {
    baryonB = baryonB + pythia->particleData.baryonNumberType(idqsetB[3]);
    //chargeB = chargeB + pythia->particleData.chargeType(idqsetB[3]);
  }

  /*
  BINI = baryonA + baryonB;
  CINI = chargeA + chargeB;
  EINI = ( plabA + plabB ).x0();
  //EINI = plabA[0] + plabB[0];
  pxINI = ( plabA + plabB ).x1();
  //pxINI = plabA[1] + plabB[1];
  pyINI = ( plabA + plabB ).x2();
  //pyINI = plabA[2] + plabB[2];
  pzINI = ( plabA + plabB ).x3();
  //pzINI = plabA[3] + plabB[3];
  */

  //ret = sigmaTot.calc(PDGidA, PDGidB, sqrtsAB);
  //XSecTot = sigmaTot.sigmaTot();
  //XSecEl = sigmaTot.sigmaEl();
  //XSecAX = sigmaTot.sigmaAX();
  //XSecXB = sigmaTot.sigmaXB();
  //XSecXX = sigmaTot.sigmaXX();
  //XSecAXB = sigmaTot.sigmaAXB();
  //XSecInel = XSecTot - XSecEl;
  //XSecND = XSecInel - XSecAX - XSecXB - XSecXX;

  //XSecSummed[0] = 0.;
  //XSecSummed[1] = XSecAX;
  //XSecSummed[2] = XSecSummed[1] + XSecXB;
  //XSecSummed[3] = XSecSummed[2] + XSecXX;
  //XSecSummed[4] = XSecSummed[3] + XSecND;

  ret = true;
  return ret;
}

bool StringProcess::init_lab(PdgCode &idAIn, PdgCode &idBIn, double massAIn, double massBIn,
                       Pythia8::Vec4 plabAIn, Pythia8::Vec4 plabBIn) {
  bool ret;
  //int imu, inu;
  std::array<int, 3> qcontent;

  //double E, px, py, pz;
  //double ex, ey, et;
  //double theta, phi;

  PDGidA = idAIn.get_decimal();
  PDGidB = idBIn.get_decimal();
  massA = massAIn;
  massB = massBIn;

  plabA.set_x0( plabAIn.e() );
  plabA.set_x1( plabAIn.px() );
  plabA.set_x2( plabAIn.py() );
  plabA.set_x3( plabAIn.pz() );
  //plabA[0] = plabAIn.e();
  //plabA[1] = plabAIn.px();
  //plabA[2] = plabAIn.py();
  //plabA[3] = plabAIn.pz();

  plabB = FourVector(plabBIn.e(), plabBIn.px(), plabBIn.py(), plabBIn.pz());
  plabB.set_x0( plabBIn.e() );
  plabB.set_x1( plabBIn.px() );
  plabB.set_x2( plabBIn.py() );
  plabB.set_x3( plabBIn.pz() );
  //plabB[0] = plabBIn.e();
  //plabB[1] = plabBIn.px();
  //plabB[2] = plabBIn.py();
  //plabB[3] = plabBIn.pz();

  //E = plabA[0] + plabB[0];
  //px = plabA[1] + plabB[1];
  //py = plabA[2] + plabB[2];
  //pz = plabA[3] + plabB[3];
  sqrtsAB = ( plabA + plabB ).abs();
  //sqrtsAB = sqrt(pow(fabs(E), 2.) - pow(fabs(px), 2.) - pow(fabs(py), 2.) -
  //               pow(fabs(pz), 2.));
  pabscomAB = pCM(sqrtsAB, massA, massB);
  //pabscomAB = sqrt((pow(fabs(sqrtsAB), 2.) - pow(fabs(massA + massB), 2.)) *
  //                 (pow(fabs(sqrtsAB), 2.) - pow(fabs(massA - massB), 2.))) /
  //            (2. * sqrtsAB);
  ucomAB = ( plabA + plabB )/sqrtsAB;
  //ucomAB[0] = E / sqrtsAB;
  //ucomAB[1] = px / sqrtsAB;
  //ucomAB[2] = py / sqrtsAB;
  //ucomAB[3] = pz / sqrtsAB;
  vcomAB = ucomAB.velocity();

  pcomA = plabA.LorentzBoost(vcomAB);
  //lorentz->Boost1_LabToRest(3, ucomAB, plabA, pcomA);
  pcomB = plabB.LorentzBoost(vcomAB);
  //lorentz->Boost1_LabToRest(3, ucomAB, plabB, pcomB);

  make_orthonormal_basis();
  /*
  if (fabs(pcomA.x3()) < (1. - 1.0e-8) * pabscomAB) {
    evecBasisAB[3] = pcomA.threevec() / pabscomAB;
    //for (imu = 1; imu <= 3; imu++) {
    //  evecBasisAB[3][imu] = pcomA[imu] / pabscomAB;
    //}
    theta = acos(evecBasisAB[3].x3());
    //theta = acos(evecBasisAB[3][3]);

    ex = evecBasisAB[3].x1();
    ey = evecBasisAB[3].x2();
    //ex = evecBasisAB[3][1];
    //ey = evecBasisAB[3][2];
    et = std::sqrt(ex*ex + ey*ey);
    if (ey > 0.) {
      phi = acos(ex / et);
    } else {
      phi = -acos(ex / et);
    }

    evecBasisAB[1].set_x1( cos(theta) * cos(phi) );
    evecBasisAB[1].set_x2( cos(theta) * sin(phi) );
    evecBasisAB[1].set_x3( -sin(theta) );
    //evecBasisAB[1][1] = cos(theta) * cos(phi);
    //evecBasisAB[1][2] = cos(theta) * sin(phi);
    //evecBasisAB[1][3] = -sin(theta);

    evecBasisAB[2].set_x1( -sin(phi) );
    evecBasisAB[2].set_x2( cos(phi) );
    evecBasisAB[2].set_x3( 0. );
    //evecBasisAB[2][1] = -sin(phi);
    //evecBasisAB[2][2] = cos(phi);
    //evecBasisAB[2][3] = 0.;
  } else {
    if (pcomA.x3() > 0.) {
      evecBasisAB[1] = ThreeVector(1., 0., 0.);
      //evecBasisAB[1].set_x1( 1. );
      //evecBasisAB[1].set_x2( 0. );
      //evecBasisAB[1].set_x3( 0. );

      evecBasisAB[2] = ThreeVector(0., 1., 0.);
      //evecBasisAB[2].set_x1( 0. );
      //evecBasisAB[2].set_x2( 1. );
      //evecBasisAB[2].set_x3( 0. );

      evecBasisAB[3] = ThreeVector(0., 0., 1.);
      //evecBasisAB[3].set_x1( 0. );
      //evecBasisAB[3].set_x2( 0. );
      //evecBasisAB[3].set_x3( 1. );
      //for (imu = 1; imu <= 3; imu++) {
      //  for (inu = 1; inu <= 3; inu++) {
      //    if (imu == inu) {
      //      evecBasisAB[imu][inu] = 1.;
      //    } else {
      //      evecBasisAB[imu][inu] = 0.;
      //    }
      //  }
      //}
    } else {
      evecBasisAB[1] = ThreeVector(0., 1., 0.);
      //evecBasisAB[1].set_x1( 0. );
      //evecBasisAB[1].set_x2( 1. );
      //evecBasisAB[1].set_x3( 0. );
      //evecBasisAB[1][1] = 0.;
      //evecBasisAB[1][2] = 1.;
      //evecBasisAB[1][3] = 0.;

      evecBasisAB[2] = ThreeVector(1., 0., 0.);
      //evecBasisAB[2].set_x1( 1. );
      //evecBasisAB[2].set_x2( 0. );
      //evecBasisAB[2].set_x3( 0. );
      //evecBasisAB[2][1] = 1.;
      //evecBasisAB[2][2] = 0.;
      //evecBasisAB[2][3] = 0.;

      evecBasisAB[3] = ThreeVector(0., 0., -1.);
      //evecBasisAB[3].set_x1( 0. );
      //evecBasisAB[3].set_x2( 0. );
      //evecBasisAB[3].set_x3( -1. );
      //evecBasisAB[3][1] = 0.;
      //evecBasisAB[3][2] = 0.;
      //evecBasisAB[3][3] = -1.;
    }
  }
  */
  compute_incoming_lightcone_momenta();

  // fprintf(stderr,"  StringProcess::init : evecBasisAB1 = (%e, %e, %e)\n",
  //	evecBasisAB[1][1], evecBasisAB[1][2], evecBasisAB[1][3]);
  // fprintf(stderr,"  StringProcess::init : evecBasisAB2 = (%e, %e, %e)\n",
  //	evecBasisAB[2][1], evecBasisAB[2][2], evecBasisAB[2][3]);
  // fprintf(stderr,"  StringProcess::init : evecBasisAB3 = (%e, %e, %e)\n",
  //	evecBasisAB[3][1], evecBasisAB[3][2], evecBasisAB[3][3]);

  qcontent = idAIn.quark_content();
  idqsetA[0] = static_cast<int>( idAIn.spin_degeneracy() );
  idqsetA[3] = qcontent[0];
  idqsetA[2] = qcontent[1];
  idqsetA[1] = qcontent[2];
  //PDGid2idqset(PDGidA, idqsetA);
  qcontent = idBIn.quark_content();
  idqsetB[0] = static_cast<int>( idBIn.spin_degeneracy() );
  idqsetB[3] = qcontent[0];
  idqsetB[2] = qcontent[1];
  idqsetB[1] = qcontent[2];
  //PDGid2idqset(PDGidB, idqsetB);

  xfracMin = pLightConeMin / sqrtsAB;

  // quantum numbers of hadron A
  baryonA = pythia->particleData.baryonNumberType(idqsetA[1]) +
            pythia->particleData.baryonNumberType(idqsetA[2]);
  //chargeA = pythia->particleData.chargeType(idqsetA[1]) +
  //          pythia->particleData.chargeType(idqsetA[2]);
  if (idqsetA[3] != 0) {
    baryonA = baryonA + pythia->particleData.baryonNumberType(idqsetA[3]);
    //chargeA = chargeA + pythia->particleData.chargeType(idqsetA[3]);
  }
  // quantum numbers of hadron B
  baryonB = pythia->particleData.baryonNumberType(idqsetB[1]) +
            pythia->particleData.baryonNumberType(idqsetB[2]);
  //chargeB = pythia->particleData.chargeType(idqsetB[1]) +
  //          pythia->particleData.chargeType(idqsetB[2]);
  if (idqsetB[3] != 0) {
    baryonB = baryonB + pythia->particleData.baryonNumberType(idqsetB[3]);
    //chargeB = chargeB + pythia->particleData.chargeType(idqsetB[3]);
  }

  /*
  BINI = baryonA + baryonB;
  CINI = chargeA + chargeB;
  EINI = ( plabA + plabB ).x0();
  //EINI = plabA[0] + plabB[0];
  pxINI = ( plabA + plabB ).x1();
  //pxINI = plabA[1] + plabB[1];
  pyINI = ( plabA + plabB ).x2();
  //pyINI = plabA[2] + plabB[2];
  pzINI = ( plabA + plabB ).x3();
  //pzINI = plabA[3] + plabB[3];
  */

  //ret = sigmaTot.calc(PDGidA, PDGidB, sqrtsAB);
  //XSecTot = sigmaTot.sigmaTot();
  //XSecEl = sigmaTot.sigmaEl();
  //XSecAX = sigmaTot.sigmaAX();
  //XSecXB = sigmaTot.sigmaXB();
  //XSecXX = sigmaTot.sigmaXX();
  //XSecAXB = sigmaTot.sigmaAXB();
  //XSecInel = XSecTot - XSecEl;
  //XSecND = XSecInel - XSecAX - XSecXB - XSecXX;

  //XSecSummed[0] = 0.;
  //XSecSummed[1] = XSecAX;
  //XSecSummed[2] = XSecSummed[1] + XSecXB;
  //XSecSummed[3] = XSecSummed[2] + XSecXX;
  //XSecSummed[4] = XSecSummed[3] + XSecND;

  ret = true;
  return ret;
}

bool StringProcess::init_com(PdgCode &idAIn, PdgCode &idBIn, double massAIn, double massBIn,
                       double sqrtsABIn) {
  bool ret;
  //int imu, inu;
  std::array<int, 3> qcontent;

  //double E, px, py, pz;

  PDGidA = idAIn.get_decimal();
  PDGidB = idBIn.get_decimal();
  massA = massAIn;
  massB = massBIn;
  sqrtsAB = sqrtsABIn;
  pabscomAB = pCM(sqrtsAB, massA, massB);
  //pabscomAB = sqrt((pow(fabs(sqrtsAB), 2.) - pow(fabs(massA + massB), 2.)) *
  //                 (pow(fabs(sqrtsAB), 2.) - pow(fabs(massA - massB), 2.))) /
  //            (2. * sqrtsAB);

  plabA = FourVector( std::sqrt(massA*massA + pabscomAB*pabscomAB), 0., 0., pabscomAB );
  //plabA.set_x0( sqrt(massA*massA + pabscomAB*pabscomAB) );
  //plabA.set_x1( 0. );
  //plabA.set_x2( 0. );
  //plabA.set_x3( pabscomAB );

  plabB = FourVector( std::sqrt(massB*massB + pabscomAB*pabscomAB), 0., 0., -pabscomAB );
  //plabB.set_x0( sqrt(massB*massB + pabscomAB*pabscomAB) );
  //plabB.set_x1( 0. );
  //plabB.set_x2( 0. );
  //plabB.set_x3( -pabscomAB );

  //E = plabA[0] + plabB[0];
  //px = plabA[1] + plabB[1];
  //py = plabA[2] + plabB[2];
  //pz = plabA[3] + plabB[3];
  ucomAB = ( plabA + plabB )/sqrtsAB;
  //ucomAB[0] = E / sqrtsAB;
  //ucomAB[1] = px / sqrtsAB;
  //ucomAB[2] = py / sqrtsAB;
  //ucomAB[3] = pz / sqrtsAB;
  vcomAB = ucomAB.velocity();

  pcomA = plabA.LorentzBoost(vcomAB);
  //lorentz->Boost1_LabToRest(3, ucomAB, plabA, pcomA);
  pcomB = plabB.LorentzBoost(vcomAB);
  //lorentz->Boost1_LabToRest(3, ucomAB, plabB, pcomB);

  make_orthonormal_basis();
  /*
  evecBasisAB[1] = ThreeVector(1., 0., 0.);
  //evecBasisAB[1].set_x1( 1. );
  //evecBasisAB[1].set_x2( 0. );
  //evecBasisAB[1].set_x3( 0. );

  evecBasisAB[2] = ThreeVector(0., 1., 0.);
  //evecBasisAB[2].set_x1( 0. );
  //evecBasisAB[2].set_x2( 1. );
  //evecBasisAB[2].set_x3( 0. );

  evecBasisAB[3] = ThreeVector(0., 0., 1.);
  //evecBasisAB[3].set_x1( 0. );
  //evecBasisAB[3].set_x2( 0. );
  //evecBasisAB[3].set_x3( 1. );
  //for (imu = 1; imu <= 3; imu++) {
  //  for (inu = 1; inu <= 3; inu++) {
  //    if (imu == inu) {
  //      evecBasisAB[imu][inu] = 1.;
  //    } else {
  //      evecBasisAB[imu][inu] = 0.;
  //    }
  //  }
  //}
  */
  compute_incoming_lightcone_momenta();

  // fprintf(stderr,"  StringProcess::init : evecBasisAB1 = (%e, %e, %e)\n",
  //	evecBasisAB[1][1], evecBasisAB[1][2], evecBasisAB[1][3]);
  // fprintf(stderr,"  StringProcess::init : evecBasisAB2 = (%e, %e, %e)\n",
  //	evecBasisAB[2][1], evecBasisAB[2][2], evecBasisAB[2][3]);
  // fprintf(stderr,"  StringProcess::init : evecBasisAB3 = (%e, %e, %e)\n",
  //	evecBasisAB[3][1], evecBasisAB[3][2], evecBasisAB[3][3]);

  qcontent = idAIn.quark_content();
  idqsetA[0] = static_cast<int>( idAIn.spin_degeneracy() );
  idqsetA[3] = qcontent[0];
  idqsetA[2] = qcontent[1];
  idqsetA[1] = qcontent[2];
  //PDGid2idqset(PDGidA, idqsetA);
  qcontent = idBIn.quark_content();
  idqsetB[0] = static_cast<int>( idBIn.spin_degeneracy() );
  idqsetB[3] = qcontent[0];
  idqsetB[2] = qcontent[1];
  idqsetB[1] = qcontent[2];
  //PDGid2idqset(PDGidB, idqsetB);

  xfracMin = pLightConeMin / sqrtsAB;

  // quantum numbers of hadron A
  baryonA = pythia->particleData.baryonNumberType(idqsetA[1]) +
            pythia->particleData.baryonNumberType(idqsetA[2]);
  //chargeA = pythia->particleData.chargeType(idqsetA[1]) +
  //          pythia->particleData.chargeType(idqsetA[2]);
  if (idqsetA[3] != 0) {
    baryonA = baryonA + pythia->particleData.baryonNumberType(idqsetA[3]);
    //chargeA = chargeA + pythia->particleData.chargeType(idqsetA[3]);
  }
  // quantum numbers of hadron B
  baryonB = pythia->particleData.baryonNumberType(idqsetB[1]) +
            pythia->particleData.baryonNumberType(idqsetB[2]);
  //chargeB = pythia->particleData.chargeType(idqsetB[1]) +
  //          pythia->particleData.chargeType(idqsetB[2]);
  if (idqsetB[3] != 0) {
    baryonB = baryonB + pythia->particleData.baryonNumberType(idqsetB[3]);
    //chargeB = chargeB + pythia->particleData.chargeType(idqsetB[3]);
  }

  /*
  BINI = baryonA + baryonB;
  CINI = chargeA + chargeB;
  EINI = ( plabA + plabB ).x0();
  //EINI = plabA[0] + plabB[0];
  pxINI = ( plabA + plabB ).x1();
  //pxINI = plabA[1] + plabB[1];
  pyINI = ( plabA + plabB ).x2();
  //pyINI = plabA[2] + plabB[2];
  pzINI = ( plabA + plabB ).x3();
  //pzINI = plabA[3] + plabB[3];
  */

  //ret = sigmaTot.calc(PDGidA, PDGidB, sqrtsAB);
  //XSecTot = sigmaTot.sigmaTot();
  //XSecEl = sigmaTot.sigmaEl();
  //XSecAX = sigmaTot.sigmaAX();
  //XSecXB = sigmaTot.sigmaXB();
  //XSecXX = sigmaTot.sigmaXX();
  //XSecAXB = sigmaTot.sigmaAXB();
  //XSecInel = XSecTot - XSecEl;
  //XSecND = XSecInel - XSecAX - XSecXB - XSecXX;

  //XSecSummed[0] = 0.;
  //XSecSummed[1] = XSecAX;
  //XSecSummed[2] = XSecSummed[1] + XSecXB;
  //XSecSummed[3] = XSecSummed[2] + XSecXX;
  //XSecSummed[4] = XSecSummed[3] + XSecND;

  ret = true;
  return ret;
}

/*
bool StringProcess::next_Inel() {
  bool ret;

  double rproc;
  rproc = XSecInel * Random::uniform(0., 1.);
  if ((rproc > XSecSummed[0]) && (rproc <= XSecSummed[1])) {
    ret = next_SDiff_AX();
  } else if ((rproc > XSecSummed[1]) && (rproc <= XSecSummed[2])) {
    ret = next_SDiff_XB();
  } else if ((rproc > XSecSummed[2]) && (rproc <= XSecSummed[3])) {
    ret = next_DDiff_XX();
  } else if ((rproc > XSecSummed[3]) && (rproc < XSecSummed[4])) {
    ret = next_NDiff();
  } else {
    ret = false;
  }

  return ret;
}
*/

// single diffractive AB > AX
bool StringProcess::next_SDiff_AX() {
  bool ret;
  //bool conserved;

  //int imu;
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
  //double *pstrHcom;
  //double *pstrHlab;
  //double *pstrXcom;
  //double *pstrXlab;
  ThreeVector threeMomentum;

  FourVector ustrHcom;
  FourVector ustrHlab;
  FourVector ustrXcom;
  FourVector ustrXlab;
  //double *ustrHcom;
  //double *ustrHlab;
  //double *ustrXcom;
  //double *ustrXlab;

  double pabs;
  FourVector pnull;
  //double *pnull;
  FourVector prs;
  //double *prs;
  ThreeVector evec;
  //double *evec;

  //pstrHcom = static_cast<double *>(malloc(4 * sizeof(double)));
  //pstrHlab = static_cast<double *>(malloc(4 * sizeof(double)));
  //pstrXcom = static_cast<double *>(malloc(4 * sizeof(double)));
  //pstrXlab = static_cast<double *>(malloc(4 * sizeof(double)));

  //ustrHcom = static_cast<double *>(malloc(4 * sizeof(double)));
  //ustrHlab = static_cast<double *>(malloc(4 * sizeof(double)));
  //ustrXcom = static_cast<double *>(malloc(4 * sizeof(double)));
  //ustrXlab = static_cast<double *>(malloc(4 * sizeof(double)));

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

    makeStringEnds(idqsetB, &idqX1, &idqX2);

    QTrx = Random::normal(0., sigmaQperp/std::sqrt(2.) );
    QTry = Random::normal(0., sigmaQperp/std::sqrt(2.) );
    QTrn = std::sqrt(QTrx*QTrx + QTry*QTry);
    //QTrn = sample_Qperp(sigmaQperp);
    //phiQ = 2. * M_PI * Random::uniform(0., 1.);
    // fprintf(stderr,"  StringProcess::next_SDiff_AX : QTrn = %e GeV, phiQ = %e\n",
    // QTrn, phiQ);

    rmass = std::log(mstrMax / mstrMin) * Random::uniform(0., 1.);
    massX = mstrMin * exp(rmass);
    pabscomAX = pCM( sqrtsAB, massA, massX );
    //pabscomAX = sqrt((pow(fabs(sqrtsAB), 2.) - pow(fabs(massA + massX), 2.)) *
    //                 (pow(fabs(sqrtsAB), 2.) - pow(fabs(massA - massX), 2.))) /
    //            (2. * sqrtsAB);

    foundPabsX = pabscomAX > QTrn;
    foundMassX = massX > (pythia->particleData.m0(idqX1) +
                          pythia->particleData.m0(idqX2));
  }

  ret = false;
  if ((foundPabsX == true) && (foundMassX == true)) {
    //pstrHcom[0] = sqrt(pow(fabs(pabscomAX), 2.) + pow(fabs(massA), 2.));
    //pstrXcom[0] = sqrt(pow(fabs(pabscomAX), 2.) + pow(fabs(massX), 2.));
    threeMomentum = evecBasisAB[3] * std::sqrt(pabscomAX*pabscomAX - QTrn*QTrn) +
                        evecBasisAB[1] * QTrx +
                        evecBasisAB[2] * QTry;
    pstrHcom = FourVector( std::sqrt(pabscomAX*pabscomAX + massA*massA), threeMomentum );
    //pstrHcom.set_x0( std::sqrt(pabscomAX*pabscomAX + massA*massA) );
    //pstrHcom.set_x1( threeMomentum.x1() );
    //pstrHcom.set_x2( threeMomentum.x2() );
    //pstrHcom.set_x3( threeMomentum.x3() );
    threeMomentum = -evecBasisAB[3] * std::sqrt(pabscomAX*pabscomAX - QTrn*QTrn) -
                        evecBasisAB[1] * QTrx -
                        evecBasisAB[2] * QTry;
    pstrXcom = FourVector( std::sqrt(pabscomAX*pabscomAX + massX*massX), threeMomentum );
    //pstrXcom.set_x0( std::sqrt(pabscomAX*pabscomAX + massX*massX) );
    //pstrXcom.set_x1( threeMomentum.x1() );
    //pstrXcom.set_x2( threeMomentum.x2() );
    //pstrXcom.set_x3( threeMomentum.x3() );
    //for (imu = 1; imu < 4; imu++) {
    //  pstrHcom[imu] = evecBasisAB[3][imu] *
    //                      sqrt(pow(fabs(pabscomAX), 2.) - pow(fabs(QTrn), 2.)) +
    //                  evecBasisAB[1][imu] * QTrn * cos(phiQ) +
    //                  evecBasisAB[2][imu] * QTrn * sin(phiQ);
    //  pstrXcom[imu] = -evecBasisAB[3][imu] *
    //                      sqrt(pow(fabs(pabscomAX), 2.) - pow(fabs(QTrn), 2.)) -
    //                  evecBasisAB[1][imu] * QTrn * cos(phiQ) -
    //                  evecBasisAB[2][imu] * QTrn * sin(phiQ);
    //}
    // fprintf(stderr,"  StringProcess::next_SDiff_AX : pstrXcom = (%e, %e, %e, %e)
    // GeV\n",
    //	pstrXcom[0], pstrXcom[1], pstrXcom[2], pstrXcom[3]);

    // fprintf(stderr,"  StringProcess::next_SDiff_AX : idqX1 = %d, idqX2 = %d\n",
    // idqX1, idqX2);
    // fprintf(stderr,"  StringProcess::next_SDiff_AX : massX = %e GeV\n", massX);

    //pnull = static_cast<double *>(malloc(4 * sizeof(double)));
    //prs = static_cast<double *>(malloc(4 * sizeof(double)));
    //evec = static_cast<double *>(malloc(4 * sizeof(double)));

    pstrHlab = pstrHcom.LorentzBoost( -vcomAB );
    pstrXlab = pstrXcom.LorentzBoost( -vcomAB );
    //lorentz->Boost1_RestToLab(3, ucomAB, pstrHcom, pstrHlab);
    //lorentz->Boost1_RestToLab(3, ucomAB, pstrXcom, pstrXlab);
    ustrHcom = pstrHcom / massA;
    ustrXcom = pstrXcom / massX;
    ustrHlab = pstrHlab / massA;
    ustrXlab = pstrXlab / massX;
    //for (imu = 0; imu < 4; imu++) {
    //  ustrHcom[imu] = pstrHcom[imu] / massA;
    //  ustrXcom[imu] = pstrXcom[imu] / massX;

    //  ustrHlab[imu] = pstrHlab[imu] / massA;
    //  ustrXlab[imu] = pstrXlab[imu] / massX;
    //}

    threeMomentum = pstrXcom.threevec();
    pnull = FourVector( threeMomentum.abs(), threeMomentum );
    //pnull.set_x0( threeMomentum.abs() );
    //pnull.set_x1( threeMomentum.x1() );
    //pnull.set_x2( threeMomentum.x2() );
    //pnull.set_x3( threeMomentum.x3() );
    //for (imu = 1; imu < 4; imu++) {
    //  pnull[imu] = pstrXcom[imu];
    //}
    //pnull[0] = sqrt(pow(fabs(pnull[1]), 2.) + pow(fabs(pnull[2]), 2.) +
    //                pow(fabs(pnull[3]), 2.));
    prs = pnull.LorentzBoost( ustrXcom.velocity() );
    pabs = prs.threevec().abs();
    //lorentz->Boost1_LabToRest(3, ustrXcom, pnull, prs);
    //pabs = sqrt(pow(fabs(prs[1]), 2.) + pow(fabs(prs[2]), 2.) +
    //            pow(fabs(prs[3]), 2.));
    evec = prs.threevec() / pabs;
    //evec[0] = 0.;
    //for (imu = 1; imu < 4; imu++) {
    //  evec[imu] = prs[imu] / pabs;
    //}
    // fprintf(stderr,"  StringProcess::next_SDiff_AX : evec = (%e, %e, %e)\n",
    // evec[1], evec[2], evec[3]);
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

    //free(pnull);
    //free(prs);
    //free(evec);
  }

  //free(pstrHcom);
  //free(pstrHlab);
  //free(pstrXcom);
  //free(pstrXlab);

  //free(ustrHcom);
  //free(ustrHlab);
  //free(ustrXcom);
  //free(ustrXlab);

  //if (ret == true) {
  //  conserved = check_conservation();
  //  ret = ret && conserved;
  //}
  return ret;
}

// single diffractive AB > XB
bool StringProcess::next_SDiff_XB() {
  bool ret;
  //bool conserved;

  //int imu;
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
  //double *pstrHcom;
  //double *pstrHlab;
  //double *pstrXcom;
  //double *pstrXlab;
  ThreeVector threeMomentum;

  FourVector ustrHcom;
  FourVector ustrHlab;
  FourVector ustrXcom;
  FourVector ustrXlab;
  //double *ustrHcom;
  //double *ustrHlab;
  //double *ustrXcom;
  //double *ustrXlab;

  double pabs;
  FourVector pnull;
  //double *pnull;
  FourVector prs;
  //double *prs;
  ThreeVector evec;
  //double *evec;

  //pstrHcom = static_cast<double *>(malloc(4 * sizeof(double)));
  //pstrHlab = static_cast<double *>(malloc(4 * sizeof(double)));
  //pstrXcom = static_cast<double *>(malloc(4 * sizeof(double)));
  //pstrXlab = static_cast<double *>(malloc(4 * sizeof(double)));

  //ustrHcom = static_cast<double *>(malloc(4 * sizeof(double)));
  //ustrHlab = static_cast<double *>(malloc(4 * sizeof(double)));
  //ustrXcom = static_cast<double *>(malloc(4 * sizeof(double)));
  //ustrXlab = static_cast<double *>(malloc(4 * sizeof(double)));

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

    makeStringEnds(idqsetA, &idqX1, &idqX2);

    QTrx = Random::normal(0., sigmaQperp/std::sqrt(2.) );
    QTry = Random::normal(0., sigmaQperp/std::sqrt(2.) );
    QTrn = std::sqrt(QTrx*QTrx + QTry*QTry);
    //QTrn = sample_Qperp(sigmaQperp);
    //phiQ = 2. * M_PI * Random::uniform(0., 1.);
    // fprintf(stderr,"  StringProcess::next_SDiff_XB : QTrn = %e GeV, phiQ = %e\n",
    // QTrn, phiQ);

    rmass = log(mstrMax / mstrMin) * Random::uniform(0., 1.);
    massX = mstrMin * exp(rmass);
    pabscomXB = pCM( sqrtsAB, massX, massB );
    //pabscomXB = sqrt((pow(fabs(sqrtsAB), 2.) - pow(fabs(massB + massX), 2.)) *
    //                 (pow(fabs(sqrtsAB), 2.) - pow(fabs(massB - massX), 2.))) /
    //            (2. * sqrtsAB);

    foundPabsX = pabscomXB > QTrn;
    foundMassX = massX > (pythia->particleData.m0(idqX1) +
                          pythia->particleData.m0(idqX2));
  }

  ret = false;
  if ((foundPabsX == true) && (foundMassX == true)) {
    //pstrHcom[0] = sqrt(pow(fabs(pabscomXB), 2.) + pow(fabs(massB), 2.));
    //pstrXcom[0] = sqrt(pow(fabs(pabscomXB), 2.) + pow(fabs(massX), 2.));
    threeMomentum = -evecBasisAB[3] * std::sqrt(pabscomXB*pabscomXB - QTrn*QTrn) -
                        evecBasisAB[1] * QTrx -
                        evecBasisAB[2] * QTry;
    pstrHcom = FourVector( std::sqrt(pabscomXB*pabscomXB + massB*massB), threeMomentum );
    //pstrHcom.set_x0( std::sqrt(pabscomXB*pabscomXB + massB*massB) );
    //pstrHcom.set_x1( threeMomentum.x1() );
    //pstrHcom.set_x2( threeMomentum.x2() );
    //pstrHcom.set_x3( threeMomentum.x3() );
    threeMomentum = evecBasisAB[3] * std::sqrt(pabscomXB*pabscomXB - QTrn*QTrn) +
                        evecBasisAB[1] * QTrx +
                        evecBasisAB[2] * QTry;
    pstrXcom = FourVector( std::sqrt(pabscomXB*pabscomXB + massX*massX), threeMomentum );
    //pstrXcom.set_x0( std::sqrt(pabscomXB*pabscomXB + massX*massX) );
    //pstrXcom.set_x1( threeMomentum.x1() );
    //pstrXcom.set_x2( threeMomentum.x2() );
    //pstrXcom.set_x3( threeMomentum.x3() );
    //for (imu = 1; imu < 4; imu++) {
    //  pstrHcom[imu] = -evecBasisAB[3][imu] *
    //                      sqrt(pow(fabs(pabscomXB), 2.) - pow(fabs(QTrn), 2.)) -
    //                  evecBasisAB[1][imu] * QTrn * cos(phiQ) -
    //                  evecBasisAB[2][imu] * QTrn * sin(phiQ);
    //  pstrXcom[imu] = evecBasisAB[3][imu] *
    //                      sqrt(pow(fabs(pabscomXB), 2.) - pow(fabs(QTrn), 2.)) +
    //                  evecBasisAB[1][imu] * QTrn * cos(phiQ) +
    //                  evecBasisAB[2][imu] * QTrn * sin(phiQ);
    //}
    // fprintf(stderr,"  StringProcess::next_SDiff_XB : pstrXcom = (%e, %e, %e, %e)
    // GeV\n",
    //	pstrXcom[0], pstrXcom[1], pstrXcom[2], pstrXcom[3]);

    // fprintf(stderr,"  StringProcess::next_SDiff_XB : idqX1 = %d, idqX2 = %d\n",
    // idqX1, idqX2);
    // fprintf(stderr,"  StringProcess::next_SDiff_XB : massX = %e GeV\n", massX);

    //pnull = static_cast<double *>(malloc(4 * sizeof(double)));
    //prs = static_cast<double *>(malloc(4 * sizeof(double)));
    //evec = static_cast<double *>(malloc(4 * sizeof(double)));

    pstrHlab = pstrHcom.LorentzBoost( -vcomAB );
    pstrXlab = pstrXcom.LorentzBoost( -vcomAB );
    //lorentz->Boost1_RestToLab(3, ucomAB, pstrHcom, pstrHlab);
    //lorentz->Boost1_RestToLab(3, ucomAB, pstrXcom, pstrXlab);
    ustrHcom = pstrHcom / massB;
    ustrXcom = pstrXcom / massX;
    ustrHlab = pstrHlab / massB;
    ustrXlab = pstrXlab / massX;
    //for (imu = 0; imu < 4; imu++) {
    //  ustrHcom[imu] = pstrHcom[imu] / massB;
    //  ustrXcom[imu] = pstrXcom[imu] / massX;

    //  ustrHlab[imu] = pstrHlab[imu] / massB;
    //  ustrXlab[imu] = pstrXlab[imu] / massX;
    //}

    threeMomentum = pstrXcom.threevec();
    pnull = FourVector( threeMomentum.abs(), threeMomentum );
    //pnull.set_x0( threeMomentum.abs() );
    //pnull.set_x1( threeMomentum.x1() );
    //pnull.set_x2( threeMomentum.x2() );
    //pnull.set_x3( threeMomentum.x3() );
    //for (imu = 1; imu < 4; imu++) {
    //  pnull[imu] = pstrXcom[imu];
    //}
    //pnull[0] = sqrt(pow(fabs(pnull[1]), 2.) + pow(fabs(pnull[2]), 2.) +
    //                pow(fabs(pnull[3]), 2.));
    prs = pnull.LorentzBoost( ustrXcom.velocity() );
    pabs = prs.threevec().abs();
    //lorentz->Boost1_LabToRest(3, ustrXcom, pnull, prs);
    //pabs = sqrt(pow(fabs(prs[1]), 2.) + pow(fabs(prs[2]), 2.) +
    //            pow(fabs(prs[3]), 2.));
    evec = prs.threevec() / pabs;
    //evec[0] = 0.;
    //for (imu = 1; imu < 4; imu++) {
    //  evec[imu] = prs[imu] / pabs;
    //}
    // fprintf(stderr,"  StringProcess::next_SDiff_XB : evec = (%e, %e, %e)\n",
    // evec[1], evec[2], evec[3]);
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

    //free(pnull);
    //free(prs);
    //free(evec);
  }

  //free(pstrHcom);
  //free(pstrHlab);
  //free(pstrXcom);
  //free(pstrXlab);

  //free(ustrHcom);
  //free(ustrHlab);
  //free(ustrXcom);
  //free(ustrXlab);

  //if (ret == true) {
  //  conserved = check_conservation();
  //  ret = ret && conserved;
  //}
  return ret;
}

// double diffractive AB > XX
bool StringProcess::next_DDiff_XX() {
  bool ret;
  //bool conserved;

  //int imu;
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
  //double *pstr1com;
  //double *pstr1lab;
  //double *pstr2com;
  //double *pstr2lab;
  ThreeVector threeMomentum;

  FourVector ustr1com;
  FourVector ustr1lab;
  FourVector ustr2com;
  FourVector ustr2lab;
  //double *ustr1com;
  //double *ustr1lab;
  //double *ustr2com;
  //double *ustr2lab;

  double pabs;
  FourVector pnull;
  //double *pnull;
  FourVector prs;
  //double *prs;
  ThreeVector evec;
  //double *evec;

  //pstr1com = static_cast<double *>(malloc(4 * sizeof(double)));
  //pstr1lab = static_cast<double *>(malloc(4 * sizeof(double)));
  //pstr2com = static_cast<double *>(malloc(4 * sizeof(double)));
  //pstr2lab = static_cast<double *>(malloc(4 * sizeof(double)));

  //ustr1com = static_cast<double *>(malloc(4 * sizeof(double)));
  //ustr1lab = static_cast<double *>(malloc(4 * sizeof(double)));
  //ustr2com = static_cast<double *>(malloc(4 * sizeof(double)));
  //ustr2lab = static_cast<double *>(malloc(4 * sizeof(double)));

  reset_finalArray();

  // PDGid2idqset(PDGidA, idqsetA);
  // PDGid2idqset(PDGidB, idqsetB);
  //compute_incoming_lightcone_momenta();
  //PPosA = (pcomA.x0() + evecBasisAB[3].x1() * pcomA.x1() +
  //         evecBasisAB[3].x2() * pcomA.x2() + evecBasisAB[3].x3() * pcomA.x3()) /
  //        std::sqrt(2.);
  //PNegA = (pcomA.x0() - evecBasisAB[3].x1() * pcomA.x1() -
  //         evecBasisAB[3].x2() * pcomA.x2() - evecBasisAB[3].x3() * pcomA.x3()) /
  //        std::sqrt(2.);
  //PPosB = (pcomB.x0() + evecBasisAB[3].x1() * pcomB.x1() +
  //         evecBasisAB[3].x2() * pcomB.x2() + evecBasisAB[3].x3() * pcomB.x3()) /
  //        std::sqrt(2.);
  //PNegB = (pcomB.x0() - evecBasisAB[3].x1() * pcomB.x1() -
  //         evecBasisAB[3].x2() * pcomB.x2() - evecBasisAB[3].x3() * pcomB.x3()) /
  //        std::sqrt(2.);

  ntry = 0;
  foundMass1 = false;
  foundMass2 = false;
  while (((foundMass1 == false) || (foundMass2 == false)) && (ntry < 100)) {
    ntry = ntry + 1;

    makeStringEnds(idqsetA, &idq11, &idq12);
    makeStringEnds(idqsetB, &idq21, &idq22);

    xfracA = Random::beta_a0(xfracMin, betapowS + 1.);
    xfracB = Random::beta_a0(xfracMin, betapowS + 1.);
    //xfracA = sample_XSDIS(xfracMin, betapowS);
    //xfracB = sample_XSDIS(xfracMin, betapowS);
    QTrx = Random::normal(0., sigmaQperp/std::sqrt(2.) );
    QTry = Random::normal(0., sigmaQperp/std::sqrt(2.) );
    QTrn = std::sqrt(QTrx*QTrx + QTry*QTry);
    //QTrn = sample_Qperp(sigmaQperp);
    //phiQ = 2. * M_PI * Random::uniform(0., 1.);
    // fprintf(stderr,"  StringProcess::next_DDiff : xfracA = %e, xfracB = %e\n",
    // xfracA, xfracB);
    // fprintf(stderr,"  StringProcess::next_DDiff : QTrn = %e GeV, phiQ = %e\n",
    // QTrn, phiQ);

    QPos = -pow(fabs(QTrn), 2.) / (2. * xfracB * PNegB);
    QNeg = pow(fabs(QTrn), 2.) / (2. * xfracA * PPosA);

    //pstr1com[0] = (PPosA + QPos + PNegA + QNeg) / sqrt(2.);
    //pstr2com[0] = (PPosB - QPos + PNegB - QNeg) / sqrt(2.);
    //mstr1 = pow(fabs(pstr1com[0]), 2.);
    //mstr2 = pow(fabs(pstr2com[0]), 2.);
    threeMomentum = evecBasisAB[3] * (PPosA + QPos - PNegA - QNeg) / std::sqrt(2.) +
                        evecBasisAB[1] * QTrx + evecBasisAB[2] * QTry;
    pstr1com = FourVector( (PPosA + QPos + PNegA + QNeg) / std::sqrt(2.), threeMomentum );
    //pstr1com.set_x0( (PPosA + QPos + PNegA + QNeg) / std::sqrt(2.) );
    //pstr1com.set_x1( threeMomentum.x1() );
    //pstr1com.set_x2( threeMomentum.x2() );
    //pstr1com.set_x3( threeMomentum.x3() );
    mstr1 = pstr1com.sqr();
    threeMomentum = evecBasisAB[3] * (PPosB - QPos - PNegB + QNeg) / std::sqrt(2.) -
                        evecBasisAB[1] * QTrx - evecBasisAB[2] * QTry;
    pstr2com = FourVector( (PPosB - QPos + PNegB - QNeg) / std::sqrt(2.), threeMomentum );
    //pstr2com.set_x0( (PPosB - QPos + PNegB - QNeg) / std::sqrt(2.) );
    //pstr2com.set_x1( threeMomentum.x1() );
    //pstr2com.set_x2( threeMomentum.x2() );
    //pstr2com.set_x3( threeMomentum.x3() );
    mstr2 = pstr2com.sqr();
    //for (imu = 1; imu < 4; imu++) {
    //  pstr1com[imu] =
    //      evecBasisAB[3][imu] * (PPosA + QPos - PNegA - QNeg) / sqrt(2.) +
    //      evecBasisAB[1][imu] * QTrn * cos(phiQ) +
    //      evecBasisAB[2][imu] * QTrn * sin(phiQ);
    //  pstr2com[imu] =
    //      evecBasisAB[3][imu] * (PPosB - QPos - PNegB + QNeg) / sqrt(2.) -
    //      evecBasisAB[1][imu] * QTrn * cos(phiQ) -
    //      evecBasisAB[2][imu] * QTrn * sin(phiQ);
    //  mstr1 = mstr1 - pow(fabs(pstr1com[imu]), 2.);
    //  mstr2 = mstr2 - pow(fabs(pstr2com[imu]), 2.);
    //}
    // fprintf(stderr,"  StringProcess::next_DDiff : pstr1com = (%e, %e, %e, %e)
    // GeV\n",
    //	pstr1com[0], pstr1com[1], pstr1com[2], pstr1com[3]);
    // fprintf(stderr,"  StringProcess::next_DDiff : pstr2com = (%e, %e, %e, %e)
    // GeV\n",
    //	pstr2com[0], pstr2com[1], pstr2com[2], pstr2com[3]);

    mstr1 = (mstr1 > 0.) ? std::sqrt(mstr1) : 0.;
    //if (mstr1 > 0.) {
    //  mstr1 = std::sqrt(mstr1);
    //} else {
    //  mstr1 = 0.;
    //}
    foundMass1 = mstr1 > (pythia->particleData.m0(idq11) +
                          pythia->particleData.m0(idq12));

    mstr2 = (mstr2 > 0.) ? std::sqrt(mstr2) : 0.;
    //if (mstr2 > 0.) {
    //  mstr2 = std::sqrt(mstr2);
    //} else {
    //  mstr2 = 0.;
    //}
    foundMass2 = mstr2 > (pythia->particleData.m0(idq21) +
                          pythia->particleData.m0(idq22));
  }

  ret = false;
  bool both_masses_above_pythia_threshold = foundMass1 && foundMass2;
  if ( both_masses_above_pythia_threshold ) {
    // fprintf(stderr,"  StringProcess::next_DDiff : idq11 = %d, idq12 = %d\n", idq11,
    // idq12);
    // fprintf(stderr,"  StringProcess::next_DDiff : mstr1 = %e GeV\n", mstr1);

    // fprintf(stderr,"  StringProcess::next_DDiff : idq21 = %d, idq22 = %d\n", idq21,
    // idq22);
    // fprintf(stderr,"  StringProcess::next_DDiff : mstr2 = %e GeV\n", mstr2);

    //pnull = static_cast<double *>(malloc(4 * sizeof(double)));
    //prs = static_cast<double *>(malloc(4 * sizeof(double)));
    //evec = static_cast<double *>(malloc(4 * sizeof(double)));

    pstr1lab = pstr1com.LorentzBoost( -vcomAB );
    pstr2lab = pstr2com.LorentzBoost( -vcomAB );
    //lorentz->Boost1_RestToLab(3, ucomAB, pstr1com, pstr1lab);
    //lorentz->Boost1_RestToLab(3, ucomAB, pstr2com, pstr2lab);
    ustr1com = pstr1com / mstr1;
    ustr2com = pstr2com / mstr2;
    ustr1lab = pstr1lab / mstr1;
    ustr2lab = pstr2lab / mstr2;
    //for (imu = 0; imu < 4; imu++) {
    //  ustr1com[imu] = pstr1com[imu] / mstr1;
    //  ustr2com[imu] = pstr2com[imu] / mstr2;

    //  ustr1lab[imu] = pstr1lab[imu] / mstr1;
    //  ustr2lab[imu] = pstr2lab[imu] / mstr2;
    //}

    threeMomentum = pstr1com.threevec();
    pnull = FourVector( threeMomentum.abs(), threeMomentum );
    //pnull.set_x0( threeMomentum.abs() );
    //pnull.set_x1( threeMomentum.x1() );
    //pnull.set_x2( threeMomentum.x2() );
    //pnull.set_x3( threeMomentum.x3() );
    //for (imu = 1; imu < 4; imu++) {
    //  pnull[imu] = pstr1com[imu];
    //}
    //pnull[0] = sqrt(pow(fabs(pnull[1]), 2.) + pow(fabs(pnull[2]), 2.) +
    //                pow(fabs(pnull[3]), 2.));
    prs = pnull.LorentzBoost( ustr1com.velocity() );
    pabs = prs.threevec().abs();
    //lorentz->Boost1_LabToRest(3, ustr1com, pnull, prs);
    //pabs = sqrt(pow(fabs(prs[1]), 2.) + pow(fabs(prs[2]), 2.) +
    //            pow(fabs(prs[3]), 2.));
    evec = prs.threevec() / pabs;
    //evec[0] = 0.;
    //for (imu = 1; imu < 4; imu++) {
    //  evec[imu] = prs[imu] / pabs;
    //}
    // fprintf(stderr,"  StringProcess::next_DDiff : evec = (%e, %e, %e)\n", evec[1],
    // evec[2], evec[3]);
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
    //pnull.set_x0( threeMomentum.abs() );
    //pnull.set_x1( threeMomentum.x1() );
    //pnull.set_x2( threeMomentum.x2() );
    //pnull.set_x3( threeMomentum.x3() );
    //for (imu = 1; imu < 4; imu++) {
    //  pnull[imu] = pstr2com[imu];
    //}
    //pnull[0] = sqrt(pow(fabs(pnull[1]), 2.) + pow(fabs(pnull[2]), 2.) +
    //                pow(fabs(pnull[3]), 2.));
    prs = pnull.LorentzBoost( ustr2com.velocity() );
    pabs = prs.threevec().abs();
    //lorentz->Boost1_LabToRest(3, ustr2com, pnull, prs);
    //pabs = sqrt(pow(fabs(prs[1]), 2.) + pow(fabs(prs[2]), 2.) +
    //            pow(fabs(prs[3]), 2.));
    evec = prs.threevec() / pabs;
    //evec[0] = 0.;
    //for (imu = 1; imu < 4; imu++) {
    //  evec[imu] = prs[imu] / pabs;
    //}
    // fprintf(stderr,"  StringProcess::next_DDiff : evec = (%e, %e, %e)\n", evec[1],
    // evec[2], evec[3]);
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

    //free(pnull);
    //free(prs);
    //free(evec);
  }

  //free(pstr1com);
  //free(pstr1lab);
  //free(pstr2com);
  //free(pstr2lab);

  //free(ustr1com);
  //free(ustr1lab);
  //free(ustr2com);
  //free(ustr2lab);

  //if (ret == true) {
  //  conserved = check_conservation();
  //  ret = ret && conserved;
  //}
  return ret;
}

// non-diffractive
bool StringProcess::next_NDiff() {
  bool ret;
  //bool conserved;

  //int imu;
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
  //double *pstr1com;
  //double *pstr1lab;
  //double *pstr2com;
  //double *pstr2lab;
  ThreeVector threeMomentum;

  FourVector ustr1com;
  FourVector ustr1lab;
  FourVector ustr2com;
  FourVector ustr2lab;
  //double *ustr1com;
  //double *ustr1lab;
  //double *ustr2com;
  //double *ustr2lab;

  double pabs;
  FourVector pnull;
  //double *pnull;
  FourVector prs;
  //double *prs;
  ThreeVector evec;
  //double *evec;

  //pstr1com = static_cast<double *>(malloc(4 * sizeof(double)));
  //pstr1lab = static_cast<double *>(malloc(4 * sizeof(double)));
  //pstr2com = static_cast<double *>(malloc(4 * sizeof(double)));
  //pstr2lab = static_cast<double *>(malloc(4 * sizeof(double)));

  //ustr1com = static_cast<double *>(malloc(4 * sizeof(double)));
  //ustr1lab = static_cast<double *>(malloc(4 * sizeof(double)));
  //ustr2com = static_cast<double *>(malloc(4 * sizeof(double)));
  //ustr2lab = static_cast<double *>(malloc(4 * sizeof(double)));

  reset_finalArray();

  // PDGid2idqset(PDGidA, idqsetA);
  // PDGid2idqset(PDGidB, idqsetB);
  //compute_incoming_lightcone_momenta();
  //PPosA = (pcomA.x0() + evecBasisAB[3].x1() * pcomA.x1() +
  //         evecBasisAB[3].x2() * pcomA.x2() + evecBasisAB[3].x3() * pcomA.x3()) /
  //        std::sqrt(2.);
  //PNegA = (pcomA.x0() - evecBasisAB[3].x1() * pcomA.x1() -
  //         evecBasisAB[3].x2() * pcomA.x2() - evecBasisAB[3].x3() * pcomA.x3()) /
  //        std::sqrt(2.);
  //PPosB = (pcomB.x0() + evecBasisAB[3].x1() * pcomB.x1() +
  //         evecBasisAB[3].x2() * pcomB.x2() + evecBasisAB[3].x3() * pcomB.x3()) /
  //        std::sqrt(2.);
  //PNegB = (pcomB.x0() - evecBasisAB[3].x1() * pcomB.x1() -
  //         evecBasisAB[3].x2() * pcomB.x2() - evecBasisAB[3].x3() * pcomB.x3()) /
  //        std::sqrt(2.);

  ntry = 0;
  foundMass1 = false;
  foundMass2 = false;
  while (((foundMass1 == false) || (foundMass2 == false)) && (ntry < 100)) {
    ntry = ntry + 1;

    makeStringEnds(idqsetA, &idqA1, &idqA2);
    makeStringEnds(idqsetB, &idqB1, &idqB2);

    // idq11 = idqB1;
    // idq12 = idqA2;
    // idq21 = idqA1;
    // idq22 = idqB2;

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
    //xfracA = sample_XVDIS(0., alphapowV, betapowV);
    //xfracB = sample_XVDIS(0., alphapowV, betapowV);
    QTrx = Random::normal(0., sigmaQperp/std::sqrt(2.) );
    QTry = Random::normal(0., sigmaQperp/std::sqrt(2.) );
    QTrn = std::sqrt(QTrx*QTrx + QTry*QTry);
    //QTrn = sample_Qperp(sigmaQperp);
    //phiQ = 2. * M_PI * Random::uniform(0., 1.);
    // fprintf(stderr,"  StringProcess::next_NDiff : xfracA = %e, xfracB = %e\n",
    // xfracA, xfracB);
    // fprintf(stderr,"  StringProcess::next_NDiff : QTrn = %e GeV, phiQ = %e\n",
    // QTrn, phiQ);

    QPos = -pow(fabs(QTrn), 2.) / (2. * xfracB * PNegB);
    QNeg = pow(fabs(QTrn), 2.) / (2. * xfracA * PPosA);

    dPPos = -xfracA * PPosA - QPos;
    dPNeg = xfracB * PNegB - QNeg;

    //pstr1com[0] = (PPosA + dPPos + PNegA + dPNeg) / sqrt(2.);
    //pstr2com[0] = (PPosB - dPPos + PNegB - dPNeg) / sqrt(2.);
    //mstr1 = pow(fabs(pstr1com[0]), 2.);
    //mstr2 = pow(fabs(pstr2com[0]), 2.);
    threeMomentum = evecBasisAB[3] * (PPosA + dPPos - PNegA - dPNeg) / std::sqrt(2.) +
                        evecBasisAB[1] * QTrx + evecBasisAB[2] * QTry;
    pstr1com = FourVector( (PPosA + dPPos + PNegA + dPNeg) / std::sqrt(2.), threeMomentum );
    //pstr1com.set_x0( (PPosA + dPPos + PNegA + dPNeg) / std::sqrt(2.) );
    //pstr1com.set_x1( threeMomentum.x1() );
    //pstr1com.set_x2( threeMomentum.x2() );
    //pstr1com.set_x3( threeMomentum.x3() );
    mstr1 = pstr1com.sqr();
    threeMomentum = evecBasisAB[3] * (PPosB - dPPos - PNegB + dPNeg) / std::sqrt(2.) -
                        evecBasisAB[1] * QTrx - evecBasisAB[2] * QTry;
    pstr2com = FourVector( (PPosB - dPPos + PNegB - dPNeg) / std::sqrt(2.), threeMomentum );
    //pstr2com.set_x0( (PPosB - dPPos + PNegB - dPNeg) / std::sqrt(2.) );
    //pstr2com.set_x1( threeMomentum.x1() );
    //pstr2com.set_x2( threeMomentum.x2() );
    //pstr2com.set_x3( threeMomentum.x3() );
    mstr2 = pstr2com.sqr();
    //for (imu = 1; imu < 4; imu++) {
    //  pstr1com[imu] =
    //      evecBasisAB[3][imu] * (PPosA + dPPos - PNegA - dPNeg) / sqrt(2.) +
    //      evecBasisAB[1][imu] * QTrn * cos(phiQ) +
    //      evecBasisAB[2][imu] * QTrn * sin(phiQ);
    //  pstr2com[imu] =
    //      evecBasisAB[3][imu] * (PPosB - dPPos - PNegB + dPNeg) / sqrt(2.) -
    //      evecBasisAB[1][imu] * QTrn * cos(phiQ) -
    //      evecBasisAB[2][imu] * QTrn * sin(phiQ);
    //  mstr1 = mstr1 - pow(fabs(pstr1com[imu]), 2.);
    //  mstr2 = mstr2 - pow(fabs(pstr2com[imu]), 2.);
    //}
    // fprintf(stderr,"  StringProcess::next_DDiff : pstr1com = (%e, %e, %e, %e)
    // GeV\n",
    //	pstr1com[0], pstr1com[1], pstr1com[2], pstr1com[3]);
    // fprintf(stderr,"  StringProcess::next_DDiff : pstr2com = (%e, %e, %e, %e)
    // GeV\n",
    //	pstr2com[0], pstr2com[1], pstr2com[2], pstr2com[3]);

    mstr1 = (mstr1 > 0.) ? std::sqrt(mstr1) : 0.;
    //if (mstr1 > 0.) {
    //  mstr1 = sqrt(mstr1);
    //} else {
    //  mstr1 = 0.;
    //}
    foundMass1 = mstr1 > (pythia->particleData.m0(idq11) +
                          pythia->particleData.m0(idq12));

    mstr2 = (mstr2 > 0.) ? std::sqrt(mstr2) : 0.;
    //if (mstr2 > 0.) {
    //  mstr2 = sqrt(mstr2);
    //} else {
    //  mstr2 = 0.;
    //}
    foundMass2 = mstr2 > (pythia->particleData.m0(idq21) +
                          pythia->particleData.m0(idq22));
  }

  ret = false;
  bool both_masses_above_pythia_threshold = foundMass1 && foundMass2;
  if ( both_masses_above_pythia_threshold ) {
    // fprintf(stderr,"  StringProcess::next_NDiff : idq11 = %d, idq12 = %d\n", idq11,
    // idq12);
    // fprintf(stderr,"  StringProcess::next_NDiff : mstr1 = %e GeV\n", mstr1);

    // fprintf(stderr,"  StringProcess::next_NDiff : idq21 = %d, idq22 = %d\n", idq21,
    // idq22);
    // fprintf(stderr,"  StringProcess::next_NDiff : mstr2 = %e GeV\n", mstr2);

    //pnull = static_cast<double *>(malloc(4 * sizeof(double)));
    //prs = static_cast<double *>(malloc(4 * sizeof(double)));
    //evec = static_cast<double *>(malloc(4 * sizeof(double)));

    pstr1lab = pstr1com.LorentzBoost( -vcomAB );
    pstr2lab = pstr2com.LorentzBoost( -vcomAB );
    //lorentz->Boost1_RestToLab(3, ucomAB, pstr1com, pstr1lab);
    //lorentz->Boost1_RestToLab(3, ucomAB, pstr2com, pstr2lab);
    ustr1com = pstr1com / mstr1;
    ustr2com = pstr2com / mstr2;
    ustr1lab = pstr1lab / mstr1;
    ustr2lab = pstr2lab / mstr2;
    //for (imu = 0; imu < 4; imu++) {
    //  ustr1com[imu] = pstr1com[imu] / mstr1;
    //  ustr2com[imu] = pstr2com[imu] / mstr2;

    //  ustr1lab[imu] = pstr1lab[imu] / mstr1;
    //  ustr2lab[imu] = pstr2lab[imu] / mstr2;
    //}

    threeMomentum = pstr1com.threevec();
    pnull = FourVector( threeMomentum.abs(), threeMomentum );
    //pnull.set_x0( threeMomentum.abs() );
    //pnull.set_x1( threeMomentum.x1() );
    //pnull.set_x2( threeMomentum.x2() );
    //pnull.set_x3( threeMomentum.x3() );
    //for (imu = 1; imu < 4; imu++) {
    //  pnull[imu] = pstr1com[imu];
    //}
    //pnull[0] = sqrt(pow(fabs(pnull[1]), 2.) + pow(fabs(pnull[2]), 2.) +
    //                pow(fabs(pnull[3]), 2.));
    prs = pnull.LorentzBoost( ustr1com.velocity() );
    pabs = prs.threevec().abs();
    //lorentz->Boost1_LabToRest(3, ustr1com, pnull, prs);
    //pabs = sqrt(pow(fabs(prs[1]), 2.) + pow(fabs(prs[2]), 2.) +
    //            pow(fabs(prs[3]), 2.));
    evec = prs.threevec() / pabs;
    //evec[0] = 0.;
    //for (imu = 1; imu < 4; imu++) {
    //  evec[imu] = prs[imu] / pabs;
    //}
    // fprintf(stderr,"  StringProcess::next_NDiff : evec = (%e, %e, %e)\n", evec[1],
    // evec[2], evec[3]);
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
    //pnull.set_x0( threeMomentum.abs() );
    //pnull.set_x1( threeMomentum.x1() );
    //pnull.set_x2( threeMomentum.x2() );
    //pnull.set_x3( threeMomentum.x3() );
    //for (imu = 1; imu < 4; imu++) {
    //  pnull[imu] = pstr2com[imu];
    //}
    //pnull[0] = sqrt(pow(fabs(pnull[1]), 2.) + pow(fabs(pnull[2]), 2.) +
    //                pow(fabs(pnull[3]), 2.));
    prs = pnull.LorentzBoost( ustr2com.velocity() );
    pabs = prs.threevec().abs();
    //lorentz->Boost1_LabToRest(3, ustr2com, pnull, prs);
    //pabs = sqrt(pow(fabs(prs[1]), 2.) + pow(fabs(prs[2]), 2.) +
    //            pow(fabs(prs[3]), 2.));
    evec = prs.threevec() / pabs;
    //evec[0] = 0.;
    //for (imu = 1; imu < 4; imu++) {
    //  evec[imu] = prs[imu] / pabs;
    //}
    // fprintf(stderr,"  StringProcess::next_NDiff : evec = (%e, %e, %e)\n", evec[1],
    // evec[2], evec[3]);
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

    //free(pnull);
    //free(prs);
    //free(evec);
  }

  //free(pstr1com);
  //free(pstr1lab);
  //free(pstr2com);
  //free(pstr2lab);

  //free(ustr1com);
  //free(ustr1lab);
  //free(ustr2com);
  //free(ustr2lab);

  //if (ret == true) {
  //  conserved = check_conservation();
  //  ret = ret && conserved;
  //}
  return ret;
}

// baryon-antibaryon annihilation
bool StringProcess::next_BBbarAnn(){
	bool ret;
	//bool conserved;

	//int imu;
	int ntry;
	bool isBBbarpair, isAnnihilating;

	int ic, jc;
	int ijc, ipr, npr;
	double rflip;
	std::vector<int> indexAnn;

	int idq11, idq12, idq12prev;
	int idq21, idq22, idq22prev;
	int nfrag1, nfrag2;
	double mstr1, mstr2;
	double mstr1Min, mstr2Min;

	FourVector ustr1lab;
	FourVector ustr2lab;
	//double *ustr1lab;
	//double *ustr2lab;

	double pabs;
	ThreeVector evec;
	//double *evec;

	ustr1lab = ucomAB;
	ustr2lab = ucomAB;
	//ustr1lab = static_cast<double *>( malloc(4*sizeof(double)) );
	//ustr2lab = static_cast<double *>( malloc(4*sizeof(double)) );
	//for(imu = 0; imu < 4; imu++){
	//	ustr1lab[imu] = ucomAB[imu];
	//	ustr2lab[imu] = ucomAB[imu];
	//}

	//indexAnn = new vector<int>;
	indexAnn.resize(0);

	reset_finalArray();

	isBBbarpair = ( (baryonA == 3) && (baryonB == -3) ) || ( (baryonA == -3) && (baryonB == 3) );
	isAnnihilating = false;

	// if it is baryon-antibaryon pair
	if( isBBbarpair == true ){ // if it is

		//PDGid2idqset(PDGidA, idqsetA);
		//PDGid2idqset(PDGidB, idqsetB);

		mstr1 = 0.5*sqrtsAB;
		mstr2 = 0.5*sqrtsAB;

		for(ic = 1; ic <= 3; ic++){
			for(jc = 1; jc <= 3; jc++){
				if( idqsetA[ic] == -idqsetB[jc] ){
					ijc = ic*10 + jc;
					indexAnn.push_back( ijc );
				}
			}
		}
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
			ipr = static_cast<int>( floor( Random::uniform(0., 1.)*static_cast<double>(npr) ) );
			ijc = indexAnn.at(ipr);
			ic = ( ijc - (ijc%10) )/10;
			jc = ijc%10;
			fprintf(stderr,"  StringProcess::next_BBarAnn : ic = %d, jc = %d chosen\n", ic, jc);
			// make two qqbar pairs to excite strings
			if( (baryonA == 3) && (baryonB == -3) ){
				idq11 = idqsetA[1 + ic%3];
				idq12 = idqsetB[1 + jc%3];
				idq21 = idqsetA[1 + (ic + 1)%3];
				idq22 = idqsetB[1 + (jc + 1)%3];
			}
			else if( (baryonA == -3) && (baryonB == 3) ){
				idq11 = idqsetB[1 + jc%3];
				idq12 = idqsetA[1 + ic%3];
				idq21 = idqsetB[1 + (jc + 1)%3];
				idq22 = idqsetA[1 + (ic + 1)%3];
			}
			// randomly choose if we flip the antiquark contents
			rflip = Random::uniform(0., 1.);
			if( rflip > 0.5 ){
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

		//evec = static_cast<double *>( malloc(4*sizeof(double)) );

		// string 1
		pabs = pcomA.threevec().abs();
		evec = pcomA.threevec() / pabs;
		//pabs = sqrt( pow(fabs(pcomA[1]),2.) + pow(fabs(pcomA[2]),2.) + pow(fabs(pcomA[3]),2.) );
		//evec[0] = 0.;
		//for(imu = 1; imu < 4; imu++){
		//	evec[imu] = pcomA[imu]/pabs;
		//}
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
		//pabs = sqrt( pow(fabs(pcomB[1]),2.) + pow(fabs(pcomB[2]),2.) + pow(fabs(pcomB[3]),2.) );
		//evec[0] = 0.;
		//for(imu = 1; imu < 4; imu++){
		//	evec[imu] = pcomB[imu]/pabs;
		//}
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

		//free(evec);
	}

	//free(ustr1lab);
	//free(ustr2lab);

	//delete indexAnn;

	//if( ret == true ){
	//	conserved = check_conservation();
	//	ret = ret && conserved;
	//}
	return ret;
}

/*
void StringProcess::PDGid2idqset(int pdgid, int *idqset) {
  int ic;
  int pdgidnew;

  if (pdgid > 0) {
    pdgidnew = pdgid % 10000;
  } else {
    pdgidnew = abs(pdgid) % 10000;
  }

  if (pdgidnew < 1000) {
    idqset[3] = 0;
  } else {
    idqset[3] = (pdgidnew - (pdgidnew % 1000)) / 1000;
  }

  pdgidnew = pdgidnew % 1000;
  idqset[2] = (pdgidnew - (pdgidnew % 100)) / 100;

  pdgidnew = pdgidnew % 100;
  idqset[1] = (pdgidnew - (pdgidnew % 10)) / 10;

  pdgidnew = pdgidnew % 10;
  idqset[0] = pdgidnew;

  // if it is meson/baryon
  if (idqset[3] == 0) {  // mesons
    if (pdgid < 0) {
      idqset[2] = -idqset[2];
    } else {
      idqset[1] = -idqset[1];
    }
  } else {  // baryons
    if (pdgid < 0) {
      for (ic = 1; ic <= 3; ic++) {
        idqset[ic] = -idqset[ic];
      }
    }
  }  // endif meson/baryon
}
*/

void StringProcess::make_orthonormal_basis(){
  if (fabs(pcomA.x3()) < (1. - 1.0e-8) * pabscomAB) {
    double ex, ey, et;
    double theta, phi;

    evecBasisAB[3] = pcomA.threevec() / pabscomAB;
    //for (imu = 1; imu <= 3; imu++) {
    //  evecBasisAB[3][imu] = pcomA[imu] / pabscomAB;
    //}
    theta = acos(evecBasisAB[3].x3());
    //theta = acos(evecBasisAB[3][3]);

    ex = evecBasisAB[3].x1();
    ey = evecBasisAB[3].x2();
    //ex = evecBasisAB[3][1];
    //ey = evecBasisAB[3][2];
    et = std::sqrt(ex*ex + ey*ey);
    if (ey > 0.) {
      phi = acos(ex / et);
    } else {
      phi = -acos(ex / et);
    }

    evecBasisAB[1].set_x1( cos(theta) * cos(phi) );
    evecBasisAB[1].set_x2( cos(theta) * sin(phi) );
    evecBasisAB[1].set_x3( -sin(theta) );
    //evecBasisAB[1][1] = cos(theta) * cos(phi);
    //evecBasisAB[1][2] = cos(theta) * sin(phi);
    //evecBasisAB[1][3] = -sin(theta);

    evecBasisAB[2].set_x1( -sin(phi) );
    evecBasisAB[2].set_x2( cos(phi) );
    evecBasisAB[2].set_x3( 0. );
    //evecBasisAB[2][1] = -sin(phi);
    //evecBasisAB[2][2] = cos(phi);
    //evecBasisAB[2][3] = 0.;
  } else {
    if (pcomA.x3() > 0.) {
      evecBasisAB[1] = ThreeVector(1., 0., 0.);
      //evecBasisAB[1].set_x1( 1. );
      //evecBasisAB[1].set_x2( 0. );
      //evecBasisAB[1].set_x3( 0. );

      evecBasisAB[2] = ThreeVector(0., 1., 0.);
      //evecBasisAB[2].set_x1( 0. );
      //evecBasisAB[2].set_x2( 1. );
      //evecBasisAB[2].set_x3( 0. );

      evecBasisAB[3] = ThreeVector(0., 0., 1.);
      //evecBasisAB[3].set_x1( 0. );
      //evecBasisAB[3].set_x2( 0. );
      //evecBasisAB[3].set_x3( 1. );
      //for (imu = 1; imu <= 3; imu++) {
      //  for (inu = 1; inu <= 3; inu++) {
      //    if (imu == inu) {
      //      evecBasisAB[imu][inu] = 1.;
      //    } else {
      //      evecBasisAB[imu][inu] = 0.;
      //    }
      //  }
      //}
    } else {
      evecBasisAB[1] = ThreeVector(0., 1., 0.);
      //evecBasisAB[1].set_x1( 0. );
      //evecBasisAB[1].set_x2( 1. );
      //evecBasisAB[1].set_x3( 0. );
      //evecBasisAB[1][1] = 0.;
      //evecBasisAB[1][2] = 1.;
      //evecBasisAB[1][3] = 0.;

      evecBasisAB[2] = ThreeVector(1., 0., 0.);
      //evecBasisAB[2].set_x1( 1. );
      //evecBasisAB[2].set_x2( 0. );
      //evecBasisAB[2].set_x3( 0. );
      //evecBasisAB[2][1] = 1.;
      //evecBasisAB[2][2] = 0.;
      //evecBasisAB[2][3] = 0.;

      evecBasisAB[3] = ThreeVector(0., 0., -1.);
      //evecBasisAB[3].set_x1( 0. );
      //evecBasisAB[3].set_x2( 0. );
      //evecBasisAB[3].set_x3( -1. );
      //evecBasisAB[3][1] = 0.;
      //evecBasisAB[3][2] = 0.;
      //evecBasisAB[3][3] = -1.;
    }
  }
}

void StringProcess::compute_incoming_lightcone_momenta(){
  PPosA = ( pcomA.x0() + evecBasisAB[3] * pcomA.threevec() ) / std::sqrt(2.);
  //PPosA = (pcomA.x0() + evecBasisAB[3].x1() * pcomA.x1() +
  //         evecBasisAB[3].x2() * pcomA.x2() + evecBasisAB[3].x3() * pcomA.x3()) /
  //        std::sqrt(2.);
  PNegA = ( pcomA.x0() - evecBasisAB[3] * pcomA.threevec() ) / std::sqrt(2.);
  //PNegA = (pcomA.x0() - evecBasisAB[3].x1() * pcomA.x1() -
  //         evecBasisAB[3].x2() * pcomA.x2() - evecBasisAB[3].x3() * pcomA.x3()) /
  //        std::sqrt(2.);
  PPosB = ( pcomB.x0() + evecBasisAB[3] * pcomB.threevec() ) / std::sqrt(2.);
  //PPosB = (pcomB.x0() + evecBasisAB[3].x1() * pcomB.x1() +
  //         evecBasisAB[3].x2() * pcomB.x2() + evecBasisAB[3].x3() * pcomB.x3()) /
  //        std::sqrt(2.);
  PNegB = ( pcomB.x0() - evecBasisAB[3] * pcomB.threevec() ) / std::sqrt(2.);
  //PNegB = (pcomB.x0() - evecBasisAB[3].x1() * pcomB.x1() -
  //         evecBasisAB[3].x2() * pcomB.x2() - evecBasisAB[3].x3() * pcomB.x3()) /
  //        std::sqrt(2.);
}

void StringProcess::makeStringEnds(std::array<int,4> &idqset, int *idq1, int *idq2) {
  int ir, ic, jc;
  int idq1tmp, idq2tmp;
  std::array<int,3> idqtmp;
  //int *idqtmp;

  double rspin;
  double rflip;

  //idqtmp = static_cast<int *>(malloc(3 * sizeof(int)));

  // if it is meson/baryon
  if (idqset[3] == 0) {  // meson
    ir = static_cast<int>(floor(1. + 2. * Random::uniform(0., 1.)));

    idq1tmp = idqset[ir];
    jc = 1 + ir % 2;
    idq2tmp = idqset[jc];
  } else {  // baryon
    ir = static_cast<int>(floor(1. + 3. * Random::uniform(0., 1.)));
    idqtmp[0] = idqset[0];

    idq1tmp = idqset[ir];
    for (ic = 0; ic < 2; ic++) {
      jc = 1 + (ir + ic) % 3;
      idqtmp[ic + 1] = abs(idqset[jc]);
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
      if (rspin < 0.25) {
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
    *idq1 = idq1tmp;
    *idq2 = idq2tmp;
  } else {
    *idq1 = idq2tmp;
    *idq2 = idq1tmp;
  }

  /* some mesons with PDG id 11X are actually mixed state of uubar and ddbar.
   * have a random selection whether we have uubar or ddbar in this case. */
  if ((idqset[3] == 0) && (*idq1 == 1) && (*idq2 == -1)) {
    rflip = Random::uniform(0., 1.);
    if (rflip > 0.5) {
      *idq1 = 2;
      *idq2 = -2;
    }
  }

  // fprintf(stderr,"  StringProcess::makeStringEnds : %d/%d/%d/%d decomposed into %d
  // and %d\n",
  //	idqset[3], idqset[2], idqset[1], idqset[0], *idq1, *idq2);

  //free(idqtmp);
}

int StringProcess::fragmentString(int idq1, int idq2, double mString,
                            ThreeVector &evecLong, bool ranrot) {
  int ret;
  bool Hforced;

  int ipart;
  int status;
  int col, acol;
  double pCMquark;
  double m1, m2;
  double phi, theta;
  double rswap;
  Pythia8::Vec4 pquark;

  ThreeVector p3vec;
  FourVector pvRS;
  //double *p3vec;
  //double *pvRS;

  //p3vec = static_cast<double *>(malloc(4 * sizeof(double)));
  //pvRS = static_cast<double *>(malloc(4 * sizeof(double)));

  pythia->event.reset();
  ret = 0;

  m1 = pythia->particleData.m0(idq1);
  m2 = pythia->particleData.m0(idq2);
  pCMquark = pCM( mString, m1, m2 );
  //pCMquark = sqrt((pow(fabs(mString), 2.) - pow(fabs(m1 + m2), 2.)) *
  //           (pow(fabs(mString), 2.) - pow(fabs(m1 - m2), 2.))) /
  //      (2. * mString);

  //p3vec[0] = 0.;
  if (ranrot == true) {
    theta = acos(1. - 2. * Random::uniform(0., 1.));
    phi = 2. * M_PI * Random::uniform(0., 1.);
    rswap = 0.;

    p3vec.set_x1( pCMquark * sin(theta) * cos(phi) );
    //p3vec[1] = pCMquark * sin(theta) * cos(phi);
    p3vec.set_x2( pCMquark * sin(theta) * sin(phi) );
    //p3vec[2] = pCMquark * sin(theta) * sin(phi);
    p3vec.set_x3( pCMquark * cos(theta) );
    //p3vec[3] = pCMquark * cos(theta);
  } else {
    rswap = Random::uniform(0., 1.);
    if (rswap < 0.5) {
      evecLong.set_x1( -evecLong.x1() );
      evecLong.set_x2( -evecLong.x2() );
      evecLong.set_x3( -evecLong.x3() );
    }

    //if (evecLong == NULL) {
    //  p3vec[1] = 0.;
    //  p3vec[2] = 0.;
    //  p3vec[3] = pCMquark;
    //} else {
    //  if (rswap < 0.5) {
    //    evecLong[1] = -evecLong[1];
    //    evecLong[2] = -evecLong[2];
    //    evecLong[3] = -evecLong[3];
    //  }
    //}
  }

  if (m1 + m2 < mString) {
    status = 1;
    col = 1;
    acol = 0;
    if (ranrot == true) {
      pvRS = FourVector(0., -p3vec);
      //pvRS[1] = -p3vec[1];
      //pvRS[2] = -p3vec[2];
      //pvRS[3] = -p3vec[3];
    } else {
      pvRS = FourVector(0., -pCMquark * evecLong);
      //pvRS[1] = -pCMquark * evecLong.x1();
      //pvRS[2] = -pCMquark * evecLong.x2();
      //pvRS[3] = -pCMquark * evecLong.x3();
      //if (evecLong == NULL) {
      //  pvRS[1] = -p3vec[1];
      //  pvRS[2] = -p3vec[2];
      //  pvRS[3] = -p3vec[3];
      //} else {
      //  pvRS[1] = -pCMquark * evecLong[1];
      //  pvRS[2] = -pCMquark * evecLong[2];
      //  pvRS[3] = -pCMquark * evecLong[3];
      //}
    }
    pvRS.set_x0( std::sqrt(m1*m1 + pCMquark*pCMquark) );
    //pvRS[0] = std::sqrt(pow(fabs(m1), 2.) + pow(fabs(pCMquark), 2.));

    pquark.e( pvRS.x0() );
    //pquark.e(pvRS[0]);
    pquark.px( pvRS.x1() );
    //pquark.px(pvRS[1]);
    pquark.py( pvRS.x2() );
    //pquark.py(pvRS[2]);
    pquark.pz( pvRS.x3() );
    //pquark.pz(pvRS[3]);

    pythia->event.append(idq1, status, col, acol, pquark, m1);

    status = 1;
    col = 0;
    acol = 1;
    if (ranrot == true) {
      pvRS = FourVector(0., p3vec);
      //pvRS[1] = p3vec[1];
      //pvRS[2] = p3vec[2];
      //pvRS[3] = p3vec[3];
    } else {
      pvRS = FourVector(0., pCMquark * evecLong);
      //pvRS[1] = pCMquark * evecLong.x1();
      //pvRS[2] = pCMquark * evecLong.x2();
      //pvRS[3] = pCMquark * evecLong.x3();
      //if (evecLong == NULL) {
      //  pvRS[1] = p3vec[1];
      //  pvRS[2] = p3vec[2];
      //  pvRS[3] = p3vec[3];
      //} else {
      //  pvRS[1] = pCMquark * evecLong[1];
      //  pvRS[2] = pCMquark * evecLong[2];
      //  pvRS[3] = pCMquark * evecLong[3];
      //}
    }
    pvRS.set_x0( std::sqrt(m2*m2 + pCMquark*pCMquark) );
    //pvRS[0] = std::sqrt(pow(fabs(m2), 2.) + pow(fabs(pCMquark), 2.));

    pquark.e( pvRS.x0() );
    //pquark.e(pvRS[0]);
    pquark.px( pvRS.x1() );
    //pquark.px(pvRS[1]);
    pquark.py( pvRS.x2() );
    //pquark.py(pvRS[2]);
    pquark.pz( pvRS.x3() );
    //pquark.pz(pvRS[3]);

    pythia->event.append(idq2, status, col, acol, pquark, m2);

    Hforced = pythia->forceHadronLevel();

    if (Hforced == true) {
      for (ipart = 0; ipart < pythia->event.size(); ipart++) {
        if (pythia->event[ipart].isFinal()) {
          ret = ret + 1;
        }
      }
    }
  } else {
    fprintf(stderr, "  StringProcess::fragmentString failure : m1 + m2 >= mString\n");
  }

  //free(p3vec);
  //free(pvRS);

  return ret;
}

/*
double StringProcess::sample_XSDIS(double xmin, double b) {
  bool accepted;

  double xfrac = 0.;
  double rx, ra;
  double pdf, env;

  accepted = false;
  while (accepted == false) {
    rx = log(1. / xmin) * Random::uniform(0., 1.);
    xfrac = xmin * exp(rx);

    env = 1. / xfrac;
    pdf = pow(fabs(1. - xfrac), 1. + b) / xfrac;

    ra = env * Random::uniform(0., 1.);
    if (ra < pdf) {
      accepted = true;
    }
  }

  return xfrac;
}
*/

/*
double StringProcess::sample_XVDIS(double xmin, double a, double b) {
  bool accepted;

  double xfrac = 0.;
  double rx, ra;
  double pdf, env;

  accepted = false;
  while (accepted == false) {
    rx = (1. - pow(fabs(xmin), a)) * Random::uniform(0., 1.);
    xfrac = pow(fabs(rx + pow(fabs(xmin), a)), 1. / a);

    env = pow(fabs(xfrac), a - 1.) * pow(fabs(1. - xmin), b - 1.);
    pdf = pow(fabs(xfrac), a - 1.) * pow(fabs(1. - xfrac), b - 1.);

    ra = env * Random::uniform(0., 1.);
    if (ra < pdf) {
      accepted = true;
    }
  }

  return xfrac;
}
*/

/*
double StringProcess::sample_Qperp(double sigQ) {
  double Qperp = 0.;
  double rq;

  rq = Random::uniform(0., 1.);
  Qperp = sigQ * sqrt(fabs(log(1. / (1. - rq))));

  return Qperp;
}
*/

}  // namespace Smash
