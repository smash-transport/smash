#include "include/spmerge.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// constructor
Lorentz::Lorentz() {}

// destructor
Lorentz::~Lorentz() {}

void Lorentz::Boost1_RestToLab(int dim, double *U, double *Vin, double *Vout) {
  int i, ii, j;
  // double **Lambda;
  Lambda = static_cast<double **>(malloc((dim + 1) * sizeof(double *)));
  for (i = 0; i <= dim; i++) {
    Lambda[i] = static_cast<double *>(malloc((dim + 1) * sizeof(double)));
  }
  Lambda[0][0] = U[0];
  for (i = 1; i <= dim; i++) {
    Lambda[i][0] = U[i];
    Lambda[0][i] = U[i];
    for (j = 1; j <= dim; j++) {
      if (i == j) {
        Lambda[i][j] = 1. + U[i] * U[j] / (U[0] + 1.);
      } else {
        Lambda[i][j] = U[i] * U[j] / (U[0] + 1.);
      }
    }
  }
  for (i = 0; i <= dim; i++) {
    Vout[i] = 0.;
    for (ii = 0; ii <= dim; ii++) {
      Vout[i] = Vout[i] + Lambda[i][ii] * Vin[ii];
    }
  }
  for (i = 0; i <= dim; i++) {
    free(Lambda[i]);
  }
  free(Lambda);
}

void Lorentz::Boost1_LabToRest(int dim, double *U, double *Vin, double *Vout) {
  int i, ii, j;
  // double **Lambda;
  Lambda = static_cast<double **>(malloc((dim + 1) * sizeof(double *)));
  for (i = 0; i <= dim; i++) {
    Lambda[i] = static_cast<double *>(malloc((dim + 1) * sizeof(double)));
  }
  Lambda[0][0] = U[0];
  for (i = 1; i <= dim; i++) {
    Lambda[i][0] = -U[i];
    Lambda[0][i] = -U[i];
    for (j = 1; j <= dim; j++) {
      if (i == j) {
        Lambda[i][j] = 1. + U[i] * U[j] / (U[0] + 1.);
      } else {
        Lambda[i][j] = U[i] * U[j] / (U[0] + 1.);
      }
    }
  }
  for (i = 0; i <= dim; i++) {
    Vout[i] = 0.;
    for (ii = 0; ii <= dim; ii++) {
      Vout[i] = Vout[i] + Lambda[i][ii] * Vin[ii];
    }
  }
  for (i = 0; i <= dim; i++) {
    free(Lambda[i]);
  }
  free(Lambda);
}

void Lorentz::Boost2_RestToLab(int dim, double *U, double **Tin,
                               double **Tout) {
  int i, ii, j, jj;
  // double **Lambda;
  Lambda = static_cast<double **>(malloc((dim + 1) * sizeof(double *)));
  for (i = 0; i <= dim; i++) {
    Lambda[i] = static_cast<double *>(malloc((dim + 1) * sizeof(double)));
  }
  Lambda[0][0] = U[0];
  for (i = 1; i <= dim; i++) {
    Lambda[i][0] = U[i];
    Lambda[0][i] = U[i];
    for (j = 1; j <= dim; j++) {
      if (i == j) {
        Lambda[i][j] = 1. + U[i] * U[j] / (U[0] + 1.);
      } else {
        Lambda[i][j] = U[i] * U[j] / (U[0] + 1.);
      }
    }
  }
  for (i = 0; i <= dim; i++) {
    for (j = 0; j <= dim; j++) {
      Tout[i][j] = 0.;
      for (ii = 0; ii <= dim; ii++) {
        for (jj = 0; jj <= dim; jj++) {
          Tout[i][j] = Tout[i][j] + Lambda[i][ii] * Lambda[j][jj] * Tin[ii][jj];
        }
      }
    }
  }
  for (i = 0; i <= dim; i++) {
    free(Lambda[i]);
  }
  free(Lambda);
}

void Lorentz::Boost2_LabToRest(int dim, double *U, double **Tin,
                               double **Tout) {
  int i, ii, j, jj;
  // double **Lambda;
  Lambda = static_cast<double **>(malloc((dim + 1) * sizeof(double *)));
  for (i = 0; i <= dim; i++) {
    Lambda[i] = static_cast<double *>(malloc((dim + 1) * sizeof(double)));
  }
  Lambda[0][0] = U[0];
  for (i = 1; i <= dim; i++) {
    Lambda[i][0] = -U[i];
    Lambda[0][i] = -U[i];
    for (j = 1; j <= dim; j++) {
      if (i == j) {
        Lambda[i][j] = 1. + U[i] * U[j] / (U[0] + 1.);
      } else {
        Lambda[i][j] = U[i] * U[j] / (U[0] + 1.);
      }
    }
  }
  for (i = 0; i <= dim; i++) {
    for (j = 0; j <= dim; j++) {
      Tout[i][j] = 0.;
      for (ii = 0; ii <= dim; ii++) {
        for (jj = 0; jj <= dim; jj++) {
          Tout[i][j] = Tout[i][j] + Lambda[i][ii] * Lambda[j][jj] * Tin[ii][jj];
        }
      }
    }
  }
  for (i = 0; i <= dim; i++) {
    free(Lambda[i]);
  }
  free(Lambda);
}

void Lorentz::TransRotation1(int dim, double phi, double *Vin, double *Vout) {
  int i, ii, j;
  // double **Mrot;
  Mrot = static_cast<double **>(malloc((dim + 1) * sizeof(double *)));
  for (i = 0; i <= dim; i++) {
    Mrot[i] = static_cast<double *>(malloc((dim + 1) * sizeof(double)));
    for (j = 0; j <= dim; j++) {
      Mrot[i][j] = 0.;
    }
  }

  Mrot[0][0] = 1.;
  for (i = 1; i <= dim; i++) {
    if (i > 2) {
      Mrot[i][i] = 1.;
    }
  }
  Mrot[1][1] = cos(phi);
  Mrot[1][2] = -sin(phi);
  Mrot[2][1] = sin(phi);
  Mrot[2][2] = cos(phi);

  for (i = 0; i <= dim; i++) {
    Vout[i] = 0.;
    for (ii = 0; ii <= dim; ii++) {
      Vout[i] = Vout[i] + Mrot[i][ii] * Vin[ii];
    }
  }

  for (i = 0; i <= dim; i++) {
    free(Mrot[i]);
  }
  free(Mrot);
}

void Lorentz::TransRotation2(int dim, double phi, double **Tin, double **Tout) {
  int i, ii, j, jj;
  // double **Mrot;
  Mrot = static_cast<double **>(malloc((dim + 1) * sizeof(double *)));
  for (i = 0; i <= dim; i++) {
    Mrot[i] = static_cast<double *>(malloc((dim + 1) * sizeof(double)));
    for (j = 0; j <= dim; j++) {
      Mrot[i][j] = 0.;
    }
  }

  Mrot[0][0] = 1.;
  for (i = 1; i <= dim; i++) {
    if (i > 2) {
      Mrot[i][i] = 1.;
    }
  }
  Mrot[1][1] = cos(phi);
  Mrot[1][2] = -sin(phi);
  Mrot[2][1] = sin(phi);
  Mrot[2][2] = cos(phi);

  for (i = 0; i <= dim; i++) {
    for (j = 0; j <= dim; j++) {
      Tout[i][j] = 0.;
      for (ii = 0; ii <= dim; ii++) {
        for (jj = 0; jj <= dim; jj++) {
          Tout[i][j] = Tout[i][j] + Mrot[i][ii] * Mrot[j][jj] * Tin[ii][jj];
        }
      }
    }
  }

  for (i = 0; i <= dim; i++) {
    free(Mrot[i]);
  }
  free(Mrot);
}

// constructor
SPmerge::SPmerge() {
  int ic, imu, inu, iproc;
  idqsetA = static_cast<int *>(malloc(4 * sizeof(int)));
  idqsetB = static_cast<int *>(malloc(4 * sizeof(int)));
  for (ic = 0; ic < 4; ic++) {
    idqsetA[ic] = 0;
    idqsetB[ic] = 0;
  }

  plabA = static_cast<double *>(malloc(4 * sizeof(double)));
  plabB = static_cast<double *>(malloc(4 * sizeof(double)));
  pcomA = static_cast<double *>(malloc(4 * sizeof(double)));
  pcomB = static_cast<double *>(malloc(4 * sizeof(double)));
  ucomAB = static_cast<double *>(malloc(4 * sizeof(double)));
  evecBasisAB = static_cast<double **>(malloc(4 * sizeof(double *)));
  for (imu = 0; imu < 4; imu++) {
    plabA[imu] = 0.;
    plabB[imu] = 0.;
    pcomA[imu] = 0.;
    pcomB[imu] = 0.;
    ucomAB[imu] = 0.;
    evecBasisAB[imu] = static_cast<double *>(malloc(4 * sizeof(double)));
    for (inu = 0; inu < 4; inu++) {
      evecBasisAB[imu][inu] = 0.;
    }
  }

  //XSecSummed = static_cast<double *>(malloc(5 * sizeof(double)));
  //for (iproc = 0; iproc < 5; iproc++) {
  //  XSecSummed[iproc] = 0.;
  //}

  PDGidA = PDGidB = 0;
  baryonA = baryonB = 0;
  chargeA = chargeB = 0;
  massA = massB = 0.;
  sqrtsAB = 0.;
  pabscomAB = 0.;

  pLightConeMin = 0.001;
  betapowS = 0.5;

  alphapowV = 1.;
  betapowV = 2.5;

  sigmaQperp = 1.2;
  kappaString = 1.;

  BINI = BFIN = 0;
  CINI = CFIN = 0;
  EINI = 0.;
  EFIN = 0.;

  lorentz = new Lorentz();

  final_PDGid = new vector<int>[2];
  final_pvec = new vector<double>[5];
  final_tform = new vector<double>[2];

  reset_finalArray();
}

// destructor
SPmerge::~SPmerge() {
  int imu;

  free(idqsetA);
  free(idqsetB);

  free(plabA);
  free(plabB);
  free(pcomA);
  free(pcomB);
  free(ucomAB);
  for (imu = 0; imu < 4; imu++) {
    free(evecBasisAB[imu]);
  }
  free(evecBasisAB);

  //free(XSecSummed);

  delete lorentz;

  delete[] final_PDGid;
  delete[] final_pvec;
  delete[] final_tform;
}

void SPmerge::set_pythia(Pythia *pythiaIn) {
  pythia = pythiaIn;

  //sigmaTot.init(&pythia->info, pythia->settings, &pythia->particleData);
}

void SPmerge::reset_finalArray() {
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
int SPmerge::append_finalArray(double *uString, double *evecLong) {
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

  int *idfrag;
  double *Efrag;
  double *pxfrag;
  double *pyfrag;
  double *pzfrag;
  double *mfrag;

  int *indY;
  // double zparallel;
  double *pparallel;
  double *Yparallel;
  double *XVertexPos;
  double *XVertexNeg;

  double *pvRS;
  double *pvCM;

  double *xfRS;
  double *xfCM;

  pvRS = static_cast<double *>(malloc(4 * sizeof(double)));
  pvCM = static_cast<double *>(malloc(4 * sizeof(double)));

  xfRS = static_cast<double *>(malloc(4 * sizeof(double)));
  xfCM = static_cast<double *>(malloc(4 * sizeof(double)));

  nfrag = 0;
  bstring = 0;
  for (ipyth = 0; ipyth < pythia->event.size(); ipyth++) {
    if (pythia->event[ipyth].isFinal()) {
      nfrag = nfrag + 1;
      id = pythia->event[ipyth].id();
      bstring = bstring + pythia->particleData.baryonNumberType(id);
    }
  }

  // fprintf(stderr,"  SPmerge::append_finalArray nfrag = %d\n", nfrag);

  idfrag = static_cast<int *>(malloc(nfrag * sizeof(int)));
  Efrag = static_cast<double *>(malloc(nfrag * sizeof(double)));
  pxfrag = static_cast<double *>(malloc(nfrag * sizeof(double)));
  pyfrag = static_cast<double *>(malloc(nfrag * sizeof(double)));
  pzfrag = static_cast<double *>(malloc(nfrag * sizeof(double)));
  mfrag = static_cast<double *>(malloc(nfrag * sizeof(double)));

  indY = static_cast<int *>(malloc(nfrag * sizeof(int)));
  pparallel = static_cast<double *>(malloc(nfrag * sizeof(double)));
  Yparallel = static_cast<double *>(malloc(nfrag * sizeof(double)));
  XVertexPos = static_cast<double *>(malloc((nfrag + 1) * sizeof(double)));
  XVertexNeg = static_cast<double *>(malloc((nfrag + 1) * sizeof(double)));

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

      if (evecLong == NULL) {
        pparallel[ipstr] = pzfrag[ipstr];
        // fprintf(stderr,"  SPmerge::append_finalArray : NULL evecLong\n");
      } else {
        pparallel[ipstr] = pxfrag[ipstr] * evecLong[1] +
                           pyfrag[ipstr] * evecLong[2] +
                           pzfrag[ipstr] * evecLong[3];
        // fprintf(stderr,"  SPmerge::append_finalArray : non-zero evecLong\n");
      }
      Yparallel[ipstr] = 0.5 * log((Efrag[ipstr] + pparallel[ipstr]) /
                                   (Efrag[ipstr] - pparallel[ipstr]));

      pPosTot = pPosTot + (Efrag[ipstr] + pparallel[ipstr]) / sqrt(2.);
      pNegTot = pNegTot + (Efrag[ipstr] - pparallel[ipstr]) / sqrt(2.);

      ipstr = ipstr + 1;
    }
  }
  // fprintf(stderr,"  SPmerge::append_finalArray pPosTot = %e\n", pPosTot);
  // fprintf(stderr,"  SPmerge::append_finalArray pNegTot = %e\n", pNegTot);

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
        (Efrag[ipstr] + pparallel[ipstr]) / (kappaString * sqrt(2.));
  }
  // fprintf(stderr,"  SPmerge::append_finalArray XVertexPos[nfrag] = %e\n",
  // XVertexPos[nfrag]);

  XVertexNeg[nfrag] = pNegTot / kappaString;
  for (kpstr = nfrag - 1; kpstr >= 0; kpstr--) {
    ipstr = indY[kpstr];

    XVertexNeg[kpstr] =
        XVertexNeg[kpstr + 1] -
        (Efrag[ipstr] - pparallel[ipstr]) / (kappaString * sqrt(2.));
  }
  // fprintf(stderr,"  SPmerge::append_finalArray XVertexNeg[0] = %e\n",
  // XVertexNeg[0]);

  ret = 0;
  foundFW = false;
  foundBW = false;
  for (kpstr = 0; kpstr < nfrag; kpstr++) {
    ipstr = indY[kpstr];

    id = idfrag[ipstr];
    baryon = pythia->particleData.baryonNumberType(id);
    pvRS[0] = Efrag[ipstr];
    pvRS[1] = pxfrag[ipstr];
    pvRS[2] = pyfrag[ipstr];
    pvRS[3] = pzfrag[ipstr];
    mass = mfrag[ipstr];

    xfRS[0] = (XVertexPos[kpstr] + XVertexNeg[kpstr + 1]) / sqrt(2.);
    if (evecLong == NULL) {
      xfRS[1] = 0.;
      xfRS[2] = 0.;
      xfRS[3] = (XVertexPos[kpstr] - XVertexNeg[kpstr + 1]) / sqrt(2.);
    } else {
      xfRS[1] =
          evecLong[1] * (XVertexPos[kpstr] - XVertexNeg[kpstr + 1]) / sqrt(2.);
      xfRS[2] =
          evecLong[2] * (XVertexPos[kpstr] - XVertexNeg[kpstr + 1]) / sqrt(2.);
      xfRS[3] =
          evecLong[3] * (XVertexPos[kpstr] - XVertexNeg[kpstr + 1]) / sqrt(2.);
    }

    tProd = (XVertexPos[kpstr] + XVertexNeg[kpstr + 1]) / sqrt(2.);
    // zparallel = (XVertexPos[kpstr] - XVertexNeg[kpstr + 1])/sqrt(2.);
    // fprintf(stderr,"  SPmerge::append_finalArray tProd[kpstr = %d] = %e\n",
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
                  "  SPmerge::append_finalArray warning : particle is not "
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
                  "  SPmerge::append_finalArray warning : particle is not "
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
              "  SPmerge::append_finalArray warning : string is neither "
              "mesonic nor baryonic.\n");
      islead = 0;
      xtotfac = 0.;
    }

    // fprintf(stderr,"  SPmerge::append_finalArray adding to final arrays\n");

    final_PDGid[0].push_back(id);
    final_PDGid[1].push_back(islead);

    if (uString == NULL) {
      E = pvRS[0];
      px = pvRS[1];
      py = pvRS[2];
      pz = pvRS[3];

      tProd = xfRS[0];
      // fprintf(stderr,"  SPmerge::append_finalArray : NULL uString\n");
    } else {
      lorentz->Boost1_RestToLab(3, uString, pvRS, pvCM);
      E = pvCM[0];
      px = pvCM[1];
      py = pvCM[2];
      pz = pvCM[3];

      lorentz->Boost1_RestToLab(3, uString, xfRS, xfCM);
      tProd = xfCM[0];
      // fprintf(stderr,"  SPmerge::append_finalArray : non-zero uString\n");
    }

    final_pvec[0].push_back(E);
    final_pvec[1].push_back(px);
    final_pvec[2].push_back(py);
    final_pvec[3].push_back(pz);
    final_pvec[4].push_back(mass);

    final_tform[0].push_back(tProd);
    final_tform[1].push_back(xtotfac);

    ret = ret + 1;
  }

  free(idfrag);
  free(Efrag);
  free(pxfrag);
  free(pyfrag);
  free(pzfrag);
  free(mfrag);

  free(indY);
  free(pparallel);
  free(Yparallel);
  free(XVertexPos);
  free(XVertexNeg);

  free(pvRS);
  free(pvCM);

  free(xfRS);
  free(xfCM);

  return ret;
}

// energy and charge conservation check
bool SPmerge::check_conservation() {
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
            "  SPmerge::check_conservation : charge violation - CINI = %d, "
            "CFIN = %d\n",
            CINI, CFIN);
    ret = false;
  }

  if (BFIN != BINI) {
    fprintf(stderr,
            "  SPmerge::check_conservation : baryon number violation - BINI = "
            "%d, BFIN = %d\n",
            BINI, BFIN);
    ret = false;
  }

  if (fabs(EFIN - EINI) > 1.0e-8 * EINI) {
    fprintf(stderr,
            "  SPmerge::check_conservation : energy violation - EINI = %e GeV, "
            "EFIN = %e GeV\n",
            EINI, EFIN);
    ret = false;
  }

  if (fabs(pxFIN - pxINI) > 1.0e-8 * EINI) {
    fprintf(stderr,
            "  SPmerge::check_conservation : momentum x violation - pxINI = %e "
            "GeV, pxFIN = %e GeV\n",
            pxINI, pxFIN);
    ret = false;
  }

  if (fabs(pyFIN - pyINI) > 1.0e-8 * EINI) {
    fprintf(stderr,
            "  SPmerge::check_conservation : momentum y violation - pyINI = %e "
            "GeV, pyFIN = %e GeV\n",
            pyINI, pyFIN);
    ret = false;
  }

  if (fabs(pzFIN - pzINI) > 1.0e-8 * EINI) {
    fprintf(stderr,
            "  SPmerge::check_conservation : momentum z violation - pzINI = %e "
            "GeV, pzFIN = %e GeV\n",
            pzINI, pzFIN);
    ret = false;
  }

  return ret;
}

bool SPmerge::init_lab(int idAIn, int idBIn, double massAIn, double massBIn,
                       Vec4 plabAIn, Vec4 plabBIn) {
  bool ret;
  int imu, inu;

  double E, px, py, pz;
  double ex, ey, et;
  double theta, phi;

  PDGidA = idAIn;
  PDGidB = idBIn;
  massA = massAIn;
  massB = massBIn;

  plabA[0] = plabAIn.e();
  plabA[1] = plabAIn.px();
  plabA[2] = plabAIn.py();
  plabA[3] = plabAIn.pz();

  plabB[0] = plabBIn.e();
  plabB[1] = plabBIn.px();
  plabB[2] = plabBIn.py();
  plabB[3] = plabBIn.pz();

  E = plabA[0] + plabB[0];
  px = plabA[1] + plabB[1];
  py = plabA[2] + plabB[2];
  pz = plabA[3] + plabB[3];
  sqrtsAB = sqrt(pow(fabs(E), 2.) - pow(fabs(px), 2.) - pow(fabs(py), 2.) -
                 pow(fabs(pz), 2.));
  pabscomAB = sqrt((pow(fabs(sqrtsAB), 2.) - pow(fabs(massA + massB), 2.)) *
                   (pow(fabs(sqrtsAB), 2.) - pow(fabs(massA - massB), 2.))) /
              (2. * sqrtsAB);
  ucomAB[0] = E / sqrtsAB;
  ucomAB[1] = px / sqrtsAB;
  ucomAB[2] = py / sqrtsAB;
  ucomAB[3] = pz / sqrtsAB;

  lorentz->Boost1_LabToRest(3, ucomAB, plabA, pcomA);
  lorentz->Boost1_LabToRest(3, ucomAB, plabB, pcomB);

  if (fabs(pcomA[3]) < (1. - 1.0e-8) * pabscomAB) {
    for (imu = 1; imu <= 3; imu++) {
      evecBasisAB[3][imu] = pcomA[imu] / pabscomAB;
    }
    theta = acos(evecBasisAB[3][3]);

    ex = evecBasisAB[3][1];
    ey = evecBasisAB[3][2];
    et = sqrt(pow(fabs(ex), 2.) + pow(fabs(ey), 2.));
    if (ey > 0.) {
      phi = acos(ex / et);
    } else {
      phi = -acos(ex / et);
    }

    evecBasisAB[1][1] = cos(theta) * cos(phi);
    evecBasisAB[1][2] = cos(theta) * sin(phi);
    evecBasisAB[1][3] = -sin(theta);

    evecBasisAB[2][1] = -sin(phi);
    evecBasisAB[2][2] = cos(phi);
    evecBasisAB[2][3] = 0.;
  } else {
    if (pcomA[3] > 0.) {
      for (imu = 1; imu <= 3; imu++) {
        for (inu = 1; inu <= 3; inu++) {
          if (imu == inu) {
            evecBasisAB[imu][inu] = 1.;
          } else {
            evecBasisAB[imu][inu] = 0.;
          }
        }
      }
    } else {
      evecBasisAB[1][1] = 0.;
      evecBasisAB[1][2] = 1.;
      evecBasisAB[1][3] = 0.;

      evecBasisAB[2][1] = 1.;
      evecBasisAB[2][2] = 0.;
      evecBasisAB[2][3] = 0.;

      evecBasisAB[3][1] = 0.;
      evecBasisAB[3][2] = 0.;
      evecBasisAB[3][3] = -1.;
    }
  }

  // fprintf(stderr,"  SPmerge::init : evecBasisAB1 = (%e, %e, %e)\n",
  //	evecBasisAB[1][1], evecBasisAB[1][2], evecBasisAB[1][3]);
  // fprintf(stderr,"  SPmerge::init : evecBasisAB2 = (%e, %e, %e)\n",
  //	evecBasisAB[2][1], evecBasisAB[2][2], evecBasisAB[2][3]);
  // fprintf(stderr,"  SPmerge::init : evecBasisAB3 = (%e, %e, %e)\n",
  //	evecBasisAB[3][1], evecBasisAB[3][2], evecBasisAB[3][3]);

  PDGid2idqset(PDGidA, idqsetA);
  PDGid2idqset(PDGidB, idqsetB);

  xfracMin = pLightConeMin / sqrtsAB;

  // quantum numbers of hadron A
  baryonA = pythia->particleData.baryonNumberType(idqsetA[1]) +
            pythia->particleData.baryonNumberType(idqsetA[2]);
  chargeA = pythia->particleData.chargeType(idqsetA[1]) +
            pythia->particleData.chargeType(idqsetA[2]);
  if (idqsetA[3] != 0) {
    baryonA = baryonA + pythia->particleData.baryonNumberType(idqsetA[3]);
    chargeA = chargeA + pythia->particleData.chargeType(idqsetA[3]);
  }
  // quantum numbers of hadron B
  baryonB = pythia->particleData.baryonNumberType(idqsetB[1]) +
            pythia->particleData.baryonNumberType(idqsetB[2]);
  chargeB = pythia->particleData.chargeType(idqsetB[1]) +
            pythia->particleData.chargeType(idqsetB[2]);
  if (idqsetB[3] != 0) {
    baryonB = baryonB + pythia->particleData.baryonNumberType(idqsetB[3]);
    chargeB = chargeB + pythia->particleData.chargeType(idqsetB[3]);
  }

  BINI = baryonA + baryonB;
  CINI = chargeA + chargeB;
  EINI = plabA[0] + plabB[0];
  pxINI = plabA[1] + plabB[1];
  pyINI = plabA[2] + plabB[2];
  pzINI = plabA[3] + plabB[3];

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

bool SPmerge::init_com(int idAIn, int idBIn, double massAIn, double massBIn,
                       double sqrtsABIn) {
  bool ret;
  int imu, inu;

  double E, px, py, pz;

  PDGidA = idAIn;
  PDGidB = idBIn;
  massA = massAIn;
  massB = massBIn;
  sqrtsAB = sqrtsABIn;
  pabscomAB = sqrt((pow(fabs(sqrtsAB), 2.) - pow(fabs(massA + massB), 2.)) *
                   (pow(fabs(sqrtsAB), 2.) - pow(fabs(massA - massB), 2.))) /
              (2. * sqrtsAB);

  plabA[0] = sqrt(pow(fabs(massA), 2.) + pow(fabs(pabscomAB), 2.));
  plabA[1] = 0.;
  plabA[2] = 0.;
  plabA[3] = pabscomAB;

  plabB[0] = sqrt(pow(fabs(massB), 2.) + pow(fabs(pabscomAB), 2.));
  plabB[1] = 0.;
  plabB[2] = 0.;
  plabB[3] = -pabscomAB;

  E = sqrtsAB;
  px = 0.;
  py = 0.;
  pz = 0.;

  ucomAB[0] = E / sqrtsAB;
  ucomAB[1] = px / sqrtsAB;
  ucomAB[2] = py / sqrtsAB;
  ucomAB[3] = pz / sqrtsAB;

  lorentz->Boost1_LabToRest(3, ucomAB, plabA, pcomA);
  lorentz->Boost1_LabToRest(3, ucomAB, plabB, pcomB);

  for (imu = 1; imu <= 3; imu++) {
    for (inu = 1; inu <= 3; inu++) {
      if (imu == inu) {
        evecBasisAB[imu][inu] = 1.;
      } else {
        evecBasisAB[imu][inu] = 0.;
      }
    }
  }

  // fprintf(stderr,"  SPmerge::init : evecBasisAB1 = (%e, %e, %e)\n",
  //	evecBasisAB[1][1], evecBasisAB[1][2], evecBasisAB[1][3]);
  // fprintf(stderr,"  SPmerge::init : evecBasisAB2 = (%e, %e, %e)\n",
  //	evecBasisAB[2][1], evecBasisAB[2][2], evecBasisAB[2][3]);
  // fprintf(stderr,"  SPmerge::init : evecBasisAB3 = (%e, %e, %e)\n",
  //	evecBasisAB[3][1], evecBasisAB[3][2], evecBasisAB[3][3]);

  PDGid2idqset(PDGidA, idqsetA);
  PDGid2idqset(PDGidB, idqsetB);

  xfracMin = pLightConeMin / sqrtsAB;

  // quantum numbers of hadron A
  baryonA = pythia->particleData.baryonNumberType(idqsetA[1]) +
            pythia->particleData.baryonNumberType(idqsetA[2]);
  chargeA = pythia->particleData.chargeType(idqsetA[1]) +
            pythia->particleData.chargeType(idqsetA[2]);
  if (idqsetA[3] != 0) {
    baryonA = baryonA + pythia->particleData.baryonNumberType(idqsetA[3]);
    chargeA = chargeA + pythia->particleData.chargeType(idqsetA[3]);
  }
  // quantum numbers of hadron B
  baryonB = pythia->particleData.baryonNumberType(idqsetB[1]) +
            pythia->particleData.baryonNumberType(idqsetB[2]);
  chargeB = pythia->particleData.chargeType(idqsetB[1]) +
            pythia->particleData.chargeType(idqsetB[2]);
  if (idqsetB[3] != 0) {
    baryonB = baryonB + pythia->particleData.baryonNumberType(idqsetB[3]);
    chargeB = chargeB + pythia->particleData.chargeType(idqsetB[3]);
  }

  BINI = baryonA + baryonB;
  CINI = chargeA + chargeB;
  EINI = plabA[0] + plabB[0];
  pxINI = plabA[1] + plabB[1];
  pyINI = plabA[2] + plabB[2];
  pzINI = plabA[3] + plabB[3];

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
bool SPmerge::next_Inel() {
  bool ret;

  double rproc;
  rproc = XSecInel * gsl_rng_uniform(rng);
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
bool SPmerge::next_SDiff_AX() {
  bool ret, conserved;

  int imu;
  int ntry;
  bool foundPabsX, foundMassX;

  double mstrMin, mstrMax;
  double pabscomAX, massX, rmass;
  double QTrn, phiQ;

  int nfrag;
  int idqX1, idqX2;

  double *pstrHcom;
  double *pstrHlab;
  double *pstrXcom;
  double *pstrXlab;

  double *ustrHcom;
  double *ustrHlab;
  double *ustrXcom;
  double *ustrXlab;

  double pabs;
  double *pnull;
  double *prs;
  double *evec;

  pstrHcom = static_cast<double *>(malloc(4 * sizeof(double)));
  pstrHlab = static_cast<double *>(malloc(4 * sizeof(double)));
  pstrXcom = static_cast<double *>(malloc(4 * sizeof(double)));
  pstrXlab = static_cast<double *>(malloc(4 * sizeof(double)));

  ustrHcom = static_cast<double *>(malloc(4 * sizeof(double)));
  ustrHlab = static_cast<double *>(malloc(4 * sizeof(double)));
  ustrXcom = static_cast<double *>(malloc(4 * sizeof(double)));
  ustrXlab = static_cast<double *>(malloc(4 * sizeof(double)));

  reset_finalArray();

  ntry = 0;
  foundPabsX = false;
  foundMassX = false;
  mstrMin = massB;
  mstrMax = sqrtsAB - massA;
  if (mstrMin > mstrMax) {
    fprintf(stderr, "  SPmerge::next_SDiff_AX : mstrMin > mstrMax\n");
    fprintf(stderr,
            "  SPmerge::next_SDiff_AX : mstrMin = %e GeV, mstrMax = %e GeV\n",
            mstrMin, mstrMax);
    ntry = 100;
  }
  while (((foundPabsX == false) || (foundMassX == false)) && (ntry < 100)) {
    ntry = ntry + 1;

    makeStringEnds(idqsetB, &idqX1, &idqX2);

    QTrn = sample_Qperp(sigmaQperp);
    phiQ = 2. * M_PI * gsl_rng_uniform(rng);
    // fprintf(stderr,"  SPmerge::next_SDiff_AX : QTrn = %e GeV, phiQ = %e\n",
    // QTrn, phiQ);

    rmass = log(mstrMax / mstrMin) * gsl_rng_uniform(rng);
    massX = mstrMin * exp(rmass);
    pabscomAX = sqrt((pow(fabs(sqrtsAB), 2.) - pow(fabs(massA + massX), 2.)) *
                     (pow(fabs(sqrtsAB), 2.) - pow(fabs(massA - massX), 2.))) /
                (2. * sqrtsAB);

    foundPabsX = pabscomAX > QTrn;
    foundMassX = massX > (pythia->particleData.m0(idqX1) +
                          pythia->particleData.m0(idqX2));
  }

  ret = false;
  if ((foundPabsX == true) && (foundMassX == true)) {
    pstrHcom[0] = sqrt(pow(fabs(pabscomAX), 2.) + pow(fabs(massA), 2.));
    pstrXcom[0] = sqrt(pow(fabs(pabscomAX), 2.) + pow(fabs(massX), 2.));
    for (imu = 1; imu < 4; imu++) {
      pstrHcom[imu] = evecBasisAB[3][imu] *
                          sqrt(pow(fabs(pabscomAX), 2.) - pow(fabs(QTrn), 2.)) +
                      evecBasisAB[1][imu] * QTrn * cos(phiQ) +
                      evecBasisAB[2][imu] * QTrn * sin(phiQ);
      pstrXcom[imu] = -evecBasisAB[3][imu] *
                          sqrt(pow(fabs(pabscomAX), 2.) - pow(fabs(QTrn), 2.)) -
                      evecBasisAB[1][imu] * QTrn * cos(phiQ) -
                      evecBasisAB[2][imu] * QTrn * sin(phiQ);
    }
    // fprintf(stderr,"  SPmerge::next_SDiff_AX : pstrXcom = (%e, %e, %e, %e)
    // GeV\n",
    //	pstrXcom[0], pstrXcom[1], pstrXcom[2], pstrXcom[3]);

    // fprintf(stderr,"  SPmerge::next_SDiff_AX : idqX1 = %d, idqX2 = %d\n",
    // idqX1, idqX2);
    // fprintf(stderr,"  SPmerge::next_SDiff_AX : massX = %e GeV\n", massX);

    pnull = static_cast<double *>(malloc(4 * sizeof(double)));
    prs = static_cast<double *>(malloc(4 * sizeof(double)));
    evec = static_cast<double *>(malloc(4 * sizeof(double)));

    lorentz->Boost1_RestToLab(3, ucomAB, pstrHcom, pstrHlab);
    lorentz->Boost1_RestToLab(3, ucomAB, pstrXcom, pstrXlab);
    for (imu = 0; imu < 4; imu++) {
      ustrHcom[imu] = pstrHcom[imu] / massA;
      ustrXcom[imu] = pstrXcom[imu] / massX;

      ustrHlab[imu] = pstrHlab[imu] / massA;
      ustrXlab[imu] = pstrXlab[imu] / massX;
    }

    for (imu = 1; imu < 4; imu++) {
      pnull[imu] = pstrXcom[imu];
    }
    pnull[0] = sqrt(pow(fabs(pnull[1]), 2.) + pow(fabs(pnull[2]), 2.) +
                    pow(fabs(pnull[3]), 2.));
    lorentz->Boost1_LabToRest(3, ustrXcom, pnull, prs);
    pabs = sqrt(pow(fabs(prs[1]), 2.) + pow(fabs(prs[2]), 2.) +
                pow(fabs(prs[3]), 2.));
    evec[0] = 0.;
    for (imu = 1; imu < 4; imu++) {
      evec[imu] = prs[imu] / pabs;
    }
    // fprintf(stderr,"  SPmerge::next_SDiff_AX : evec = (%e, %e, %e)\n",
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
    final_pvec[0].push_back(pstrHlab[0]);
    final_pvec[1].push_back(pstrHlab[1]);
    final_pvec[2].push_back(pstrHlab[2]);
    final_pvec[3].push_back(pstrHlab[3]);
    final_pvec[4].push_back(massA);
    final_tform[0].push_back(0.);
    final_tform[1].push_back(1.);

    if ((NpartString1 > 0) && (NpartString2 > 0) && (nfrag == NpartString1)) {
      NpartFinal = NpartString1 + NpartString2;
      ret = true;
    }

    free(pnull);
    free(prs);
    free(evec);
  }

  free(pstrHcom);
  free(pstrHlab);
  free(pstrXcom);
  free(pstrXlab);

  free(ustrHcom);
  free(ustrHlab);
  free(ustrXcom);
  free(ustrXlab);

  if (ret == true) {
    conserved = check_conservation();
    ret = ret && conserved;
  }
  return ret;
}

// single diffractive AB > XB
bool SPmerge::next_SDiff_XB() {
  bool ret, conserved;

  int imu;
  int ntry;
  bool foundPabsX, foundMassX;

  double mstrMin, mstrMax;
  double pabscomXB, massX, rmass;
  double QTrn, phiQ;

  int nfrag;
  int idqX1, idqX2;

  double *pstrHcom;
  double *pstrHlab;
  double *pstrXcom;
  double *pstrXlab;

  double *ustrHcom;
  double *ustrHlab;
  double *ustrXcom;
  double *ustrXlab;

  double pabs;
  double *pnull;
  double *prs;
  double *evec;

  pstrHcom = static_cast<double *>(malloc(4 * sizeof(double)));
  pstrHlab = static_cast<double *>(malloc(4 * sizeof(double)));
  pstrXcom = static_cast<double *>(malloc(4 * sizeof(double)));
  pstrXlab = static_cast<double *>(malloc(4 * sizeof(double)));

  ustrHcom = static_cast<double *>(malloc(4 * sizeof(double)));
  ustrHlab = static_cast<double *>(malloc(4 * sizeof(double)));
  ustrXcom = static_cast<double *>(malloc(4 * sizeof(double)));
  ustrXlab = static_cast<double *>(malloc(4 * sizeof(double)));

  reset_finalArray();

  ntry = 0;
  foundPabsX = false;
  foundMassX = false;
  mstrMin = massA;
  mstrMax = sqrtsAB - massB;
  if (mstrMin > mstrMax) {
    fprintf(stderr, "  SPmerge::next_SDiff_XB : mstrMin > mstrMax\n");
    fprintf(stderr,
            "  SPmerge::next_SDiff_XB : mstrMin = %e GeV, mstrMax = %e GeV\n",
            mstrMin, mstrMax);
    ntry = 100;
  }
  while (((foundPabsX == false) || (foundMassX == false)) && (ntry < 100)) {
    ntry = ntry + 1;

    makeStringEnds(idqsetA, &idqX1, &idqX2);

    QTrn = sample_Qperp(sigmaQperp);
    phiQ = 2. * M_PI * gsl_rng_uniform(rng);
    // fprintf(stderr,"  SPmerge::next_SDiff_XB : QTrn = %e GeV, phiQ = %e\n",
    // QTrn, phiQ);

    rmass = log(mstrMax / mstrMin) * gsl_rng_uniform(rng);
    massX = mstrMin * exp(rmass);
    pabscomXB = sqrt((pow(fabs(sqrtsAB), 2.) - pow(fabs(massB + massX), 2.)) *
                     (pow(fabs(sqrtsAB), 2.) - pow(fabs(massB - massX), 2.))) /
                (2. * sqrtsAB);

    foundPabsX = pabscomXB > QTrn;
    foundMassX = massX > (pythia->particleData.m0(idqX1) +
                          pythia->particleData.m0(idqX2));
  }

  ret = false;
  if ((foundPabsX == true) && (foundMassX == true)) {
    pstrHcom[0] = sqrt(pow(fabs(pabscomXB), 2.) + pow(fabs(massB), 2.));
    pstrXcom[0] = sqrt(pow(fabs(pabscomXB), 2.) + pow(fabs(massX), 2.));
    for (imu = 1; imu < 4; imu++) {
      pstrHcom[imu] = -evecBasisAB[3][imu] *
                          sqrt(pow(fabs(pabscomXB), 2.) - pow(fabs(QTrn), 2.)) -
                      evecBasisAB[1][imu] * QTrn * cos(phiQ) -
                      evecBasisAB[2][imu] * QTrn * sin(phiQ);
      pstrXcom[imu] = evecBasisAB[3][imu] *
                          sqrt(pow(fabs(pabscomXB), 2.) - pow(fabs(QTrn), 2.)) +
                      evecBasisAB[1][imu] * QTrn * cos(phiQ) +
                      evecBasisAB[2][imu] * QTrn * sin(phiQ);
    }
    // fprintf(stderr,"  SPmerge::next_SDiff_XB : pstrXcom = (%e, %e, %e, %e)
    // GeV\n",
    //	pstrXcom[0], pstrXcom[1], pstrXcom[2], pstrXcom[3]);

    // fprintf(stderr,"  SPmerge::next_SDiff_XB : idqX1 = %d, idqX2 = %d\n",
    // idqX1, idqX2);
    // fprintf(stderr,"  SPmerge::next_SDiff_XB : massX = %e GeV\n", massX);

    pnull = static_cast<double *>(malloc(4 * sizeof(double)));
    prs = static_cast<double *>(malloc(4 * sizeof(double)));
    evec = static_cast<double *>(malloc(4 * sizeof(double)));

    lorentz->Boost1_RestToLab(3, ucomAB, pstrHcom, pstrHlab);
    lorentz->Boost1_RestToLab(3, ucomAB, pstrXcom, pstrXlab);
    for (imu = 0; imu < 4; imu++) {
      ustrHcom[imu] = pstrHcom[imu] / massB;
      ustrXcom[imu] = pstrXcom[imu] / massX;

      ustrHlab[imu] = pstrHlab[imu] / massB;
      ustrXlab[imu] = pstrXlab[imu] / massX;
    }

    for (imu = 1; imu < 4; imu++) {
      pnull[imu] = pstrXcom[imu];
    }
    pnull[0] = sqrt(pow(fabs(pnull[1]), 2.) + pow(fabs(pnull[2]), 2.) +
                    pow(fabs(pnull[3]), 2.));
    lorentz->Boost1_LabToRest(3, ustrXcom, pnull, prs);
    pabs = sqrt(pow(fabs(prs[1]), 2.) + pow(fabs(prs[2]), 2.) +
                pow(fabs(prs[3]), 2.));
    evec[0] = 0.;
    for (imu = 1; imu < 4; imu++) {
      evec[imu] = prs[imu] / pabs;
    }
    // fprintf(stderr,"  SPmerge::next_SDiff_XB : evec = (%e, %e, %e)\n",
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
    final_pvec[0].push_back(pstrHlab[0]);
    final_pvec[1].push_back(pstrHlab[1]);
    final_pvec[2].push_back(pstrHlab[2]);
    final_pvec[3].push_back(pstrHlab[3]);
    final_pvec[4].push_back(massB);
    final_tform[0].push_back(0.);
    final_tform[1].push_back(1.);

    if ((NpartString1 > 0) && (NpartString2 > 0) && (nfrag == NpartString1)) {
      NpartFinal = NpartString1 + NpartString2;
      ret = true;
    }

    free(pnull);
    free(prs);
    free(evec);
  }

  free(pstrHcom);
  free(pstrHlab);
  free(pstrXcom);
  free(pstrXlab);

  free(ustrHcom);
  free(ustrHlab);
  free(ustrXcom);
  free(ustrXlab);

  if (ret == true) {
    conserved = check_conservation();
    ret = ret && conserved;
  }
  return ret;
}

// double diffractive AB > XX
bool SPmerge::next_DDiff_XX() {
  bool ret, conserved;

  int imu;
  int ntry;
  bool foundMass1, foundMass2;

  double xfracA, xfracB;
  double PPosA, PNegA;
  double PPosB, PNegB;
  double QPos, QNeg, QTrn, phiQ;

  int nfrag1, nfrag2;
  int idq11, idq12;
  int idq21, idq22;
  double mstr1, mstr2;

  double *pstr1com;
  double *pstr1lab;
  double *pstr2com;
  double *pstr2lab;

  double *ustr1com;
  double *ustr1lab;
  double *ustr2com;
  double *ustr2lab;

  double pabs;
  double *pnull;
  double *prs;
  double *evec;

  pstr1com = static_cast<double *>(malloc(4 * sizeof(double)));
  pstr1lab = static_cast<double *>(malloc(4 * sizeof(double)));
  pstr2com = static_cast<double *>(malloc(4 * sizeof(double)));
  pstr2lab = static_cast<double *>(malloc(4 * sizeof(double)));

  ustr1com = static_cast<double *>(malloc(4 * sizeof(double)));
  ustr1lab = static_cast<double *>(malloc(4 * sizeof(double)));
  ustr2com = static_cast<double *>(malloc(4 * sizeof(double)));
  ustr2lab = static_cast<double *>(malloc(4 * sizeof(double)));

  reset_finalArray();

  // PDGid2idqset(PDGidA, idqsetA);
  // PDGid2idqset(PDGidB, idqsetB);
  PPosA = (pcomA[0] + evecBasisAB[3][1] * pcomA[1] +
           evecBasisAB[3][2] * pcomA[2] + evecBasisAB[3][3] * pcomA[3]) /
          sqrt(2.);
  PNegA = (pcomA[0] - evecBasisAB[3][1] * pcomA[1] -
           evecBasisAB[3][2] * pcomA[2] - evecBasisAB[3][3] * pcomA[3]) /
          sqrt(2.);
  PPosB = (pcomB[0] + evecBasisAB[3][1] * pcomB[1] +
           evecBasisAB[3][2] * pcomB[2] + evecBasisAB[3][3] * pcomB[3]) /
          sqrt(2.);
  PNegB = (pcomB[0] - evecBasisAB[3][1] * pcomB[1] -
           evecBasisAB[3][2] * pcomB[2] - evecBasisAB[3][3] * pcomB[3]) /
          sqrt(2.);

  ntry = 0;
  foundMass1 = false;
  foundMass2 = false;
  while (((foundMass1 == false) || (foundMass2 == false)) && (ntry < 100)) {
    ntry = ntry + 1;

    makeStringEnds(idqsetA, &idq11, &idq12);
    makeStringEnds(idqsetB, &idq21, &idq22);

    xfracA = sample_XSDIS(xfracMin, betapowS);
    xfracB = sample_XSDIS(xfracMin, betapowS);
    QTrn = sample_Qperp(sigmaQperp);
    phiQ = 2. * M_PI * gsl_rng_uniform(rng);
    // fprintf(stderr,"  SPmerge::next_DDiff : xfracA = %e, xfracB = %e\n",
    // xfracA, xfracB);
    // fprintf(stderr,"  SPmerge::next_DDiff : QTrn = %e GeV, phiQ = %e\n",
    // QTrn, phiQ);

    QPos = -pow(fabs(QTrn), 2.) / (2. * xfracB * PNegB);
    QNeg = pow(fabs(QTrn), 2.) / (2. * xfracA * PPosA);

    pstr1com[0] = (PPosA + QPos + PNegA + QNeg) / sqrt(2.);
    pstr2com[0] = (PPosB - QPos + PNegB - QNeg) / sqrt(2.);
    mstr1 = pow(fabs(pstr1com[0]), 2.);
    mstr2 = pow(fabs(pstr2com[0]), 2.);
    for (imu = 1; imu < 4; imu++) {
      pstr1com[imu] =
          evecBasisAB[3][imu] * (PPosA + QPos - PNegA - QNeg) / sqrt(2.) +
          evecBasisAB[1][imu] * QTrn * cos(phiQ) +
          evecBasisAB[2][imu] * QTrn * sin(phiQ);
      pstr2com[imu] =
          evecBasisAB[3][imu] * (PPosB - QPos - PNegB + QNeg) / sqrt(2.) -
          evecBasisAB[1][imu] * QTrn * cos(phiQ) -
          evecBasisAB[2][imu] * QTrn * sin(phiQ);
      mstr1 = mstr1 - pow(fabs(pstr1com[imu]), 2.);
      mstr2 = mstr2 - pow(fabs(pstr2com[imu]), 2.);
    }
    // fprintf(stderr,"  SPmerge::next_DDiff : pstr1com = (%e, %e, %e, %e)
    // GeV\n",
    //	pstr1com[0], pstr1com[1], pstr1com[2], pstr1com[3]);
    // fprintf(stderr,"  SPmerge::next_DDiff : pstr2com = (%e, %e, %e, %e)
    // GeV\n",
    //	pstr2com[0], pstr2com[1], pstr2com[2], pstr2com[3]);

    if (mstr1 > 0.) {
      mstr1 = sqrt(mstr1);
    } else {
      mstr1 = 0.;
    }

    if (mstr2 > 0.) {
      mstr2 = sqrt(mstr2);
    } else {
      mstr2 = 0.;
    }

    foundMass1 = mstr1 > (pythia->particleData.m0(idq11) +
                          pythia->particleData.m0(idq12));
    foundMass2 = mstr2 > (pythia->particleData.m0(idq21) +
                          pythia->particleData.m0(idq22));
  }

  ret = false;
  if ((foundMass1 == true) && (foundMass2 == true)) {
    // fprintf(stderr,"  SPmerge::next_DDiff : idq11 = %d, idq12 = %d\n", idq11,
    // idq12);
    // fprintf(stderr,"  SPmerge::next_DDiff : mstr1 = %e GeV\n", mstr1);

    // fprintf(stderr,"  SPmerge::next_DDiff : idq21 = %d, idq22 = %d\n", idq21,
    // idq22);
    // fprintf(stderr,"  SPmerge::next_DDiff : mstr2 = %e GeV\n", mstr2);

    pnull = static_cast<double *>(malloc(4 * sizeof(double)));
    prs = static_cast<double *>(malloc(4 * sizeof(double)));
    evec = static_cast<double *>(malloc(4 * sizeof(double)));

    lorentz->Boost1_RestToLab(3, ucomAB, pstr1com, pstr1lab);
    lorentz->Boost1_RestToLab(3, ucomAB, pstr2com, pstr2lab);
    for (imu = 0; imu < 4; imu++) {
      ustr1com[imu] = pstr1com[imu] / mstr1;
      ustr2com[imu] = pstr2com[imu] / mstr2;

      ustr1lab[imu] = pstr1lab[imu] / mstr1;
      ustr2lab[imu] = pstr2lab[imu] / mstr2;
    }

    for (imu = 1; imu < 4; imu++) {
      pnull[imu] = pstr1com[imu];
    }
    pnull[0] = sqrt(pow(fabs(pnull[1]), 2.) + pow(fabs(pnull[2]), 2.) +
                    pow(fabs(pnull[3]), 2.));
    lorentz->Boost1_LabToRest(3, ustr1com, pnull, prs);
    pabs = sqrt(pow(fabs(prs[1]), 2.) + pow(fabs(prs[2]), 2.) +
                pow(fabs(prs[3]), 2.));
    evec[0] = 0.;
    for (imu = 1; imu < 4; imu++) {
      evec[imu] = prs[imu] / pabs;
    }
    // fprintf(stderr,"  SPmerge::next_DDiff : evec = (%e, %e, %e)\n", evec[1],
    // evec[2], evec[3]);
    nfrag1 = fragmentString(idq11, idq12, mstr1, evec, false);
    if (nfrag1 > 0) {
      NpartString1 = append_finalArray(ustr1lab, evec);
    } else {
      nfrag1 = 0;
      NpartString1 = 0;
      ret = false;
    }

    for (imu = 1; imu < 4; imu++) {
      pnull[imu] = pstr2com[imu];
    }
    pnull[0] = sqrt(pow(fabs(pnull[1]), 2.) + pow(fabs(pnull[2]), 2.) +
                    pow(fabs(pnull[3]), 2.));
    lorentz->Boost1_LabToRest(3, ustr2com, pnull, prs);
    pabs = sqrt(pow(fabs(prs[1]), 2.) + pow(fabs(prs[2]), 2.) +
                pow(fabs(prs[3]), 2.));
    evec[0] = 0.;
    for (imu = 1; imu < 4; imu++) {
      evec[imu] = prs[imu] / pabs;
    }
    // fprintf(stderr,"  SPmerge::next_DDiff : evec = (%e, %e, %e)\n", evec[1],
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

    free(pnull);
    free(prs);
    free(evec);
  }

  free(pstr1com);
  free(pstr1lab);
  free(pstr2com);
  free(pstr2lab);

  free(ustr1com);
  free(ustr1lab);
  free(ustr2com);
  free(ustr2lab);

  if (ret == true) {
    conserved = check_conservation();
    ret = ret && conserved;
  }
  return ret;
}

// non-diffractive
bool SPmerge::next_NDiff() {
  bool ret, conserved;

  int imu;
  int ntry;
  bool foundMass1, foundMass2;

  double xfracA, xfracB;
  double PPosA, PNegA;
  double PPosB, PNegB;
  double QPos, QNeg, QTrn, phiQ;
  double dPPos, dPNeg;

  int nfrag1, nfrag2;
  int idqA1, idqA2;
  int idqB1, idqB2;
  int idq11, idq12;
  int idq21, idq22;
  double mstr1, mstr2;

  double *pstr1com;
  double *pstr1lab;
  double *pstr2com;
  double *pstr2lab;

  double *ustr1com;
  double *ustr1lab;
  double *ustr2com;
  double *ustr2lab;

  double pabs;
  double *pnull;
  double *prs;
  double *evec;

  pstr1com = static_cast<double *>(malloc(4 * sizeof(double)));
  pstr1lab = static_cast<double *>(malloc(4 * sizeof(double)));
  pstr2com = static_cast<double *>(malloc(4 * sizeof(double)));
  pstr2lab = static_cast<double *>(malloc(4 * sizeof(double)));

  ustr1com = static_cast<double *>(malloc(4 * sizeof(double)));
  ustr1lab = static_cast<double *>(malloc(4 * sizeof(double)));
  ustr2com = static_cast<double *>(malloc(4 * sizeof(double)));
  ustr2lab = static_cast<double *>(malloc(4 * sizeof(double)));

  reset_finalArray();

  // PDGid2idqset(PDGidA, idqsetA);
  // PDGid2idqset(PDGidB, idqsetB);
  PPosA = (pcomA[0] + evecBasisAB[3][1] * pcomA[1] +
           evecBasisAB[3][2] * pcomA[2] + evecBasisAB[3][3] * pcomA[3]) /
          sqrt(2.);
  PNegA = (pcomA[0] - evecBasisAB[3][1] * pcomA[1] -
           evecBasisAB[3][2] * pcomA[2] - evecBasisAB[3][3] * pcomA[3]) /
          sqrt(2.);
  PPosB = (pcomB[0] + evecBasisAB[3][1] * pcomB[1] +
           evecBasisAB[3][2] * pcomB[2] + evecBasisAB[3][3] * pcomB[3]) /
          sqrt(2.);
  PNegB = (pcomB[0] - evecBasisAB[3][1] * pcomB[1] -
           evecBasisAB[3][2] * pcomB[2] - evecBasisAB[3][3] * pcomB[3]) /
          sqrt(2.);

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
              "  SPmerge::next_NDiff : incorrect baryon number of incoming "
              "hadrons.\n");
      fprintf(stderr, "  SPmerge::next_NDiff : baryonA = %d, baryonB = %d\n",
              baryonA, baryonB);
      exit(1);
    }

    xfracA = sample_XVDIS(xfracMin, alphapowV, betapowV);
    xfracB = sample_XVDIS(xfracMin, alphapowV, betapowV);
    QTrn = sample_Qperp(sigmaQperp);
    phiQ = 2. * M_PI * gsl_rng_uniform(rng);
    // fprintf(stderr,"  SPmerge::next_NDiff : xfracA = %e, xfracB = %e\n",
    // xfracA, xfracB);
    // fprintf(stderr,"  SPmerge::next_NDiff : QTrn = %e GeV, phiQ = %e\n",
    // QTrn, phiQ);

    QPos = -pow(fabs(QTrn), 2.) / (2. * xfracB * PNegB);
    QNeg = pow(fabs(QTrn), 2.) / (2. * xfracA * PPosA);

    dPPos = -xfracA * PPosA - QPos;
    dPNeg = xfracB * PNegB - QNeg;

    pstr1com[0] = (PPosA + dPPos + PNegA + dPNeg) / sqrt(2.);
    pstr2com[0] = (PPosB - dPPos + PNegB - dPNeg) / sqrt(2.);
    mstr1 = pow(fabs(pstr1com[0]), 2.);
    mstr2 = pow(fabs(pstr2com[0]), 2.);
    for (imu = 1; imu < 4; imu++) {
      pstr1com[imu] =
          evecBasisAB[3][imu] * (PPosA + dPPos - PNegA - dPNeg) / sqrt(2.) -
          evecBasisAB[1][imu] * QTrn * cos(phiQ) -
          evecBasisAB[2][imu] * QTrn * sin(phiQ);
      pstr2com[imu] =
          evecBasisAB[3][imu] * (PPosB - dPPos - PNegB + dPNeg) / sqrt(2.) +
          evecBasisAB[1][imu] * QTrn * cos(phiQ) +
          evecBasisAB[2][imu] * QTrn * sin(phiQ);
      mstr1 = mstr1 - pow(fabs(pstr1com[imu]), 2.);
      mstr2 = mstr2 - pow(fabs(pstr2com[imu]), 2.);
    }
    // fprintf(stderr,"  SPmerge::next_NDiff : pstr1com = (%e, %e, %e, %e)
    // GeV\n",
    //	pstr1com[0], pstr1com[1], pstr1com[2], pstr1com[3]);
    // fprintf(stderr,"  SPmerge::next_NDiff : pstr2com = (%e, %e, %e, %e)
    // GeV\n",
    //	pstr2com[0], pstr2com[1], pstr2com[2], pstr2com[3]);

    if (mstr1 > 0.) {
      mstr1 = sqrt(mstr1);
    } else {
      mstr1 = 0.;
    }

    if (mstr2 > 0.) {
      mstr2 = sqrt(mstr2);
    } else {
      mstr2 = 0.;
    }

    foundMass1 = mstr1 > (pythia->particleData.m0(idq11) +
                          pythia->particleData.m0(idq12));
    foundMass2 = mstr2 > (pythia->particleData.m0(idq21) +
                          pythia->particleData.m0(idq22));
  }

  ret = false;
  if ((foundMass1 == true) && (foundMass2 == true)) {
    // fprintf(stderr,"  SPmerge::next_NDiff : idq11 = %d, idq12 = %d\n", idq11,
    // idq12);
    // fprintf(stderr,"  SPmerge::next_NDiff : mstr1 = %e GeV\n", mstr1);

    // fprintf(stderr,"  SPmerge::next_NDiff : idq21 = %d, idq22 = %d\n", idq21,
    // idq22);
    // fprintf(stderr,"  SPmerge::next_NDiff : mstr2 = %e GeV\n", mstr2);

    pnull = static_cast<double *>(malloc(4 * sizeof(double)));
    prs = static_cast<double *>(malloc(4 * sizeof(double)));
    evec = static_cast<double *>(malloc(4 * sizeof(double)));

    lorentz->Boost1_RestToLab(3, ucomAB, pstr1com, pstr1lab);
    lorentz->Boost1_RestToLab(3, ucomAB, pstr2com, pstr2lab);
    for (imu = 0; imu < 4; imu++) {
      ustr1com[imu] = pstr1com[imu] / mstr1;
      ustr2com[imu] = pstr2com[imu] / mstr2;

      ustr1lab[imu] = pstr1lab[imu] / mstr1;
      ustr2lab[imu] = pstr2lab[imu] / mstr2;
    }

    for (imu = 1; imu < 4; imu++) {
      pnull[imu] = pstr1com[imu];
    }
    pnull[0] = sqrt(pow(fabs(pnull[1]), 2.) + pow(fabs(pnull[2]), 2.) +
                    pow(fabs(pnull[3]), 2.));
    lorentz->Boost1_LabToRest(3, ustr1com, pnull, prs);
    pabs = sqrt(pow(fabs(prs[1]), 2.) + pow(fabs(prs[2]), 2.) +
                pow(fabs(prs[3]), 2.));
    evec[0] = 0.;
    for (imu = 1; imu < 4; imu++) {
      evec[imu] = prs[imu] / pabs;
    }
    // fprintf(stderr,"  SPmerge::next_NDiff : evec = (%e, %e, %e)\n", evec[1],
    // evec[2], evec[3]);
    nfrag1 = fragmentString(idq11, idq12, mstr1, evec, false);
    if (nfrag1 > 0) {
      NpartString1 = append_finalArray(ustr1lab, evec);
    } else {
      nfrag1 = 0;
      NpartString1 = 0;
      ret = false;
    }

    for (imu = 1; imu < 4; imu++) {
      pnull[imu] = pstr2com[imu];
    }
    pnull[0] = sqrt(pow(fabs(pnull[1]), 2.) + pow(fabs(pnull[2]), 2.) +
                    pow(fabs(pnull[3]), 2.));
    lorentz->Boost1_LabToRest(3, ustr2com, pnull, prs);
    pabs = sqrt(pow(fabs(prs[1]), 2.) + pow(fabs(prs[2]), 2.) +
                pow(fabs(prs[3]), 2.));
    evec[0] = 0.;
    for (imu = 1; imu < 4; imu++) {
      evec[imu] = prs[imu] / pabs;
    }
    // fprintf(stderr,"  SPmerge::next_NDiff : evec = (%e, %e, %e)\n", evec[1],
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

    free(pnull);
    free(prs);
    free(evec);
  }

  free(pstr1com);
  free(pstr1lab);
  free(pstr2com);
  free(pstr2lab);

  free(ustr1com);
  free(ustr1lab);
  free(ustr2com);
  free(ustr2lab);

  if (ret == true) {
    conserved = check_conservation();
    ret = ret && conserved;
  }
  return ret;
}

/*
// baryon-antibaryon annihilation
bool SPmerge::next_BBarAnn(){
        bool ret, conserved;
        bool ranrot = false;

        int idq11, idq12;
        int idq21, idq22;
        double Mstring1, Mstring2;
        double mstringMin;

        reset_finalArray();

        // if it is baryon-antibaryon pair
        if( (baryonA == 3) && (baryonB == -3) ){ // if it is

                ret = true;

                //PDGid2idqset(PDGidA, idqsetA);
                //PDGid2idqset(PDGidB, idqsetB);

                mstringMin = pythia->particleData.m0(idq11) +
pythia->particleData.m0(idq12);
                Mstring1 = 0.5*sqrtsAB;

                if(Mstring1 < mstringMin){
                        fprintf(stderr,"  SPmerge::next_BBarAnn failure :
Mstring1 < mstringMin.\n");
                        ret = false;
                }

                mstringMin = pythia->particleData.m0(idq21) +
pythia->particleData.m0(idq22);
                Mstring2 = 0.5*sqrtsAB;

                if(Mstring2 < mstringMin){
                        fprintf(stderr,"  SPmerge::next_BBarAnn failure :
Mstring2 < mstringMin.\n");
                        ret = false;
                }

                if(ret == true){
                        NpartString1 = fragmentString(idq11, idq12, Mstring1,
NULL, ranrot);
                        NpartString2 = fragmentString(idq21, idq22, Mstring2,
NULL, ranrot);
                }

        }
        else{ // if it is not
                fprintf(stderr,"  SPmerge::next_BBarAnn failure : it is not
BBbar pair.\n");
                ret = false;
        }
        // endif baryon-antibaryon pair

        if( ret == true ){
                conserved = check_conservation();
                ret = ret && conserved;
        }
        return ret;
}
*/

void SPmerge::PDGid2idqset(int pdgid, int *idqset) {
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

void SPmerge::makeStringEnds(int *idqset, int *idq1, int *idq2) {
  int ir, ic, jc;
  int idq1tmp, idq2tmp;
  int *idqtmp;

  double rspin;
  double rflip;

  idqtmp = static_cast<int *>(malloc(3 * sizeof(int)));

  // if it is meson/baryon
  if (idqset[3] == 0) {  // meson
    ir = static_cast<int>(floor(1. + 2. * gsl_rng_uniform(rng)));

    idq1tmp = idqset[ir];
    jc = 1 + ir % 2;
    idq2tmp = idqset[jc];
  } else {  // baryon
    ir = static_cast<int>(floor(1. + 3. * gsl_rng_uniform(rng)));
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

      rspin = gsl_rng_uniform(rng);
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
    rflip = gsl_rng_uniform(rng);
    if (rflip > 0.5) {
      *idq1 = 2;
      *idq2 = -2;
    }
  }

  // fprintf(stderr,"  SPmerge::makeStringEnds : %d/%d/%d/%d decomposed into %d
  // and %d\n",
  //	idqset[3], idqset[2], idqset[1], idqset[0], *idq1, *idq2);

  free(idqtmp);
}

int SPmerge::fragmentString(int idq1, int idq2, double mString,
                            double *evecLong, bool ranrot) {
  int ret;
  bool Hforced;

  int ipart;
  int status;
  int col, acol;
  double pCM;
  double m1, m2;
  double phi, theta;
  double rswap;
  Vec4 pquark;

  double *p3vec;
  double *pvRS;

  p3vec = static_cast<double *>(malloc(4 * sizeof(double)));
  pvRS = static_cast<double *>(malloc(4 * sizeof(double)));

  pythia->event.reset();
  ret = 0;

  m1 = pythia->particleData.m0(idq1);
  m2 = pythia->particleData.m0(idq2);
  pCM = sqrt((pow(fabs(mString), 2.) - pow(fabs(m1 + m2), 2.)) *
             (pow(fabs(mString), 2.) - pow(fabs(m1 - m2), 2.))) /
        (2. * mString);

  p3vec[0] = 0.;
  if (ranrot == true) {
    theta = acos(1. - 2. * gsl_rng_uniform(rng));
    phi = 2. * M_PI * gsl_rng_uniform(rng);
    rswap = 0.;

    p3vec[1] = pCM * sin(theta) * cos(phi);
    p3vec[2] = pCM * sin(theta) * sin(phi);
    p3vec[3] = pCM * cos(theta);
  } else {
    rswap = gsl_rng_uniform(rng);

    if (evecLong == NULL) {
      p3vec[1] = 0.;
      p3vec[2] = 0.;
      p3vec[3] = pCM;
    } else {
      if (rswap < 0.5) {
        evecLong[1] = -evecLong[1];
        evecLong[2] = -evecLong[2];
        evecLong[3] = -evecLong[3];
      }
    }
  }

  if (m1 + m2 < mString) {
    status = 1;
    col = 1;
    acol = 0;
    if (ranrot == true) {
      pvRS[1] = -p3vec[1];
      pvRS[2] = -p3vec[2];
      pvRS[3] = -p3vec[3];
    } else {
      if (evecLong == NULL) {
        pvRS[1] = -p3vec[1];
        pvRS[2] = -p3vec[2];
        pvRS[3] = -p3vec[3];
      } else {
        pvRS[1] = -pCM * evecLong[1];
        pvRS[2] = -pCM * evecLong[2];
        pvRS[3] = -pCM * evecLong[3];
      }
    }
    pvRS[0] = sqrt(pow(fabs(m1), 2.) + pow(fabs(pCM), 2.));

    pquark.e(pvRS[0]);
    pquark.px(pvRS[1]);
    pquark.py(pvRS[2]);
    pquark.pz(pvRS[3]);

    pythia->event.append(idq1, status, col, acol, pquark, m1);

    status = 1;
    col = 0;
    acol = 1;
    if (ranrot == true) {
      pvRS[1] = p3vec[1];
      pvRS[2] = p3vec[2];
      pvRS[3] = p3vec[3];
    } else {
      if (evecLong == NULL) {
        pvRS[1] = p3vec[1];
        pvRS[2] = p3vec[2];
        pvRS[3] = p3vec[3];
      } else {
        pvRS[1] = pCM * evecLong[1];
        pvRS[2] = pCM * evecLong[2];
        pvRS[3] = pCM * evecLong[3];
      }
    }
    pvRS[0] = sqrt(pow(fabs(m2), 2.) + pow(fabs(pCM), 2.));

    pquark.e(pvRS[0]);
    pquark.px(pvRS[1]);
    pquark.py(pvRS[2]);
    pquark.pz(pvRS[3]);

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
    fprintf(stderr, "  SPmerge::fragmentString failure : m1 + m2 >= mString\n");
  }

  free(p3vec);
  free(pvRS);

  return ret;
}

double SPmerge::sample_XSDIS(double xmin, double beta) {
  bool accepted;

  double xfrac = 0.;
  double rx, ra;
  double pdf, env;

  accepted = false;
  while (accepted == false) {
    rx = log(1. / xmin) * gsl_rng_uniform(rng);
    xfrac = xmin * exp(rx);

    env = 1. / xfrac;
    pdf = pow(fabs(1. - xfrac), 1. + beta) / xfrac;

    ra = env * gsl_rng_uniform(rng);
    if (ra < pdf) {
      accepted = true;
    }
  }

  return xfrac;
}

double SPmerge::sample_XVDIS(double xmin, double alpha, double beta) {
  bool accepted;

  double xfrac = 0.;
  double rx, ra;
  double pdf, env;

  accepted = false;
  while (accepted == false) {
    rx = (1. - pow(fabs(xmin), alpha)) * gsl_rng_uniform(rng);
    xfrac = pow(fabs(rx + pow(fabs(xmin), alpha)), 1. / alpha);

    env = pow(fabs(xfrac), alpha - 1.) * pow(fabs(1. - xmin), beta - 1.);
    pdf = pow(fabs(xfrac), alpha - 1.) * pow(fabs(1. - xfrac), beta - 1.);

    ra = env * gsl_rng_uniform(rng);
    if (ra < pdf) {
      accepted = true;
    }
  }

  return xfrac;
}

double SPmerge::sample_Qperp(double sigQ) {
  double Qperp = 0.;
  double rq;

  rq = gsl_rng_uniform(rng);
  Qperp = sigQ * sqrt(fabs(log(1. / (1. - rq))));

  return Qperp;
}
