// GammaKinematics.cc is a part of the PYTHIA event generator.
// Copyright (C) 2017 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the GammaKinematics
// class.

#include "Pythia8/GammaKinematics.h"

namespace Pythia8 {

//==========================================================================

// The GammaKinematics class.
// Generates the kinematics of emitted photons according to phase space limits.

//--------------------------------------------------------------------------

// Initialize phase space limits.

bool GammaKinematics::init(Info* infoPtrIn, Settings* settingsPtrIn,
  Rndm* rndmPtrIn, BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn){

  // Store input pointers for future use.
  infoPtr       = infoPtrIn;
  settingsPtr   = settingsPtrIn;
  rndmPtr       = rndmPtrIn;
  beamAPtr      = beamAPtrIn;
  beamBPtr      = beamBPtrIn;

  // Save the applied cuts.
  Q2maxGamma    = settingsPtr->parm("Photon:Q2max");
  Wmin          = settingsPtr->parm("Photon:Wmin");
  Wmax          = settingsPtr->parm("Photon:Wmax");
  theta1Max     = settingsPtr->parm("Photon:thetaAMax");
  theta2Max     = settingsPtr->parm("Photon:thetaBMax");

  // Direct or resolved photons.
  gammaMode     = settingsPtr->mode("Photon:ProcessType");

  // Get the masses and collision energy and derive useful ratios.
  eCM           = infoPtr->eCM();
  sCM           = pow2( eCM);
  m2BeamA       = pow2( beamAPtr->m() );
  m2BeamB       = pow2( beamBPtr->m() );
  m2sA          = 4. * m2BeamA / sCM;
  m2sB          = 4. * m2BeamB / sCM;
  sHatNew       = 0.;

  // If Wmax below Wmin (negative by default) use the total invariant mass.
  if ( Wmax < Wmin ) Wmax = eCM;

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Sample the Q2 values and phi angles for each beam and derive kT according
// to sampled x_gamma. Check that sampled values within required limits.

bool GammaKinematics::sampleKTgamma(){

  // Get the x_gamma values set in PartonDistributions for hard processes
  // and in PhaseSpace for soft processes.
  xGamma1 = beamAPtr->xGamma();
  xGamma2 = beamBPtr->xGamma();

  // Calculate Q2 limit for given x_gamma.
  Q2min1 = 2. * m2BeamA * pow2(xGamma1) / ( 1. - xGamma1 - m2sA
         + sqrt(1. - m2sA) * sqrt( pow2(1. - xGamma1) - m2sA ) );
  Q2min2 = 2. * m2BeamB * pow2(xGamma2) / ( 1. - xGamma2 - m2sB
         + sqrt(1. - m2sB) * sqrt( pow2(1. - xGamma2) - m2sB ) );

  // Check if allowed x_gamma. May fail for direct processes.
  double xGamMaxA = Q2maxGamma / (2. * m2BeamA) * ( sqrt(
    (1. + 4. * m2BeamA / Q2maxGamma) * (1. - 4. * m2BeamA / sCM) ) - 1. );
  double xGamMaxB = Q2maxGamma / (2. * m2BeamB) * ( sqrt(
    (1. + 4. * m2BeamB / Q2maxGamma) * (1. - 4. * m2BeamB / sCM) ) - 1. );
  if ( xGamma1 > xGamMaxA || xGamma2 > xGamMaxB ) return false;

  // Sample Q2_gamma values for each beam.
  Q2gamma1 = Q2min1 * pow( Q2maxGamma / Q2min1, rndmPtr->flat() );
  Q2gamma2 = Q2min2 * pow( Q2maxGamma / Q2min2, rndmPtr->flat() );

  // Sample azimuthal angles from flat [0,2*pi[.
  phi1     = 2. * M_PI * rndmPtr->flat();
  phi2     = 2. * M_PI * rndmPtr->flat();
  double cosPhi12 = cos(phi1 - phi2);

  // Calculate the CM-energy of incoming leptons.
  eCM2A = 0.25 * pow2( sCM + m2BeamA - m2BeamB ) / sCM;
  eCM2B = 0.25 * pow2( sCM - m2BeamA + m2BeamB ) / sCM;

  // Calculate kT^2 for photons from massive leptons.
  double kT2gamma1 = ( ( 1. - xGamma1 - 0.25 * Q2gamma1 / eCM2A ) * Q2gamma1
    - m2BeamA * ( Q2gamma1 / eCM2A + pow2(xGamma1) ) ) / (1.- m2BeamA / eCM2A);
  double kT2gamma2 = ( ( 1. - xGamma2 - 0.25 * Q2gamma2 / eCM2B ) * Q2gamma2
    - m2BeamB * ( Q2gamma2 / eCM2B + pow2(xGamma2) ) ) / (1.- m2BeamB / eCM2B);

  // Check that physical values for kT's (very rarely fails if ever but may
  // cause numerical issues).
  if ( kT2gamma1 < 0. || kT2gamma2 < 0. ) {
    infoPtr->errorMsg("Error in gammaKinematics::sampleKTgamma: "
                      "unphysical kT value.");
    return false;
  }

  // Calculate the kT's.
  kT1 = sqrt( kT2gamma1 );
  kT2 = sqrt( kT2gamma2 );

  // Calculate the lepton scattering angles.
  theta1 = atan( sqrt( eCM2A* ( Q2gamma1 * ( 1. - xGamma1 )
    - m2BeamA * pow2(xGamma1) ) - m2BeamA * Q2gamma1 - pow2( 0.5 * Q2gamma1) )
    / ( eCM2A * ( 1. - xGamma1) - m2BeamA - 0.5 * Q2gamma1 ) );
  theta2 = atan( sqrt( eCM2B* ( Q2gamma2 * ( 1. - xGamma2 )
    - m2BeamB * pow2(xGamma2) ) - m2BeamB * Q2gamma2 - pow2( 0.5 * Q2gamma2) )
    / ( eCM2B * ( 1. - xGamma2) - m2BeamB - 0.5 * Q2gamma2 ) );

  // Reject kinematics if a scattering angle above cut.
  if ( theta1Max > 0 && ( theta1 > theta1Max ) ) return false;
  if ( theta2Max > 0 && ( theta2 > theta2Max ) ) return false;

  // Calculate the k_z for virtual photons including masses.
  double kz1 = (xGamma1 * eCM2A + 0.5 * Q2gamma1) / ( sqrt(eCM2A - m2BeamA) );
  double kz2 = (xGamma2 * eCM2B + 0.5 * Q2gamma2) / ( sqrt(eCM2B - m2BeamB) );

  // Calculate invariant mass for gamma-gamma pair with kT.
  m2GmGm = 2. * sqrt( eCM2A * eCM2B) * xGamma1 * xGamma2 - Q2gamma1 - Q2gamma2
         + 2. * kz1 * kz2 - 2. * kT1 * kT2 * cosPhi12;

  // Check if derived value within bounds set by user.
  if ( ( m2GmGm < pow2(Wmin) ) || ( m2GmGm > pow2(Wmax) ) ) return false;

  // Calculate invariant mass now that the square is positive.
  mGmGm = sqrt(m2GmGm);

  return true;
}

//--------------------------------------------------------------------------

// Calculates the new sHat for direct-direct and direct-resolved processes.

double GammaKinematics::calcNewSHat(double sHatOld){
  if      (gammaMode == 4) sHatNew = m2GmGm;
  else if (gammaMode == 2 || gammaMode == 3)
    sHatNew = sHatOld * m2GmGm / ( xGamma1 * xGamma2 * sCM);
  else sHatNew = 0.;
  return sHatNew;
}

//--------------------------------------------------------------------------

// Save the accepted values for further use.

bool GammaKinematics::finalize(){

  // Propagate the sampled values for beam particles.
  beamAPtr->newGammaKTPhi(kT1, phi1);
  beamBPtr->newGammaKTPhi(kT2, phi2);
  beamAPtr->Q2Gamma(Q2gamma1);
  beamBPtr->Q2Gamma(Q2gamma2);

  // Set the sampled values also to Info object.
  infoPtr->setQ2Gamma1(Q2gamma1);
  infoPtr->setQ2Gamma2(Q2gamma2);
  infoPtr->setX1Gamma(xGamma1);
  infoPtr->setX2Gamma(xGamma2);
  infoPtr->setTheta1(theta1);
  infoPtr->setTheta2(theta2);
  infoPtr->setECMsub(mGmGm);
  infoPtr->setsHatNew(sHatNew);

  // Done.
  return true;
}

//==========================================================================

} // end namespace Pythia8
