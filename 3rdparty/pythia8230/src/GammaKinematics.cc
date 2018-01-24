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

  // Initial choice for the process type and whether to use external flux.
  gammaMode     = settingsPtr->mode("Photon:ProcessType");
  externalFlux  = settingsPtr->mode("PDF:lepton2gammaSet") == 2;

  // Flag from virtuality sampling.
  sampleQ2      = settingsPtr->flag("Photon:sampleQ2");

  // Check if photons from both beams or only from one beam.
  hasGammaA = beamAPtr->isLepton();
  hasGammaB = beamBPtr->isLepton();

  // Get the masses and collision energy and derive useful ratios.
  eCM           = infoPtr->eCM();
  sCM           = pow2( eCM);
  m2BeamA       = pow2( beamAPtr->m() );
  m2BeamB       = pow2( beamBPtr->m() );
  sHatNew       = 0.;

  // Calculate the CM-energies of incoming beams.
  eCM2A = 0.25 * pow2( sCM + m2BeamA - m2BeamB ) / sCM;
  eCM2B = 0.25 * pow2( sCM - m2BeamA + m2BeamB ) / sCM;

  // Derive ratios used often.
  m2eA  = m2BeamA / eCM2A;
  m2eB  = m2BeamB / eCM2B;

  // Derive the kinematic limits.
  xGammaMax1 = 2. * ( 1. - 0.25 * Q2maxGamma / eCM2A - m2eA)
    / ( 1. + sqrt((1. + 4. * m2BeamA / Q2maxGamma) * (1. - m2eA)) );
  xGammaMax2 = 2. * ( 1. - 0.25 * Q2maxGamma / eCM2B - m2eB)
    / ( 1. + sqrt((1. + 4. * m2BeamB / Q2maxGamma) * (1. - m2eB)) );

  // No limits for xGamma if Q2-integrated flux.
  if (!sampleQ2) {
    xGammaMax1 = 1.;
    xGammaMax2 = 1.;
  }

  // If Wmax below Wmin (negative by default) use the total invariant mass.
  if ( Wmax < Wmin ) Wmax = eCM;

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Sample kinematics of one or two photon beams from the original beams.

bool GammaKinematics::sampleKTgamma(){

  // Get the sampled x_gamma values from beams.
  xGamma1 = beamAPtr->xGamma();
  xGamma2 = beamBPtr->xGamma();

  // Type of current process.
  gammaMode = infoPtr->photonMode();

  // Reject already sampled x_gamma values outside kinematic bounds.
  if ( hasGammaA && (!externalFlux || ( externalFlux
    && (gammaMode == 3 || gammaMode == 4) ) ) && (xGamma1 > xGammaMax1) )
       return false;
  if ( hasGammaB && (!externalFlux || ( externalFlux
    && (gammaMode == 2 || gammaMode == 4) ) ) && (xGamma2 > xGammaMax2) )
       return false;

  // Sample virtuality for photon A.
  if ( hasGammaA ) {

    // Sample the x_gamma value if needed and check that value is valid.
    if ( externalFlux && (gammaMode == 1 || gammaMode == 2) ) {
      xGamma1 = beamAPtr->sampleXgamma();
      if ( xGamma1 > xGammaMax1 ) return false;
    }

    // Derive the accurate lower Q2 limit and sample value.
    Q2min1 = 2. * m2BeamA * pow2(xGamma1) / ( 1. - xGamma1 - m2eA
           + sqrt(1. - m2eA) * sqrt( pow2(1. - xGamma1) - m2eA ) );

    // Sample the Q2 if requested, otherwise pick use the maximum value.
    if (sampleQ2) Q2gamma1 = beamAPtr->sampleQ2gamma(Q2min1);
    else Q2gamma1 = 0.;

    // Reject sampled values outside limits (relevant for external flux).
    if ( sampleQ2 && (Q2gamma1 < Q2min1) ) return false;
  }

  // Sample virtuality for photon B.
  if ( hasGammaB ) {

    // Sample the x_gamma value if needed and check that value is valid.
    if ( externalFlux && (gammaMode == 1 || gammaMode == 3) ) {
      xGamma2 = beamBPtr->sampleXgamma();
      if ( xGamma2 > xGammaMax2 ) return false;
    }

    // Derive the accurate lower Q2 limit and sample value.
    Q2min2 = 2. * m2BeamB * pow2(xGamma2) / ( 1. - xGamma2 - m2eB
           + sqrt(1. - m2eB) * sqrt( pow2(1. - xGamma2) - m2eB ) );

    // Sample the Q2 if requested, otherwise pick use the maximum value.
    if (sampleQ2) Q2gamma2 = beamBPtr->sampleQ2gamma(Q2min2);
    else Q2gamma2 = 0.;

    // Reject sampled values outside limits (relevant for external flux).
    if ( sampleQ2 && (Q2gamma2 < Q2min2) ) return false;
  }

  // Derive the full photon momenta from the sampled values.
  if ( hasGammaA) {
    if ( !deriveKin(xGamma1, Q2gamma1, m2BeamA, eCM2A) ) return false;
    kT1      = kT;
    kz1      = kz;
    phi1     = phi;
    theta1   = theta;

    // Reject kinematics if a scattering angle above cut.
    if ( theta1Max > 0 && ( theta1 > theta1Max ) ) return false;
  }
  if ( hasGammaB) {
    if ( !deriveKin(xGamma2, Q2gamma2, m2BeamB, eCM2B) ) return false;
    kT2      = kT;
    kz2      = kz;
    phi2     = phi;
    theta2   = theta;

    // Reject kinematics if a scattering angle above cut.
    if ( theta2Max > 0 && ( theta2 > theta2Max ) ) return false;
  }

  // Invariant mass of photon-photon system.
  if ( hasGammaA && hasGammaB) {

    // Derive the invariant mass and check the kinematic limits.
    double cosPhi12 = cos(phi1 - phi2);
    m2GmGm = 2. * sqrt(eCM2A * eCM2B) * xGamma1 * xGamma2 - Q2gamma1 - Q2gamma2
           + 2. * kz1 * kz2 - 2. * kT1 * kT2 * cosPhi12;

    // Check if derived value within bounds set by user.
    if ( ( m2GmGm < pow2(Wmin) ) || ( m2GmGm > pow2(Wmax) ) ) return false;

    // Calculate invariant mass now that the square is positive.
    mGmGm = sqrt(m2GmGm);

    return true;

  // Invariant mass of photon-hadron system.
  } else if (hasGammaA || hasGammaB) {

    // Derive the invariant mass and check the limits.
    // Solve the longitudinal momentum of beam particles in CM frame.
    double pz2 = ( pow2(sCM - m2BeamA - m2BeamB) - 4. * m2BeamA * m2BeamB )
               * 0.25 / sCM;
    double pz  = sqrtpos( pz2);

    // Pick the correct beam mass and photon kinematics.
    double m2Beam  = hasGammaA ? m2BeamB : m2BeamA;
    double xGamma  = hasGammaA ? xGamma1 : xGamma2;
    double Q2gamma = hasGammaA ? Q2gamma1 : Q2gamma2;

    // Calculate the invariant mass of the photon-hadron pair and check limits.
    m2GmGm     = m2Beam - Q2gamma + 2. * ( xGamma * sqrt(eCM2A) * sqrt(eCM2B)
               + kz * pz );
    if ( ( m2GmGm < pow2(Wmin) ) || ( m2GmGm > pow2(Wmax) ) ) return false;
    mGmGm      = sqrt(m2GmGm);

    return true;
  }

  else return false;

}

//--------------------------------------------------------------------------

// Sample the Q2 values and phi angles for each beam and derive kT according
// to sampled x_gamma. Check that sampled values within required limits.

bool GammaKinematics::deriveKin(double xGamma, double Q2gamma,
                                double m2Beam, double eCM2) {

  // Sample azimuthal angle from flat [0,2*pi[.
  phi = 2. * M_PI * rndmPtr->flat();

  // Calculate kT^2 for photon from particle with non-zero mass.
  double kT2gamma = ( ( 1. - xGamma - 0.25 * Q2gamma / eCM2 ) * Q2gamma
    - m2Beam * ( Q2gamma / eCM2 + pow2(xGamma) ) ) / (1.- m2Beam / eCM2);

  // If no virtuality sampled set transverse momentum to zero.
  if ( !sampleQ2 ) kT2gamma = 0.;

  // Check that physical values for kT (very rarely fails if ever but may
  // cause numerical issues).
  if ( kT2gamma < 0. ) {
    infoPtr->errorMsg("Error in gammaKinematics::sampleKTgamma: "
                      "unphysical kT value.");
    return false;
  }

  // Calculate the transverse and longitudinal momenta and scattering angle
  // of the beam particle.
  kT = sqrt( kT2gamma );
  theta = atan( sqrt( eCM2 * ( Q2gamma * ( 1. - xGamma )
        - m2Beam * pow2(xGamma) ) - m2Beam * Q2gamma - pow2( 0.5 * Q2gamma) )
        / ( eCM2 * ( 1. - xGamma) - m2Beam - 0.5 * Q2gamma ) );
  kz = (xGamma * eCM2 + 0.5 * Q2gamma) / ( sqrt(eCM2 - m2Beam) );

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Calculates the new sHat for direct-direct and direct-resolved processes.

double GammaKinematics::calcNewSHat(double sHatOld){

  // Need to recalculate only if two photons.
  if ( hasGammaA && hasGammaB) {

    // Calculate the new sHat for direct-resolved system.
    gammaMode = infoPtr->photonMode();
    if      (gammaMode == 4) sHatNew = m2GmGm;
    else if (gammaMode == 2 || gammaMode == 3)
      sHatNew = sHatOld * m2GmGm / ( xGamma1 * xGamma2 * sCM);
  }

  // Otherwise no need for a new value.
  else sHatNew = sHatOld;

  return sHatNew;
}

//--------------------------------------------------------------------------

// Calculate weight from oversampling with approximated flux.

double GammaKinematics::fluxWeight() {

  // Initially unit weight.
  double wt = 1.;

  // Calculate the weight according to accurate flux.
  if ( sampleQ2) {
    if (hasGammaA) wt *= beamAPtr->xfFlux(22, xGamma1, Q2gamma1) /
                         beamAPtr->xfApprox(22, xGamma1, Q2gamma1);
    if (hasGammaB) wt *= beamBPtr->xfFlux(22, xGamma2, Q2gamma2) /
                         beamBPtr->xfApprox(22, xGamma2, Q2gamma2);

  // When no sampling of virtuality use the Q2-integrated flux.
  } else {
    if (hasGammaA) wt *= beamAPtr->xfFlux(22, xGamma1, Q2gamma1) /
                         beamAPtr->xf(22, xGamma1, Q2gamma1);
    if (hasGammaB) wt *= beamBPtr->xfFlux(22, xGamma2, Q2gamma2) /
                         beamBPtr->xf(22, xGamma2, Q2gamma2);
  }

  // Done.
  return wt;
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
