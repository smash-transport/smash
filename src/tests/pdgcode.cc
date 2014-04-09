/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "tests/unittest.h"
#include "include/pdgcode.h"

using namespace Smash;

// non-hadrons:
PDGCode electron( 0x11);
PDGCode   antimu(-0x13);
PDGCode   photon( 0x22);
// mesons:
PDGCode   pion( 0x211);
PDGCode   kaon( 0x311);
PDGCode kminus(-0x321);
PDGCode dminus(-0x411);
PDGCode bnulls( 0x531);
PDGCode bPcbar(-0x541);
PDGCode eta_pr( 0x331);
PDGCode  j_psi( 0x443);
// baryons:
PDGCode      proton( 0x2212);
PDGCode   antidelta(-0x2224);
PDGCode       sigma( 0x3222);
PDGCode      lambda( 0x3122);
PDGCode      antixi(-0x3312);
PDGCode   omega_bar(-0x3334);
PDGCode    lambda_c( 0x4122);
PDGCode sigma_c_bar(-0x4114);
PDGCode        xi_c( 0x4322);
PDGCode omega_c_bar(-0x4332);
PDGCode   xi_cc_bar(-0x4422);
PDGCode    omega_bc( 0x5432);

TEST(write_codes) {
  printf("################# Non-Hadrons:\n");
  printf("e^-:       %8x 0x%08x\n",   electron.code(),   electron.dump());
  printf("μ^+:       %8x 0x%08x\n",     antimu.code(),     antimu.dump());
  printf("γ:         %8x 0x%08x\n",     photon.code(),     photon.dump());
  printf("##################### Mesons:\n");
  printf("π^+:       %8x 0x%08x\n",       pion.code(),       pion.dump());
  printf("K^0:       %8x 0x%08x\n",       kaon.code(),       kaon.dump());
  printf("K^-:       %8x 0x%08x\n",     kminus.code(),     kminus.dump());
  printf("D^-:       %8x 0x%08x\n",     dminus.code(),     dminus.dump());
  printf("B^0_s:     %8x 0x%08x\n",     bnulls.code(),     bnulls.dump());
  printf("bar B^+_c: %8x 0x%08x\n",     bPcbar.code(),     bPcbar.dump());
  printf("η^-:       %8x 0x%08x\n",     eta_pr.code(),     eta_pr.dump());
  printf("J/Ψ:       %8x 0x%08x\n",      j_psi.code(),      j_psi.dump());
  printf("##################### Baryons:\n");
  printf("p:         %8x 0x%08x\n",      proton.code(),      proton.dump());
  printf("bar Δ:     %8x 0x%08x\n",   antidelta.code(),   antidelta.dump());
  printf("Σ:         %8x 0x%08x\n",       sigma.code(),       sigma.dump());
  printf("Λ:         %8x 0x%08x\n",      lambda.code(),      lambda.dump());
  printf("bar Ξ:     %8x 0x%08x\n",      antixi.code(),      antixi.dump());
  printf("bar Ω:     %8x 0x%08x\n",   omega_bar.code(),   omega_bar.dump());
  printf("Λ_c:       %8x 0x%08x\n",    lambda_c.code(),    lambda_c.dump());
  printf("bar Σ_c:   %8x 0x%08x\n", sigma_c_bar.code(), sigma_c_bar.dump());
  printf("Ξ_c:       %8x 0x%08x\n",        xi_c.code(),        xi_c.dump());
  printf("bar Ω_c:   %8x 0x%08x\n", omega_c_bar.code(), omega_c_bar.dump());
  printf("bar Ξ_cc:  %8x 0x%08x\n",   xi_cc_bar.code(),   xi_cc_bar.dump());
  printf("Ω_bc:      %8x 0x%08x\n",    omega_bc.code(),    omega_bc.dump());
}
TEST(code) {
  COMPARE( electron.code(),  0x11);
  COMPARE(   antimu.code(),  0xffffffed);
  COMPARE(   photon.code(),  0x22);
  COMPARE(     pion.code(),  0x211);
  COMPARE(     kaon.code(),  0x311);
  COMPARE(   proton.code(),  0x2212);
  COMPARE(antidelta.code(),  0xffffdddc);
  COMPARE(   lambda.code(),  0x3122);
  COMPARE(   antixi.code(),  0xffffccee);
}
TEST(dump) {
  COMPARE( electron.dump(),  0x11);
  COMPARE(   antimu.dump(),  0x80000013u);
  COMPARE(   photon.dump(),  0x22);
  COMPARE(     pion.dump(),  0x211);
  COMPARE(     kaon.dump(),  0x311);
  COMPARE(   proton.dump(),  0x2212);
  COMPARE(antidelta.dump(),  0x80002224u);
  COMPARE(   lambda.dump(),  0x3122);
  COMPARE(   antixi.dump(),  0x80003312u);
}

TEST(hadron) {
  VERIFY(! electron.is_hadron());
  VERIFY(      pion.is_hadron());
  VERIFY(    proton.is_hadron());
  VERIFY( antidelta.is_hadron());
}
TEST(baryon_number) {
  COMPARE(   electron.baryon_number(),  0);
  COMPARE(     antimu.baryon_number(),  0);
  COMPARE(     photon.baryon_number(),  0);
  COMPARE(       pion.baryon_number(),  0);
  COMPARE(       kaon.baryon_number(),  0);
  COMPARE(     kminus.baryon_number(),  0);
  COMPARE(     dminus.baryon_number(),  0);
  COMPARE(     bnulls.baryon_number(),  0);
  COMPARE(     bPcbar.baryon_number(),  0);
  COMPARE(     eta_pr.baryon_number(),  0);
  COMPARE(      j_psi.baryon_number(),  0);
  COMPARE(     proton.baryon_number(),  1);
  COMPARE(  antidelta.baryon_number(), -1);
  COMPARE(      sigma.baryon_number(),  1);
  COMPARE(     lambda.baryon_number(),  1);
  COMPARE(     antixi.baryon_number(), -1);
  COMPARE(  omega_bar.baryon_number(), -1);
  COMPARE(   lambda_c.baryon_number(),  1);
  COMPARE(sigma_c_bar.baryon_number(), -1);
  COMPARE(       xi_c.baryon_number(),  1);
  COMPARE(omega_c_bar.baryon_number(), -1);
  COMPARE(  xi_cc_bar.baryon_number(), -1);
  COMPARE(   omega_bc.baryon_number(),  1);
}
TEST(isospin3) {
  COMPARE(   electron.isospin3(),  0);
  COMPARE(     antimu.isospin3(),  0);
  COMPARE(     photon.isospin3(),  0);
  COMPARE(       pion.isospin3(), +2);
  COMPARE(       kaon.isospin3(), -1);
  COMPARE(     kminus.isospin3(), -1);
  COMPARE(     dminus.isospin3(), -1);
  COMPARE(     bnulls.isospin3(),  0);
  COMPARE(     bPcbar.isospin3(),  0);
  COMPARE(     eta_pr.isospin3(),  0);
  COMPARE(      j_psi.isospin3(),  0);
  COMPARE(     proton.isospin3(),  1);
  COMPARE(  antidelta.isospin3(), -3);
  COMPARE(      sigma.isospin3(), +2);
  COMPARE(     lambda.isospin3(),  0);
  COMPARE(     antixi.isospin3(), +1);
  COMPARE(  omega_bar.isospin3(),  0);
  COMPARE(   lambda_c.isospin3(),  0);
  COMPARE(sigma_c_bar.isospin3(), +2);
  COMPARE(       xi_c.isospin3(), +1);
  COMPARE(omega_c_bar.isospin3(),  0);
  COMPARE(  xi_cc_bar.isospin3(), -1);
  COMPARE(   omega_bc.isospin3(),  0);
}
TEST(strangeness) {
  COMPARE(   electron.strangeness(),  0);
  COMPARE(     antimu.strangeness(),  0);
  COMPARE(     photon.strangeness(),  0);
  COMPARE(       pion.strangeness(),  0);
  COMPARE(       kaon.strangeness(),  1);
  COMPARE(     kminus.strangeness(), -1);
  COMPARE(     dminus.strangeness(),  0);
  COMPARE(     bnulls.strangeness(), -1);
  COMPARE(     bPcbar.strangeness(),  0);
  COMPARE(     eta_pr.strangeness(),  0);
  COMPARE(      j_psi.strangeness(),  0);
  COMPARE(     proton.strangeness(),  0);
  COMPARE(  antidelta.strangeness(),  0);
  COMPARE(      sigma.strangeness(), -1);
  COMPARE(     lambda.strangeness(), -1);
  COMPARE(     antixi.strangeness(), +2);
  COMPARE(  omega_bar.strangeness(), +3);
  COMPARE(   lambda_c.strangeness(),  0);
  COMPARE(sigma_c_bar.strangeness(),  0);
  COMPARE(       xi_c.strangeness(), -1);
  COMPARE(omega_c_bar.strangeness(), +2);
  COMPARE(  xi_cc_bar.strangeness(),  0);
  COMPARE(   omega_bc.strangeness(), -1);
}
TEST(charmness) {
  COMPARE(   electron.charmness(),  0);
  COMPARE(     antimu.charmness(),  0);
  COMPARE(     photon.charmness(),  0);
  COMPARE(       pion.charmness(),  0);
  COMPARE(       kaon.charmness(),  0);
  COMPARE(     kminus.charmness(),  0);
  COMPARE(     dminus.charmness(), -1);
  COMPARE(     bnulls.charmness(),  0);
  COMPARE(     bPcbar.charmness(), -1);
  COMPARE(     eta_pr.charmness(),  0);
  COMPARE(      j_psi.charmness(),  0);
  COMPARE(     proton.charmness(),  0);
  COMPARE(  antidelta.charmness(),  0);
  COMPARE(      sigma.charmness(),  0);
  COMPARE(     lambda.charmness(),  0);
  COMPARE(     antixi.charmness(),  0);
  COMPARE(  omega_bar.charmness(),  0);
  COMPARE(   lambda_c.charmness(), +1);
  COMPARE(sigma_c_bar.charmness(), -1);
  COMPARE(       xi_c.charmness(), +1);
  COMPARE(omega_c_bar.charmness(), -1);
  COMPARE(  xi_cc_bar.charmness(), -2);
  COMPARE(   omega_bc.charmness(), +1);
}
TEST(bottomness) {
  COMPARE(   electron.bottomness(),  0);
  COMPARE(     antimu.bottomness(),  0);
  COMPARE(     photon.bottomness(),  0);
  COMPARE(       pion.bottomness(),  0);
  COMPARE(       kaon.bottomness(),  0);
  COMPARE(     kminus.bottomness(),  0);
  COMPARE(     dminus.bottomness(),  0);
  COMPARE(     bnulls.bottomness(),  1);
  COMPARE(     bPcbar.bottomness(), -1);
  COMPARE(     eta_pr.bottomness(),  0);
  COMPARE(      j_psi.bottomness(),  0);
  COMPARE(     proton.bottomness(),  0);
  COMPARE(  antidelta.bottomness(),  0);
  COMPARE(      sigma.bottomness(),  0);
  COMPARE(     lambda.bottomness(),  0);
  COMPARE(     antixi.bottomness(),  0);
  COMPARE(  omega_bar.bottomness(),  0);
  COMPARE(   lambda_c.bottomness(),  0);
  COMPARE(sigma_c_bar.bottomness(),  0);
  COMPARE(       xi_c.bottomness(),  0);
  COMPARE(omega_c_bar.bottomness(),  0);
  COMPARE(  xi_cc_bar.bottomness(),  0);
  COMPARE(   omega_bc.bottomness(), -1);
}
TEST(charge) {
  COMPARE(   electron.charge(), -1);
  COMPARE(     antimu.charge(), +1);
  COMPARE(     photon.charge(),  0);
  COMPARE(       pion.charge(), +1);
  COMPARE(       kaon.charge(),  0);
  COMPARE(     kminus.charge(), -1);
  COMPARE(     dminus.charge(), -1);
  COMPARE(     bnulls.charge(),  0);
  COMPARE(     bPcbar.charge(), -1);
  COMPARE(     eta_pr.charge(),  0);
  COMPARE(      j_psi.charge(),  0);
  COMPARE(     proton.charge(), +1);
  COMPARE(  antidelta.charge(), -2);
  COMPARE(      sigma.charge(), +1);
  COMPARE(     lambda.charge(),  0);
  COMPARE(     antixi.charge(), +1);
  COMPARE(  omega_bar.charge(), +1);
  COMPARE(   lambda_c.charge(), +1);
  COMPARE(sigma_c_bar.charge(),  0);
  COMPARE(       xi_c.charge(), +1);
  COMPARE(omega_c_bar.charge(),  0);
  COMPARE(  xi_cc_bar.charge(), -2);
  COMPARE(   omega_bc.charge(),  0);
}

TEST_CATCH(set_invalid_code, PDGCode::InvalidPDGCode) {
  PDGCode validparticle(211);
}
TEST_CATCH(set_invalid_code_hex, PDGCode::InvalidPDGCode) {
  PDGCode invalidparticle(0xfedcba98);
}
