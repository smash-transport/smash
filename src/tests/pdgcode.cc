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
PdgCode electron( 0x11);
PdgCode   antimu(-0x13);
PdgCode   photon( 0x22);
// mesons:
PdgCode   pion( 0x211);
PdgCode   kaon( 0x311);
PdgCode kminus(-0x321);
PdgCode dminus(-0x411);
PdgCode bnulls( 0x531);
PdgCode bPcbar(-0x541);
PdgCode eta_pr( 0x331);
PdgCode  j_psi( 0x443);
// baryons:
PdgCode      proton( 0x2212);
PdgCode   antidelta(-0x2224);
PdgCode       sigma( 0x3222);
PdgCode      lambda( 0x3122);
PdgCode      antixi(-0x3312);
PdgCode   omega_bar(-0x3334);
PdgCode    lambda_c( 0x4122);
PdgCode sigma_c_bar(-0x4114);
PdgCode        xi_c( 0x4322);
PdgCode omega_c_bar(-0x4332);
PdgCode   xi_cc_bar(-0x4422);
PdgCode    omega_bc( 0x5432);

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
TEST(size) {
  COMPARE(sizeof(PdgCode), sizeof(std::uint32_t));
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
TEST(string) {
  COMPARE( electron.string(),    "11");
  COMPARE(   antimu.string(),   "-13");
  COMPARE(   photon.string(),    "22");
  COMPARE(     pion.string(),   "211");
  COMPARE(     kaon.string(),   "311");
  COMPARE(   proton.string(),  "2212");
  COMPARE(antidelta.string(), "-2224");
  COMPARE(   lambda.string(),  "3122");
  COMPARE(   antixi.string(), "-3312");
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
TEST(isospin_total) {
  PdgCode eta(0x221);
  COMPARE(   electron.isospin_total(),  0);
  COMPARE(     antimu.isospin_total(),  0);
  COMPARE(     photon.isospin_total(),  0);
  COMPARE(       pion.isospin_total(), +2);
  COMPARE(        eta.isospin_total(),  0);
  COMPARE(       kaon.isospin_total(),  1);
  COMPARE(     kminus.isospin_total(),  1);
  COMPARE(     dminus.isospin_total(),  1);
  COMPARE(     bnulls.isospin_total(),  0);
  COMPARE(     bPcbar.isospin_total(),  0);
  COMPARE(     eta_pr.isospin_total(),  0);
  COMPARE(      j_psi.isospin_total(),  0);
  COMPARE(     proton.isospin_total(),  1);
  COMPARE(  antidelta.isospin_total(),  3);
  COMPARE(      sigma.isospin_total(), +2);
  COMPARE(     lambda.isospin_total(),  0);
  COMPARE(     antixi.isospin_total(), +1);
  COMPARE(  omega_bar.isospin_total(),  0);
  COMPARE(   lambda_c.isospin_total(),  0);
  COMPARE(sigma_c_bar.isospin_total(), +2);
  COMPARE(       xi_c.isospin_total(), +1);
  COMPARE(omega_c_bar.isospin_total(),  0);
  COMPARE(  xi_cc_bar.isospin_total(), +1);
  COMPARE(   omega_bc.isospin_total(),  0);
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
TEST(quarks) {
  COMPARE(   electron.quarks(), 0x000);
  COMPARE(     antimu.quarks(), 0x000);
  COMPARE(     photon.quarks(), 0x000);
  COMPARE(       pion.quarks(), 0x021);
  COMPARE(       kaon.quarks(), 0x031);
  COMPARE(     kminus.quarks(), 0x032);
  COMPARE(     dminus.quarks(), 0x041);
  COMPARE(     bnulls.quarks(), 0x053);
  COMPARE(     bPcbar.quarks(), 0x054);
  COMPARE(     eta_pr.quarks(), 0x033);
  COMPARE(      j_psi.quarks(), 0x044);
  COMPARE(     proton.quarks(), 0x221);
  COMPARE(  antidelta.quarks(), 0x222);
  COMPARE(      sigma.quarks(), 0x322);
  COMPARE(     lambda.quarks(), 0x312);
  COMPARE(     antixi.quarks(), 0x331);
  COMPARE(  omega_bar.quarks(), 0x333);
  COMPARE(   lambda_c.quarks(), 0x412);
  COMPARE(sigma_c_bar.quarks(), 0x411);
  COMPARE(       xi_c.quarks(), 0x432);
  COMPARE(omega_c_bar.quarks(), 0x433);
  COMPARE(  xi_cc_bar.quarks(), 0x442);
  COMPARE(   omega_bc.quarks(), 0x543);
}
TEST(multiplet) {
  COMPARE(   electron.multiplet(),  0x0);
  COMPARE(     antimu.multiplet(),  0x0);
  COMPARE(     photon.multiplet(),  0x0);
  COMPARE(       pion.multiplet(),  0x1);
  COMPARE(       kaon.multiplet(),  0x1);
  COMPARE(     kminus.multiplet(),  0x1);
  COMPARE(     dminus.multiplet(),  0x1);
  COMPARE(     bnulls.multiplet(),  0x1);
  COMPARE(     bPcbar.multiplet(),  0x1);
  COMPARE(     eta_pr.multiplet(),  0x1);
  COMPARE(      j_psi.multiplet(),  0x3);
  COMPARE(     proton.multiplet(),  0x10002);
  COMPARE(  antidelta.multiplet(), -0x10004);
  COMPARE(      sigma.multiplet(),  0x10002);
  COMPARE(     lambda.multiplet(),  0x10002);
  COMPARE(     antixi.multiplet(), -0x10002);
  COMPARE(  omega_bar.multiplet(), -0x10004);
  COMPARE(   lambda_c.multiplet(),  0x10002);
  COMPARE(sigma_c_bar.multiplet(), -0x10004);
  COMPARE(       xi_c.multiplet(),  0x10002);
  COMPARE(omega_c_bar.multiplet(), -0x10002);
  COMPARE(  xi_cc_bar.multiplet(), -0x10002);
  COMPARE(   omega_bc.multiplet(),  0x10002);
}
TEST(iso_multiplet) {
  PdgCode pinull(0x111);
  PdgCode rhominus(-0x213);
  PdgCode omega(0x223);
  COMPARE(   electron.iso_multiplet(),  0x0000);
  COMPARE(     antimu.iso_multiplet(),  0x0000);
  COMPARE(     photon.iso_multiplet(),  0x0000);
  COMPARE(       pion.iso_multiplet(),  0x0111);
  COMPARE(     pinull.iso_multiplet(),  0x0111);
  COMPARE(   rhominus.iso_multiplet(),  0x0113);
  COMPARE(      omega.iso_multiplet(),  0x0223);
  COMPARE(       kaon.iso_multiplet(),  0x0311);
  COMPARE(     kminus.iso_multiplet(), -0x0311);
  COMPARE(     dminus.iso_multiplet(), -0x0411);
  COMPARE(     bnulls.iso_multiplet(),  0x0531);
  COMPARE(     bPcbar.iso_multiplet(), -0x0541);
  COMPARE(     eta_pr.iso_multiplet(),  0x0331);
  COMPARE(      j_psi.iso_multiplet(),  0x0443);
  COMPARE(     proton.iso_multiplet(),  0x1112);
  COMPARE(  antidelta.iso_multiplet(), -0x1114);
  COMPARE(      sigma.iso_multiplet(),  0x3112);
  COMPARE(     lambda.iso_multiplet(),  0x3122);
  COMPARE(     antixi.iso_multiplet(), -0x3312);
  COMPARE(  omega_bar.iso_multiplet(), -0x3334);
  COMPARE(   lambda_c.iso_multiplet(),  0x4122);
  COMPARE(sigma_c_bar.iso_multiplet(), -0x4114);
  COMPARE(       xi_c.iso_multiplet(),  0x4312);
  COMPARE(omega_c_bar.iso_multiplet(), -0x4332);
  COMPARE(  xi_cc_bar.iso_multiplet(), -0x4412);
  COMPARE(   omega_bc.iso_multiplet(),  0x5432);
}
TEST(spin) {
  COMPARE(   electron.spin(), 1);
  COMPARE(     antimu.spin(), 1);
  COMPARE(     photon.spin(), 2);
  COMPARE(       pion.spin(), 0);
  COMPARE(       kaon.spin(), 0);
  COMPARE(     kminus.spin(), 0);
  COMPARE(     dminus.spin(), 0);
  COMPARE(     bnulls.spin(), 0);
  COMPARE(     bPcbar.spin(), 0);
  COMPARE(     eta_pr.spin(), 0);
  COMPARE(      j_psi.spin(), 2);
  COMPARE(     proton.spin(), 1);
  COMPARE(  antidelta.spin(), 3);
  COMPARE(      sigma.spin(), 1);
  COMPARE(     lambda.spin(), 1);
  COMPARE(     antixi.spin(), 1);
  COMPARE(  omega_bar.spin(), 3);
  COMPARE(   lambda_c.spin(), 1);
  COMPARE(sigma_c_bar.spin(), 3);
  COMPARE(       xi_c.spin(), 1);
  COMPARE(omega_c_bar.spin(), 1);
  COMPARE(  xi_cc_bar.spin(), 1);
  COMPARE(   omega_bc.spin(), 1);
  PdgCode higgs(0x25);
  EXPECT_ASSERT_FAILURE(assert(higgs.spin() == 0));
}
TEST(spin_degeneracy) {
  COMPARE(   electron.spin_degeneracy(), 2);
  COMPARE(     antimu.spin_degeneracy(), 2);
  COMPARE(     photon.spin_degeneracy(), 3);
  COMPARE(       pion.spin_degeneracy(), 1);
  COMPARE(       kaon.spin_degeneracy(), 1);
  COMPARE(     kminus.spin_degeneracy(), 1);
  COMPARE(     dminus.spin_degeneracy(), 1);
  COMPARE(     bnulls.spin_degeneracy(), 1);
  COMPARE(     bPcbar.spin_degeneracy(), 1);
  COMPARE(     eta_pr.spin_degeneracy(), 1);
  COMPARE(      j_psi.spin_degeneracy(), 3);
  COMPARE(     proton.spin_degeneracy(), 2);
  COMPARE(  antidelta.spin_degeneracy(), 4);
  COMPARE(      sigma.spin_degeneracy(), 2);
  COMPARE(     lambda.spin_degeneracy(), 2);
  COMPARE(     antixi.spin_degeneracy(), 2);
  COMPARE(  omega_bar.spin_degeneracy(), 4);
  COMPARE(   lambda_c.spin_degeneracy(), 2);
  COMPARE(sigma_c_bar.spin_degeneracy(), 4);
  COMPARE(       xi_c.spin_degeneracy(), 2);
  COMPARE(omega_c_bar.spin_degeneracy(), 2);
  COMPARE(  xi_cc_bar.spin_degeneracy(), 2);
  COMPARE(   omega_bc.spin_degeneracy(), 2);
  PdgCode higgs(0x25);
  EXPECT_ASSERT_FAILURE(assert(higgs.spin() == 1));
}

TEST_CATCH(set_invalid_code, PdgCode::InvalidPdgCode) {
  PdgCode invalidparticle(211);
}
TEST_CATCH(set_invalid_code_hex, PdgCode::InvalidPdgCode) {
  PdgCode invalidparticle(0xfedcba98);
}
TEST(initialize_from_string) {
  PdgCode particle1("+1234567");
  COMPARE(particle1.dump(), 0x1234567);
  PdgCode particle2("-211");
  COMPARE(particle2.dump(), 0x80000211);
  PdgCode particle3("1234");
  COMPARE(particle3.dump(), 0x1234);
}
TEST_CATCH(empty_string, PdgCode::InvalidPdgCode) {
  PdgCode particle("");
}
TEST_CATCH(long_string, PdgCode::InvalidPdgCode) {
  PdgCode particle("+12345678");
}
TEST_CATCH(plus_string, PdgCode::InvalidPdgCode) {
  PdgCode particle("+");
}
TEST_CATCH(minus_string, PdgCode::InvalidPdgCode) {
  PdgCode particle("-");
}
// this tests characters with bitmasks 0x3. (of which digits are a
// subset)
TEST_CATCH(invalid_digits_colon, PdgCode::InvalidPdgCode) {
  PdgCode particle(":");
}
TEST_CATCH(invalid_digits_semi, PdgCode::InvalidPdgCode) {
  PdgCode particle(";");
}
TEST_CATCH(invalid_digits_less, PdgCode::InvalidPdgCode) {
  PdgCode particle("<");
}
TEST_CATCH(invalid_digits_equal, PdgCode::InvalidPdgCode) {
  PdgCode particle("=");
}
TEST_CATCH(invalid_digits_greater, PdgCode::InvalidPdgCode) {
  PdgCode particle(">");
}
TEST_CATCH(invalid_digits_question, PdgCode::InvalidPdgCode) {
  PdgCode particle("?");
}
// this is for the other characters.
TEST_CATCH(invalid_characters, PdgCode::InvalidPdgCode) {
  PdgCode particle("abcdef");
}
TEST(stream) {
  PdgCode particle1;
  std::istringstream sourcestream("-1234567 +1234567 1234567 +123 -214");
  sourcestream >> particle1;
  COMPARE(particle1.code(), -0x1234567);
  COMPARE(particle1.dump(), 0x81234567);
  sourcestream >> particle1;
  COMPARE(particle1.code(), 0x1234567);
  COMPARE(particle1.dump(), 0x1234567);
  sourcestream >> particle1;
  COMPARE(particle1.code(), 0x1234567);
  COMPARE(particle1.dump(), 0x1234567);
  sourcestream >> particle1;
  COMPARE(particle1.dump(), 0x123);
  sourcestream >> particle1;
  COMPARE(particle1.dump(), 0x80000214);
}
TEST(stream_fail) {
  PdgCode particle1;
  std::istringstream sourcestream("1234567 abcdefg");
  sourcestream >> particle1;
  sourcestream >> particle1;
}
TEST(stream_fail_colon_etc) {
  PdgCode particle1;
  std::istringstream sourcestream(":;<=>?");
  sourcestream >> particle1;
}
TEST(equal) {
  VERIFY(pion != eta_pr);
  PdgCode pion2(0x211);
  VERIFY(pion == pion2);
  VERIFY(pion2 < omega_bc);
}
TEST(antiparticle) {
  PdgCode antipion(-0x211);
  VERIFY(pion.is_antiparticle_of(antipion));
}
