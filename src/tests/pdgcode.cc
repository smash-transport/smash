/*
 *
 *    Copyright (c) 2014-2023,2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/pdgcode.h"

#include "setup.h"
#include "smash/pdgcode_constants.h"

using namespace smash;

// mesons:
PdgCode pinull(0x111);
PdgCode pion(0x211);
PdgCode piminus(-0x211);
PdgCode eta(0x221);
PdgCode rhominus(-0x213);
PdgCode omega(0x223);
PdgCode phi(0x333);
PdgCode K0(0x311);
PdgCode Kplus(0x321);
PdgCode Kminus(-0x321);
PdgCode K0L(0x130);
PdgCode K0S(0x310);
PdgCode dminus(-0x411);
PdgCode bnulls(0x531);
PdgCode bPcbar(-0x541);
PdgCode eta_pr(0x331);
PdgCode j_psi(0x443);
// baryons:
PdgCode neutron(0x2112);
PdgCode proton(0x2212);
PdgCode antiproton(-0x2212);
PdgCode delta(0x2224);
PdgCode nstar(0x202112);       // N*(1440)^+
PdgCode antidelta(-0x122224);  // Δ(1700)
PdgCode sigmanull(0x3212);
PdgCode sigma(0x3222);
PdgCode lambda(0x3122);
PdgCode antixi(-0x103312);  // Anti-Ξ(1820)
PdgCode omega_bar(-0x3334);
PdgCode lambda_c(0x4122);
PdgCode sigma_c_bar(-0x4114);
PdgCode xi_c(0x4322);
PdgCode omega_c_bar(-0x4332);
PdgCode xi_cc_bar(-0x4422);
PdgCode omega_bc(0x5432);
PdgCode omega_bb(0x5532);
PdgCode deuteron("1000010020");
PdgCode antideutron("-1000010020");
PdgCode triton("1000010030");
PdgCode antitriton("-1000010030");
PdgCode he3("1000020030");
PdgCode antihe3("-1000020030");
PdgCode He4("1000020040");
PdgCode H3L("1010010030");
PdgCode antiH3L("-1010010030");
// non-hadrons:
// leptons
PdgCode electron(0x11);
PdgCode antimu(-0x13);
// bosons
PdgCode photon(0x22);
PdgCode higgs(0x25);

TEST(create_nuclei_using_hexadecimal_codes) {
  VERIFY(deuteron == PdgCode{0x1000010020});
  VERIFY(antideutron == PdgCode{-0x1000010020});
  VERIFY(triton == PdgCode{0x1000010030});
  VERIFY(antitriton == PdgCode{-0x1000010030});
  VERIFY(he3 == PdgCode{0x1000020030});
  VERIFY(antihe3 == PdgCode{-0x1000020030});
  VERIFY(He4 == PdgCode{0x1000020040});
  VERIFY(H3L == PdgCode{0x1010010030});
  VERIFY(antiH3L == PdgCode{-0x1010010030});
}

TEST(write_codes) {
  std::printf("######################### Non-Hadrons:\n");
  std::printf("e^-:       %8s %8x 0x%08x\n", electron.string().c_str(),
              electron.code(), electron.dump());
  std::printf("μ^+:       %8s %8x 0x%08x\n", antimu.string().c_str(),
              antimu.code(), antimu.dump());
  std::printf("γ:         %8s %8x 0x%08x\n", photon.string().c_str(),
              photon.code(), photon.dump());
  std::printf("############################## Mesons:\n");
  std::printf("π^+:       %8s %8x 0x%08x\n", pion.string().c_str(), pion.code(),
              pion.dump());
  std::printf("K^0:       %8s %8x 0x%08x\n", K0.string().c_str(), K0.code(),
              K0.dump());
  std::printf("K^0L:      %8s %8x 0x%08x\n", K0L.string().c_str(), K0L.code(),
              K0L.dump());
  std::printf("K^0S:      %8s %8x 0x%08x\n", K0S.string().c_str(), K0S.code(),
              K0S.dump());
  std::printf("K^-:       %8s %8x 0x%08x\n", Kminus.string().c_str(),
              Kminus.code(), Kminus.dump());
  std::printf("D^-:       %8s %8x 0x%08x\n", dminus.string().c_str(),
              dminus.code(), dminus.dump());
  std::printf("B^0_s:     %8s %8x 0x%08x\n", bnulls.string().c_str(),
              bnulls.code(), bnulls.dump());
  std::printf("bar B^+_c: %8s %8x 0x%08x\n", bPcbar.string().c_str(),
              bPcbar.code(), bPcbar.dump());
  std::printf("η^-:       %8s %8x 0x%08x\n", eta_pr.string().c_str(),
              eta_pr.code(), eta_pr.dump());
  std::printf("J/Ψ:       %8s %8x 0x%08x\n", j_psi.string().c_str(),
              j_psi.code(), j_psi.dump());
  std::printf("############################# Baryons:\n");
  std::printf("p:         %8s %8x 0x%08x\n", proton.string().c_str(),
              proton.code(), proton.dump());
  std::printf("bar Δ(1700)%8s %8x 0x%08x\n", antidelta.string().c_str(),
              antidelta.code(), antidelta.dump());
  std::printf("Σ:         %8s %8x 0x%08x\n", sigma.string().c_str(),
              sigma.code(), sigma.dump());
  std::printf("Λ:         %8s %8x 0x%08x\n", lambda.string().c_str(),
              lambda.code(), lambda.dump());
  std::printf("bar Ξ(1820)%8s %8x 0x%08x\n", antixi.string().c_str(),
              antixi.code(), antixi.dump());
  std::printf("bar Ω:     %8s %8x 0x%08x\n", omega_bar.string().c_str(),
              omega_bar.code(), omega_bar.dump());
  std::printf("Λ_c:       %8s %8x 0x%08x\n", lambda_c.string().c_str(),
              lambda_c.code(), lambda_c.dump());
  std::printf("bar Σ_c:   %8s %8x 0x%08x\n", sigma_c_bar.string().c_str(),
              sigma_c_bar.code(), sigma_c_bar.dump());
  std::printf("Ξ_c:       %8s %8x 0x%08x\n", xi_c.string().c_str(), xi_c.code(),
              xi_c.dump());
  std::printf("bar Ω_c:   %8s %8x 0x%08x\n", omega_c_bar.string().c_str(),
              omega_c_bar.code(), omega_c_bar.dump());
  std::printf("bar Ξ_cc:  %8s %8x 0x%08x\n", xi_cc_bar.string().c_str(),
              xi_cc_bar.code(), xi_cc_bar.dump());
  std::printf("Ω_bc:      %8s %8x 0x%08x\n", omega_bc.string().c_str(),
              omega_bc.code(), omega_bc.dump());
}
TEST(code) {
  COMPARE(electron.code(), 0x11);
  COMPARE(antimu.code(), static_cast<int>(0xffffffedu));
  COMPARE(photon.code(), 0x22);
  COMPARE(pion.code(), 0x211);
  COMPARE(K0.code(), 0x311);
  COMPARE(K0L.code(), 0x130);
  COMPARE(K0S.code(), 0x310);
  COMPARE(proton.code(), 0x2212);
  COMPARE(antidelta.code(), static_cast<int>(0xffeddddcu));
  COMPARE(lambda.code(), 0x3122);
  COMPARE(antixi.code(), static_cast<int>(0xffefcceeu));
}
TEST(dump) {
  COMPARE(electron.dump(), 0x11u);
  COMPARE(antimu.dump(), 0x80000013u);
  COMPARE(photon.dump(), 0x22u);
  COMPARE(pion.dump(), 0x211u);
  COMPARE(K0.dump(), 0x311u);
  COMPARE(K0L.dump(), 0x130u);
  COMPARE(K0S.dump(), 0x310u);
  COMPARE(proton.dump(), 0x2212u);
  COMPARE(antidelta.dump(), 0x80122224u);
  COMPARE(lambda.dump(), 0x3122u);
  COMPARE(antixi.dump(), 0x80103312u);
}
TEST(pdg_string) {
  COMPARE(electron.string(), "11");
  COMPARE(antimu.string(), "-13");
  COMPARE(photon.string(), "22");
  COMPARE(pion.string(), "211");
  COMPARE(K0.string(), "311");
  COMPARE(K0L.string(), "130");
  COMPARE(K0S.string(), "310");
  COMPARE(proton.string(), "2212");
  COMPARE(antidelta.string(), "-122224");
  COMPARE(lambda.string(), "3122");
  COMPARE(antixi.string(), "-103312");
  COMPARE(deuteron.string(), "1000010020");
  COMPARE(antideutron.string(), "-1000010020");
  COMPARE(He4.string(), "1000020040");
  COMPARE(H3L.string(), "1010010030");
}
TEST(decimal) {
  COMPARE(electron.get_decimal(), 11);
  COMPARE(antimu.get_decimal(), -13);
  COMPARE(photon.get_decimal(), 22);
  COMPARE(pion.get_decimal(), 211);
  COMPARE(K0.get_decimal(), 311);
  COMPARE(proton.get_decimal(), 2212);
  COMPARE(antidelta.get_decimal(), -122224);
  COMPARE(lambda.get_decimal(), 3122);
  COMPARE(antixi.get_decimal(), -103312);
  COMPARE(deuteron.get_decimal(), 1000010020);
  COMPARE(antideutron.get_decimal(), -1000010020);
  COMPARE(He4.get_decimal(), 1000020040);
  COMPARE(H3L.get_decimal(), 1010010030);
}
TEST(hexadecimal) {
  const PdgCode lambda_2350(0x990312a);
  COMPARE(lambda_2350.code(), 0x990312a);
  COMPARE(lambda_2350.dump(), 0x990312au);
  COMPARE(lambda_2350.string(), "19903129");
  COMPARE(lambda_2350.get_decimal(), 19903129);
}

TEST(hadron) {
  VERIFY(!electron.is_hadron());
  VERIFY(pion.is_hadron());
  VERIFY(proton.is_hadron());
  VERIFY(antidelta.is_hadron());
  VERIFY(!deuteron.is_hadron());
  VERIFY(!antideutron.is_hadron());
}

TEST(nucleus) {
  VERIFY(!electron.is_nucleus());
  VERIFY(!pion.is_nucleus());
  VERIFY(!proton.is_nucleus());
  VERIFY(!antidelta.is_nucleus());
  VERIFY(deuteron.is_nucleus());
  VERIFY(antideutron.is_nucleus());
  VERIFY(He4.is_nucleus());
  VERIFY(H3L.is_nucleus());
}

TEST(lepton) {
  VERIFY(electron.is_lepton());
  VERIFY(antimu.is_lepton());

  VERIFY(!photon.is_lepton());
  VERIFY(!pion.is_lepton());
  VERIFY(!proton.is_lepton());
  VERIFY(!deuteron.is_lepton());
}

TEST(is_meson) {
  VERIFY(pion.is_meson());
  VERIFY(!deuteron.is_meson());
  VERIFY(!antideutron.is_meson());
  VERIFY(!proton.is_meson());
  VERIFY(!antiproton.is_meson());
  VERIFY(!antimu.is_meson());
}

TEST(is_baryon) {
  VERIFY(!pion.is_baryon());
  VERIFY(!deuteron.is_baryon());
  VERIFY(!antideutron.is_baryon());
  VERIFY(proton.is_baryon());
  VERIFY(antiproton.is_baryon());
  VERIFY(!antimu.is_baryon());
}

TEST(dilepton) {
  VERIFY(is_dilepton(0x11, -0x11));
  VERIFY(is_dilepton(-0x11, 0x11));
  VERIFY(is_dilepton(0x13, -0x13));
  VERIFY(!is_dilepton(0x11, 0x11));
  VERIFY(!is_dilepton(-0x13, -0x13));
  VERIFY(!is_dilepton(0x211, -0x211));
}

TEST(baryon_number) {
  COMPARE(electron.baryon_number(), 0);
  COMPARE(antimu.baryon_number(), 0);
  COMPARE(photon.baryon_number(), 0);
  COMPARE(pion.baryon_number(), 0);
  COMPARE(K0.baryon_number(), 0);
  COMPARE(K0L.baryon_number(), 0);
  COMPARE(K0S.baryon_number(), 0);
  COMPARE(Kminus.baryon_number(), 0);
  COMPARE(dminus.baryon_number(), 0);
  COMPARE(bnulls.baryon_number(), 0);
  COMPARE(bPcbar.baryon_number(), 0);
  COMPARE(eta_pr.baryon_number(), 0);
  COMPARE(j_psi.baryon_number(), 0);
  COMPARE(proton.baryon_number(), 1);
  COMPARE(antidelta.baryon_number(), -1);
  COMPARE(sigma.baryon_number(), 1);
  COMPARE(lambda.baryon_number(), 1);
  COMPARE(antixi.baryon_number(), -1);
  COMPARE(omega_bar.baryon_number(), -1);
  COMPARE(lambda_c.baryon_number(), 1);
  COMPARE(sigma_c_bar.baryon_number(), -1);
  COMPARE(xi_c.baryon_number(), 1);
  COMPARE(omega_c_bar.baryon_number(), -1);
  COMPARE(xi_cc_bar.baryon_number(), -1);
  COMPARE(omega_bc.baryon_number(), 1);
  COMPARE(deuteron.baryon_number(), 2);
  COMPARE(antideutron.baryon_number(), -2);
  COMPARE(antitriton.baryon_number(), -3);
  COMPARE(he3.baryon_number(), 3);
  COMPARE(He4.baryon_number(), 4);
  COMPARE(H3L.baryon_number(), 3);
}
TEST(isospin3) {
  COMPARE(electron.isospin3(), 0);
  COMPARE(antimu.isospin3(), 0);
  COMPARE(photon.isospin3(), 0);
  COMPARE(pion.isospin3(), +2);
  COMPARE(K0.isospin3(), -1);
  COMPARE(Kminus.isospin3(), -1);
  COMPARE(dminus.isospin3(), -1);
  COMPARE(bnulls.isospin3(), 0);
  COMPARE(bPcbar.isospin3(), 0);
  COMPARE(eta_pr.isospin3(), 0);
  COMPARE(j_psi.isospin3(), 0);
  COMPARE(proton.isospin3(), 1);
  COMPARE(antidelta.isospin3(), -3);
  COMPARE(sigma.isospin3(), +2);
  COMPARE(lambda.isospin3(), 0);
  COMPARE(antixi.isospin3(), +1);
  COMPARE(omega_bar.isospin3(), 0);
  COMPARE(lambda_c.isospin3(), 0);
  COMPARE(sigma_c_bar.isospin3(), +2);
  COMPARE(xi_c.isospin3(), +1);
  COMPARE(omega_c_bar.isospin3(), 0);
  COMPARE(xi_cc_bar.isospin3(), -1);
  COMPARE(omega_bc.isospin3(), 0);
  COMPARE(deuteron.isospin3(), 0);
  COMPARE(antideutron.isospin3(), 0);
  COMPARE(triton.isospin3(), -1);
  COMPARE(antitriton.isospin3(), 1);
  COMPARE(he3.isospin3(), 1);
  COMPARE(antihe3.isospin3(), -1);
  COMPARE(H3L.isospin3(), 0);
  COMPARE(antiH3L.isospin3(), 0);
}
TEST(strangeness) {
  COMPARE(electron.strangeness(), 0);
  COMPARE(antimu.strangeness(), 0);
  COMPARE(photon.strangeness(), 0);
  COMPARE(pion.strangeness(), 0);
  COMPARE(K0.strangeness(), 1);
  COMPARE(Kminus.strangeness(), -1);
  COMPARE(dminus.strangeness(), 0);
  COMPARE(bnulls.strangeness(), -1);
  COMPARE(bPcbar.strangeness(), 0);
  COMPARE(eta_pr.strangeness(), 0);
  COMPARE(j_psi.strangeness(), 0);
  COMPARE(proton.strangeness(), 0);
  COMPARE(antidelta.strangeness(), 0);
  COMPARE(sigma.strangeness(), -1);
  COMPARE(lambda.strangeness(), -1);
  COMPARE(antixi.strangeness(), +2);
  COMPARE(omega_bar.strangeness(), +3);
  COMPARE(lambda_c.strangeness(), 0);
  COMPARE(sigma_c_bar.strangeness(), 0);
  COMPARE(xi_c.strangeness(), -1);
  COMPARE(omega_c_bar.strangeness(), +2);
  COMPARE(xi_cc_bar.strangeness(), 0);
  COMPARE(omega_bc.strangeness(), -1);
  COMPARE(omega_bb.strangeness(), -1);
  COMPARE(deuteron.strangeness(), 0);
  COMPARE(antideutron.strangeness(), 0);
  COMPARE(H3L.strangeness(), -1);
}
TEST(charmness) {
  COMPARE(electron.charmness(), 0);
  COMPARE(antimu.charmness(), 0);
  COMPARE(photon.charmness(), 0);
  COMPARE(pion.charmness(), 0);
  COMPARE(K0.charmness(), 0);
  COMPARE(Kminus.charmness(), 0);
  COMPARE(dminus.charmness(), -1);
  COMPARE(bnulls.charmness(), 0);
  COMPARE(bPcbar.charmness(), -1);
  COMPARE(eta_pr.charmness(), 0);
  COMPARE(j_psi.charmness(), 0);
  COMPARE(proton.charmness(), 0);
  COMPARE(antidelta.charmness(), 0);
  COMPARE(sigma.charmness(), 0);
  COMPARE(lambda.charmness(), 0);
  COMPARE(antixi.charmness(), 0);
  COMPARE(omega_bar.charmness(), 0);
  COMPARE(lambda_c.charmness(), +1);
  COMPARE(sigma_c_bar.charmness(), -1);
  COMPARE(xi_c.charmness(), +1);
  COMPARE(omega_c_bar.charmness(), -1);
  COMPARE(xi_cc_bar.charmness(), -2);
  COMPARE(omega_bc.charmness(), +1);
  COMPARE(omega_bb.charmness(), 0);
  COMPARE(deuteron.charmness(), 0);
}
TEST(bottomness) {
  COMPARE(electron.bottomness(), 0);
  COMPARE(antimu.bottomness(), 0);
  COMPARE(photon.bottomness(), 0);
  COMPARE(pion.bottomness(), 0);
  COMPARE(K0.bottomness(), 0);
  COMPARE(Kminus.bottomness(), 0);
  COMPARE(dminus.bottomness(), 0);
  COMPARE(bnulls.bottomness(), 1);
  COMPARE(bPcbar.bottomness(), -1);
  COMPARE(eta_pr.bottomness(), 0);
  COMPARE(j_psi.bottomness(), 0);
  COMPARE(proton.bottomness(), 0);
  COMPARE(antidelta.bottomness(), 0);
  COMPARE(sigma.bottomness(), 0);
  COMPARE(lambda.bottomness(), 0);
  COMPARE(antixi.bottomness(), 0);
  COMPARE(omega_bar.bottomness(), 0);
  COMPARE(lambda_c.bottomness(), 0);
  COMPARE(sigma_c_bar.bottomness(), 0);
  COMPARE(xi_c.bottomness(), 0);
  COMPARE(omega_c_bar.bottomness(), 0);
  COMPARE(xi_cc_bar.bottomness(), 0);
  COMPARE(omega_bc.bottomness(), -1);
  COMPARE(omega_bb.bottomness(), -2);
  COMPARE(deuteron.bottomness(), 0);
}
TEST(frac_strange) {
  COMPARE(electron.frac_strange(), 0);
  COMPARE(antimu.frac_strange(), 0);
  COMPARE(photon.frac_strange(), 0);
  COMPARE(pion.frac_strange(), 0);
  COMPARE(K0.frac_strange(), 1. / 2);
  COMPARE(Kminus.frac_strange(), 1. / 2);
  COMPARE(dminus.frac_strange(), 0);
  COMPARE(bnulls.frac_strange(), 1. / 2);
  COMPARE(bPcbar.frac_strange(), 0);
  COMPARE(eta_pr.frac_strange(), 1);
  COMPARE(j_psi.frac_strange(), 0);
  COMPARE(proton.frac_strange(), 0);
  COMPARE(antidelta.frac_strange(), 0);
  COMPARE(sigma.frac_strange(), 1. / 3);
  COMPARE(lambda.frac_strange(), 1. / 3);
  COMPARE(antixi.frac_strange(), 2. / 3);
  COMPARE(omega_bar.frac_strange(), 1);
  COMPARE(lambda_c.frac_strange(), 0);
  COMPARE(sigma_c_bar.frac_strange(), 0);
  COMPARE(xi_c.frac_strange(), 1. / 3);
  COMPARE(omega_c_bar.frac_strange(), 2. / 3);
  COMPARE(xi_cc_bar.frac_strange(), 0);
  COMPARE(omega_bc.frac_strange(), 1. / 3);
  COMPARE(omega_bb.frac_strange(), 1. / 3);
  COMPARE(deuteron.frac_strange(), 0);
}
TEST(frac_bottom) {
  COMPARE(electron.frac_bottom(), 0);
  COMPARE(antimu.frac_bottom(), 0);
  COMPARE(photon.frac_bottom(), 0);
  COMPARE(pion.frac_bottom(), 0);
  COMPARE(K0.frac_bottom(), 0);
  COMPARE(Kminus.frac_bottom(), 0);
  COMPARE(dminus.frac_bottom(), 0);
  COMPARE(bnulls.frac_bottom(), 1. / 2);
  COMPARE(bPcbar.frac_bottom(), 1. / 2);
  COMPARE(eta_pr.frac_bottom(), 0);
  COMPARE(j_psi.frac_bottom(), 0);
  COMPARE(proton.frac_bottom(), 0);
  COMPARE(antidelta.frac_bottom(), 0);
  COMPARE(sigma.frac_bottom(), 0);
  COMPARE(lambda.frac_bottom(), 0);
  COMPARE(antixi.frac_bottom(), 0);
  COMPARE(omega_bar.frac_bottom(), 0);
  COMPARE(lambda_c.frac_bottom(), 0);
  COMPARE(sigma_c_bar.frac_bottom(), 0);
  COMPARE(xi_c.frac_bottom(), 0);
  COMPARE(omega_c_bar.frac_bottom(), 0);
  COMPARE(xi_cc_bar.frac_bottom(), 0);
  COMPARE(omega_bc.frac_bottom(), 1. / 3);
  COMPARE(omega_bb.frac_bottom(), 2. / 3);
  COMPARE(deuteron.frac_bottom(), 0);
}
TEST(frac_charm) {
  COMPARE(electron.frac_charm(), 0);
  COMPARE(antimu.frac_charm(), 0);
  COMPARE(photon.frac_charm(), 0);
  COMPARE(pion.frac_charm(), 0);
  COMPARE(K0.frac_charm(), 0);
  COMPARE(Kminus.frac_charm(), 0);
  COMPARE(dminus.frac_charm(), 1. / 2);
  COMPARE(bnulls.frac_charm(), 0);
  COMPARE(bPcbar.frac_charm(), 1. / 2);
  COMPARE(eta_pr.frac_charm(), 0);
  COMPARE(j_psi.frac_charm(), 1);
  COMPARE(proton.frac_charm(), 0);
  COMPARE(antidelta.frac_charm(), 0);
  COMPARE(sigma.frac_charm(), 0);
  COMPARE(lambda.frac_charm(), 0);
  COMPARE(antixi.frac_charm(), 0);
  COMPARE(omega_bar.frac_charm(), 0);
  COMPARE(lambda_c.frac_charm(), 1. / 3);
  COMPARE(sigma_c_bar.frac_charm(), 1. / 3);
  COMPARE(xi_c.frac_charm(), 1. / 3);
  COMPARE(omega_c_bar.frac_charm(), 1. / 3);
  COMPARE(xi_cc_bar.frac_charm(), 2. / 3);
  COMPARE(omega_bc.frac_charm(), 1. / 3);
  COMPARE(omega_bb.frac_charm(), 0);
  COMPARE(deuteron.frac_charm(), 0);
}
TEST(heavy_flavor) {
  COMPARE(electron.is_heavy_flavor(), 0);
  COMPARE(antimu.is_heavy_flavor(), 0);
  COMPARE(photon.is_heavy_flavor(), 0);
  COMPARE(pion.is_heavy_flavor(), 0);
  COMPARE(K0.is_heavy_flavor(), 0);
  COMPARE(Kminus.is_heavy_flavor(), 0);
  COMPARE(dminus.is_heavy_flavor(), 1);
  COMPARE(bnulls.is_heavy_flavor(), 1);
  COMPARE(bPcbar.is_heavy_flavor(), 1);
  COMPARE(eta_pr.is_heavy_flavor(), 0);
  COMPARE(j_psi.is_heavy_flavor(), 1);
  COMPARE(proton.is_heavy_flavor(), 0);
  COMPARE(antidelta.is_heavy_flavor(), 0);
  COMPARE(sigma.is_heavy_flavor(), 0);
  COMPARE(lambda.is_heavy_flavor(), 0);
  COMPARE(antixi.is_heavy_flavor(), 0);
  COMPARE(omega_bar.is_heavy_flavor(), 0);
  COMPARE(lambda_c.is_heavy_flavor(), 1);
  COMPARE(sigma_c_bar.is_heavy_flavor(), 1);
  COMPARE(xi_c.is_heavy_flavor(), 1);
  COMPARE(omega_c_bar.is_heavy_flavor(), 1);
  COMPARE(xi_cc_bar.is_heavy_flavor(), 1);
  COMPARE(omega_bc.is_heavy_flavor(), 1);
  COMPARE(omega_bb.is_heavy_flavor(), 1);
  COMPARE(deuteron.is_heavy_flavor(), 0);
}
TEST(charge) {
  COMPARE(electron.charge(), -1);
  COMPARE(antimu.charge(), +1);
  COMPARE(photon.charge(), 0);
  COMPARE(pion.charge(), +1);
  COMPARE(K0.charge(), 0);
  COMPARE(K0L.charge(), 0);
  COMPARE(K0S.charge(), 0);
  COMPARE(Kminus.charge(), -1);
  COMPARE(dminus.charge(), -1);
  COMPARE(bnulls.charge(), 0);
  COMPARE(bPcbar.charge(), -1);
  COMPARE(eta_pr.charge(), 0);
  COMPARE(j_psi.charge(), 0);
  COMPARE(proton.charge(), +1);
  COMPARE(antidelta.charge(), -2);
  COMPARE(sigma.charge(), +1);
  COMPARE(lambda.charge(), 0);
  COMPARE(antixi.charge(), +1);
  COMPARE(omega_bar.charge(), +1);
  COMPARE(lambda_c.charge(), +1);
  COMPARE(sigma_c_bar.charge(), 0);
  COMPARE(xi_c.charge(), +1);
  COMPARE(omega_c_bar.charge(), 0);
  COMPARE(xi_cc_bar.charge(), -2);
  COMPARE(omega_bc.charge(), 0);
  COMPARE(deuteron.charge(), 1);
  COMPARE(antideutron.charge(), -1);
  COMPARE(He4.charge(), 2);
  COMPARE(H3L.charge(), 1);
}
TEST(quarks) {
  COMPARE(electron.quarks(), 0x000);
  COMPARE(antimu.quarks(), 0x000);
  COMPARE(photon.quarks(), 0x000);
  COMPARE(pion.quarks(), 0x021);
  COMPARE(K0.quarks(), 0x031);
  COMPARE(K0L.quarks(), 0x013);
  COMPARE(K0S.quarks(), 0x031);
  COMPARE(Kminus.quarks(), 0x032);
  COMPARE(dminus.quarks(), 0x041);
  COMPARE(bnulls.quarks(), 0x053);
  COMPARE(bPcbar.quarks(), 0x054);
  COMPARE(eta_pr.quarks(), 0x033);
  COMPARE(j_psi.quarks(), 0x044);
  COMPARE(proton.quarks(), 0x221);
  COMPARE(antidelta.quarks(), 0x222);
  COMPARE(sigma.quarks(), 0x322);
  COMPARE(lambda.quarks(), 0x312);
  COMPARE(antixi.quarks(), 0x331);
  COMPARE(omega_bar.quarks(), 0x333);
  COMPARE(lambda_c.quarks(), 0x412);
  COMPARE(sigma_c_bar.quarks(), 0x411);
  COMPARE(xi_c.quarks(), 0x432);
  COMPARE(omega_c_bar.quarks(), 0x433);
  COMPARE(xi_cc_bar.quarks(), 0x442);
  COMPARE(omega_bc.quarks(), 0x543);
}
TEST(spin) {
  COMPARE(electron.spin(), 1u);
  COMPARE(antimu.spin(), 1u);
  COMPARE(photon.spin(), 2u);
  COMPARE(pion.spin(), 0u);
  COMPARE(K0.spin(), 0u);
  COMPARE(K0L.spin(), 0u);
  COMPARE(K0S.spin(), 0u);
  COMPARE(Kminus.spin(), 0u);
  COMPARE(dminus.spin(), 0u);
  COMPARE(bnulls.spin(), 0u);
  COMPARE(bPcbar.spin(), 0u);
  COMPARE(eta_pr.spin(), 0u);
  COMPARE(j_psi.spin(), 2u);
  COMPARE(proton.spin(), 1u);
  COMPARE(antidelta.spin(), 3u);
  COMPARE(sigma.spin(), 1u);
  COMPARE(lambda.spin(), 1u);
  COMPARE(antixi.spin(), 1u);
  COMPARE(omega_bar.spin(), 3u);
  COMPARE(lambda_c.spin(), 1u);
  COMPARE(sigma_c_bar.spin(), 3u);
  COMPARE(xi_c.spin(), 1u);
  COMPARE(omega_c_bar.spin(), 1u);
  COMPARE(xi_cc_bar.spin(), 1u);
  COMPARE(omega_bc.spin(), 1u);
  COMPARE(deuteron.spin(), 2u);
  COMPARE(antideutron.spin(), 2u);
  COMPARE(He4.spin(), 0u);
  COMPARE(H3L.spin(), 1u);
}
TEST(spin_higgs) {
  vir::test::expect_failure();
  COMPARE(higgs.spin(), 0u);
}
TEST(spin_degeneracy) {
  COMPARE(electron.spin_degeneracy(), 2u);
  COMPARE(antimu.spin_degeneracy(), 2u);
  COMPARE(photon.spin_degeneracy(), 3u);
  COMPARE(pion.spin_degeneracy(), 1u);
  COMPARE(K0.spin_degeneracy(), 1u);
  COMPARE(K0L.spin_degeneracy(), 1u);
  COMPARE(K0S.spin_degeneracy(), 1u);
  COMPARE(Kminus.spin_degeneracy(), 1u);
  COMPARE(dminus.spin_degeneracy(), 1u);
  COMPARE(bnulls.spin_degeneracy(), 1u);
  COMPARE(bPcbar.spin_degeneracy(), 1u);
  COMPARE(eta_pr.spin_degeneracy(), 1u);
  COMPARE(j_psi.spin_degeneracy(), 3u);
  COMPARE(proton.spin_degeneracy(), 2u);
  COMPARE(antidelta.spin_degeneracy(), 4u);
  COMPARE(sigma.spin_degeneracy(), 2u);
  COMPARE(lambda.spin_degeneracy(), 2u);
  COMPARE(antixi.spin_degeneracy(), 2u);
  COMPARE(omega_bar.spin_degeneracy(), 4u);
  COMPARE(lambda_c.spin_degeneracy(), 2u);
  COMPARE(sigma_c_bar.spin_degeneracy(), 4u);
  COMPARE(xi_c.spin_degeneracy(), 2u);
  COMPARE(omega_c_bar.spin_degeneracy(), 2u);
  COMPARE(xi_cc_bar.spin_degeneracy(), 2u);
  COMPARE(omega_bc.spin_degeneracy(), 2u);
  COMPARE(H3L.spin_degeneracy(), 2u);
}
TEST(spin_degeneracy_higgs) {
  vir::test::expect_failure();
  COMPARE(higgs.spin_degeneracy(), 1u);
}

TEST_CATCH(set_invalid_code, PdgCode::InvalidPdgCode) {
  PdgCode invalidparticle(211);
}
TEST_CATCH(set_invalid_code_hex, PdgCode::InvalidPdgCode) {
  PdgCode invalidparticle(0xfedcba98);
}
TEST_CATCH(set_invalid_code_quark, PdgCode::InvalidPdgCode) {
  PdgCode invalidparticle(0x711);
}
TEST_CATCH(set_invalid_code_nJ0_meson, PdgCode::InvalidPdgCode) {
  PdgCode invalidparticle(0x110);
}
TEST_CATCH(set_invalid_code_nJ0_baryon, PdgCode::InvalidPdgCode) {
  PdgCode invalidparticle(0x2210);
}
TEST_CATCH(set_invalid_code_nJ_meson, PdgCode::InvalidPdgCode) {
  PdgCode invalidparticle(0x112);
}
TEST_CATCH(set_invalid_code_nJ_baryon, PdgCode::InvalidPdgCode) {
  PdgCode invalidparticle(0x2211);
}
TEST_CATCH(set_invalid_code_antimeson, PdgCode::InvalidPdgCode) {
  PdgCode invalidparticle(-0x331);
}

TEST(initialize_from_string) {
  PdgCode particle1("+1234568");
  COMPARE(particle1.dump(), 0x1234568u);
  PdgCode particle2("-211");
  COMPARE(particle2.dump(), 0x80000211u);
  PdgCode particle3("1234");
  COMPARE(particle3.dump(), 0x1234u);
  // Make sure hexadecimal is supported.
  PdgCode particle4("990312a");
  COMPARE(particle4.dump(), 0x990312au);
  COMPARE(particle4, PdgCode("990312A"));
  // Make sure the alternative encoding works.
  PdgCode particle5("19903129");
  COMPARE(particle4, particle5);
}
TEST_CATCH(empty_string, PdgCode::InvalidPdgCode) { PdgCode particle(""); }
TEST_CATCH(long_string, PdgCode::InvalidPdgCode) {
  PdgCode particle("+12345678");
}
TEST_CATCH(plus_string, PdgCode::InvalidPdgCode) { PdgCode particle("+"); }
TEST_CATCH(minus_string, PdgCode::InvalidPdgCode) { PdgCode particle("-"); }
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
TEST_CATCH(invalid_digits_quark, PdgCode::InvalidPdgCode) {
  PdgCode particle("771");
}
TEST_CATCH(invalid_digits_nJ0_meson, PdgCode::InvalidPdgCode) {
  PdgCode particle("110");
}
TEST_CATCH(invalid_digits_nJ0_baryon, PdgCode::InvalidPdgCode) {
  PdgCode particle("2210");
}
TEST_CATCH(invalid_digits_nJ_meson, PdgCode::InvalidPdgCode) {
  PdgCode particle("112");
}
TEST_CATCH(invalid_digits_nJ_baryon, PdgCode::InvalidPdgCode) {
  PdgCode particle("2211");
}
TEST_CATCH(invalid_digits_antimeson, PdgCode::InvalidPdgCode) {
  PdgCode particle("-331");
}
TEST_CATCH(invalid_nucleus_10, PdgCode::InvalidPdgCode) {
  PdgCode particle("2000010020");
}
TEST_CATCH(invalid_nucleus_digits, PdgCode::InvalidPdgCode) {
  PdgCode particle("100010020");
}

TEST(stream) {
  PdgCode particle1;
  std::istringstream sourcestream("-1234568 +1234568 1234568 +123 -213");
  sourcestream >> particle1;
  COMPARE(particle1.code(), -0x1234568);
  COMPARE(particle1.dump(), 0x81234568u);
  sourcestream >> particle1;
  COMPARE(particle1.code(), 0x1234568);
  COMPARE(particle1.dump(), 0x1234568u);
  sourcestream >> particle1;
  COMPARE(particle1.code(), 0x1234568);
  COMPARE(particle1.dump(), 0x1234568u);
  sourcestream >> particle1;
  COMPARE(particle1.dump(), 0x123u);
  sourcestream >> particle1;
  COMPARE(particle1.dump(), 0x80000213u);
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
TEST(antiparticle) { VERIFY(pion.is_antiparticle_of(piminus)); }
TEST(get_antiparticle) { COMPARE(pion.get_antiparticle(), piminus); }

TEST(from_decimal) {
  COMPARE(pion, PdgCode::from_decimal(211));
  COMPARE(Kminus, PdgCode::from_decimal(-321));
  COMPARE(antixi, PdgCode::from_decimal(-103312));
  COMPARE(deuteron, PdgCode::from_decimal(1000010020));
  COMPARE(antideutron, PdgCode::from_decimal(-1000010020));
}

TEST(decimal_from_decimal_consistency) {
  Test::create_actual_particletypes();
  for (const ParticleType& t : ParticleType::list_all()) {
    int dec = t.pdgcode().get_decimal();
    COMPARE(dec, PdgCode::from_decimal(dec).get_decimal());
  }
}

TEST(antiparticles) {
  COMPARE(pion.has_antiparticle(), true);
  COMPARE(pinull.has_antiparticle(), false);
  COMPARE(K0.has_antiparticle(), true);
  COMPARE(proton.has_antiparticle(), true);
  COMPARE(delta.has_antiparticle(), true);
  COMPARE(lambda.has_antiparticle(), true);
  COMPARE(electron.has_antiparticle(), true);
  COMPARE(antimu.has_antiparticle(), true);
  COMPARE(photon.has_antiparticle(), false);
  COMPARE(deuteron.has_antiparticle(), true);
  COMPARE(antideutron.has_antiparticle(), true);
  COMPARE(He4.has_antiparticle(), true);
  COMPARE(H3L.has_antiparticle(), true);
}

TEST(pack_int) {
  VERIFY(pack(pdg::Lambda, pdg::pi_m) != pack(pdg::Sigma_z, pdg::pi_m));

  const uint32_t x = 0xfeeddead;
  const uint32_t y = 0xbeefface;
  const auto x_signed = static_cast<int32_t>(x);
  const auto y_signed = static_cast<int32_t>(y);
  const auto xy = pack(x_signed, y_signed);
  const uint64_t expected = 0xfeeddeadbeefface;
  COMPARE(xy, expected);
}

TEST(quark_content) {
  PdgCode pip(0x211), pim(-0x211), pi0(0x111), p(0x2212), n(0x2112),
      ap(-0x2212), an(-0x2112), el(0x11);
  std::array<int, 3> q;
  q = pip.quark_content();
  COMPARE(q[0], 0);
  COMPARE(q[1], 2);
  COMPARE(q[2], -1);

  q = pim.quark_content();
  COMPARE(q[0], 0);
  COMPARE(q[1], -2);
  COMPARE(q[2], 1);

  q = pi0.quark_content();
  COMPARE(q[0], 0);
  COMPARE(q[1], 1);
  COMPARE(q[2], -1);

  q = p.quark_content();
  COMPARE(q[0], 2);
  COMPARE(q[1], 2);
  COMPARE(q[2], 1);

  q = n.quark_content();
  COMPARE(q[0], 2);
  COMPARE(q[1], 1);
  COMPARE(q[2], 1);

  q = ap.quark_content();
  COMPARE(q[0], -2);
  COMPARE(q[1], -2);
  COMPARE(q[2], -1);

  q = an.quark_content();
  COMPARE(q[0], -2);
  COMPARE(q[1], -1);
  COMPARE(q[2], -1);

  q = el.quark_content();
  COMPARE(q[0], 0);
  COMPARE(q[1], 0);
  COMPARE(q[2], 0);
}

TEST(net_quark_number) {
  // pion+ has one u and one dbar
  PdgCode pip(0x211);
  VERIFY(pip.net_quark_number(1) == -1);
  VERIFY(pip.net_quark_number(2) == 1);

  // pion- has one d and one ubar
  PdgCode pim(-0x211);
  VERIFY(pim.net_quark_number(1) == 1);
  VERIFY(pim.net_quark_number(2) == -1);

  // pion0 has one vanishing net u and d quark numbers
  PdgCode pi0(0x111);
  VERIFY(pi0.net_quark_number(1) == 0);
  VERIFY(pi0.net_quark_number(2) == 0);

  // proton has two u and one d
  PdgCode p(0x2212);
  VERIFY(p.net_quark_number(1) == 1);
  VERIFY(p.net_quark_number(2) == 2);

  // neutron has one u and two d
  PdgCode n(0x2112);
  VERIFY(n.net_quark_number(1) == 2);
  VERIFY(n.net_quark_number(2) == 1);

  // antiproton has two ubar and one dbar
  PdgCode ap(-0x2212);
  VERIFY(ap.net_quark_number(1) == -1);
  VERIFY(ap.net_quark_number(2) == -2);

  // antineutron has one ubar and two dbar
  PdgCode an(-0x2112);
  VERIFY(an.net_quark_number(1) == -2);
  VERIFY(an.net_quark_number(2) == -1);
}

TEST(deexcite) {
  std::vector<PdgCode> pdg_codes = {0x321, 0x100323, -0x9902214};
  for (auto& pdg : pdg_codes) {
    pdg.deexcite();
  }
  COMPARE(pdg_codes[0], 0x321);
  COMPARE(pdg_codes[1], 0x323);
  COMPARE(pdg_codes[2], -0x2214);
}

TEST(valence_quarks) {
  // pion0, meson, baryon number of 0
  PdgCode pi0(0x111);
  VERIFY(pi0.contains_enough_valence_quarks(1));
  VERIFY(pi0.contains_enough_valence_quarks(-1));

  // antineutron, baryon, baryon number of -1
  PdgCode an(-0x2112);
  VERIFY(an.contains_enough_valence_quarks(-2));
  VERIFY(an.contains_enough_valence_quarks(-1));

  // proton, baryon, baryon number of 1
  PdgCode p(0x2212);
  VERIFY(p.contains_enough_valence_quarks(1));
  VERIFY(p.contains_enough_valence_quarks(2));

  // pion-, meson, baryon number of 0
  PdgCode pim(-0x211);
  VERIFY(!(pim.contains_enough_valence_quarks(2)));
  VERIFY(!(pim.contains_enough_valence_quarks(-2)));
}

TEST(nucleus_components) {
  VERIFY(H3L.nucleus_p() == 1);
  VERIFY(H3L.nucleus_n() == 1);
  VERIFY(H3L.nucleus_La() == 1);
  VERIFY(H3L.nucleus_ap() == 0);
  VERIFY(H3L.nucleus_an() == 0);
  VERIFY(H3L.nucleus_aLa() == 0);
}
