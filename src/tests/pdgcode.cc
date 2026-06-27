/*
 *
 *    Copyright (c) 2014-2023,2025-2026
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

// mesons
PdgCode pi0(0x111);

PdgCode pi_plus(0x211);
PdgCode pi_minus(-0x211);

PdgCode eta(0x221);
PdgCode eta_prim(0x331);

PdgCode rho0(0x113);
PdgCode rho_minus(-0x213);

PdgCode omega(0x223);
PdgCode phi(0x333);

PdgCode k0(0x311);
PdgCode k0_bar(-0x311);

PdgCode k0_short(0x310);
PdgCode k0_long(0x130);

PdgCode k_plus(0x321);
PdgCode k_minus(-0x321);

PdgCode d_plus(0x411);
PdgCode d_minus(-0x411);

PdgCode d0(0x421);
PdgCode d0_bar(-0x421);

PdgCode b0(0x511);
PdgCode b0_bar(-0x511);

PdgCode bs0(0x531);
PdgCode bs0_bar(-0x531);

PdgCode bc_minus(-0x541);

PdgCode j_psi(0x443);
PdgCode upsilon(0x553);

// baryons
PdgCode neutron(0x2112);
PdgCode neutron_bar(-0x2112);

PdgCode proton(0x2212);
PdgCode proton_bar(-0x2212);

PdgCode delta_pp(0x2224);
PdgCode delta_minus(0x1114);
PdgCode delta_bar(-0x122224);  // Δ(1700)

PdgCode n_star(0x202112);  // N*(1440)^+

PdgCode sigma(0x3212);
PdgCode sigma_plus(0x3222);
PdgCode sigma_minus(0x3112);

PdgCode lambda(0x3122);
PdgCode lambda_bar(-0x3122);

PdgCode xi_bar(-0x103312);  // Anti-Ξ(1820)
PdgCode omega_bar(-0x3334);

PdgCode lambda_c(0x4122);
PdgCode sigma_c_bar(-0x4114);
PdgCode xi_c(0x4322);
PdgCode omega_c_bar(-0x4332);
PdgCode xi_cc_bar(-0x4422);

PdgCode omega_bc(0x5432);
PdgCode omega_bb(0x5532);

// nuclei
PdgCode deuteron("1000010020");
PdgCode antideutron("-1000010020");
PdgCode triton("1000010030");
PdgCode antitriton("-1000010030");
PdgCode he3("1000020030");
PdgCode antihe3("-1000020030");
PdgCode He4("1000020040");
PdgCode H3L("1010010030");
PdgCode antiH3L("-1010010030");

// non-hadrons
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
  std::printf("π^+:       %8s %8x 0x%08x\n", pi_plus.string().c_str(),
              pi_plus.code(), pi_plus.dump());
  std::printf("K^0:       %8s %8x 0x%08x\n", k0.string().c_str(), k0.code(),
              k0.dump());
  std::printf("K^0L:      %8s %8x 0x%08x\n", k0_long.string().c_str(),
              k0_long.code(), k0_long.dump());
  std::printf("K^0S:      %8s %8x 0x%08x\n", k0_short.string().c_str(),
              k0_short.code(), k0_short.dump());
  std::printf("K^-:       %8s %8x 0x%08x\n", k_minus.string().c_str(),
              k_minus.code(), k_minus.dump());
  std::printf("D^-:       %8s %8x 0x%08x\n", d_minus.string().c_str(),
              d_minus.code(), d_minus.dump());
  std::printf("B^0_s:     %8s %8x 0x%08x\n", bs0.string().c_str(), bs0.code(),
              bs0.dump());
  std::printf("bar B^+_c: %8s %8x 0x%08x\n", bc_minus.string().c_str(),
              bc_minus.code(), bc_minus.dump());
  std::printf("η':        %8s %8x 0x%08x\n", eta_prim.string().c_str(),
              eta_prim.code(), eta_prim.dump());
  std::printf("J/Ψ:       %8s %8x 0x%08x\n", j_psi.string().c_str(),
              j_psi.code(), j_psi.dump());

  std::printf("############################# Baryons:\n");
  std::printf("p:         %8s %8x 0x%08x\n", proton.string().c_str(),
              proton.code(), proton.dump());
  std::printf("bar Δ(1700)%8s %8x 0x%08x\n", delta_bar.string().c_str(),
              delta_bar.code(), delta_bar.dump());
  std::printf("Σ:         %8s %8x 0x%08x\n", sigma_plus.string().c_str(),
              sigma_plus.code(), sigma_plus.dump());
  std::printf("Λ:         %8s %8x 0x%08x\n", lambda.string().c_str(),
              lambda.code(), lambda.dump());
  std::printf("bar Ξ(1820)%8s %8x 0x%08x\n", xi_bar.string().c_str(),
              xi_bar.code(), xi_bar.dump());
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
  COMPARE(pi_plus.code(), 0x211);
  COMPARE(k0.code(), 0x311);
  COMPARE(k0_long.code(), 0x130);
  COMPARE(k0_short.code(), 0x310);
  COMPARE(proton.code(), 0x2212);
  COMPARE(delta_bar.code(), static_cast<int>(0xffeddddcu));
  COMPARE(lambda.code(), 0x3122);
  COMPARE(xi_bar.code(), static_cast<int>(0xffefcceeu));
}
TEST(dump) {
  COMPARE(electron.dump(), 0x11u);
  COMPARE(antimu.dump(), 0x80000013u);
  COMPARE(photon.dump(), 0x22u);
  COMPARE(pi_plus.dump(), 0x211u);
  COMPARE(k0.dump(), 0x311u);
  COMPARE(k0_long.dump(), 0x130u);
  COMPARE(k0_short.dump(), 0x310u);
  COMPARE(proton.dump(), 0x2212u);
  COMPARE(delta_bar.dump(), 0x80122224u);
  COMPARE(lambda.dump(), 0x3122u);
  COMPARE(xi_bar.dump(), 0x80103312u);
}

TEST(pdg_string) {
  COMPARE(electron.string(), "11");
  COMPARE(antimu.string(), "-13");
  COMPARE(photon.string(), "22");
  COMPARE(pi_plus.string(), "211");
  COMPARE(k0.string(), "311");
  COMPARE(k0_long.string(), "130");
  COMPARE(k0_short.string(), "310");
  COMPARE(proton.string(), "2212");
  COMPARE(delta_bar.string(), "-122224");
  COMPARE(lambda.string(), "3122");
  COMPARE(xi_bar.string(), "-103312");
  COMPARE(deuteron.string(), "1000010020");
  COMPARE(antideutron.string(), "-1000010020");
  COMPARE(He4.string(), "1000020040");
  COMPARE(H3L.string(), "1010010030");
}

TEST(decimal) {
  COMPARE(electron.get_decimal(), 11);
  COMPARE(antimu.get_decimal(), -13);
  COMPARE(photon.get_decimal(), 22);
  COMPARE(pi_plus.get_decimal(), 211);
  COMPARE(k0.get_decimal(), 311);
  COMPARE(proton.get_decimal(), 2212);
  COMPARE(delta_bar.get_decimal(), -122224);
  COMPARE(lambda.get_decimal(), 3122);
  COMPARE(xi_bar.get_decimal(), -103312);
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
  VERIFY(pi_plus.is_hadron());
  VERIFY(proton.is_hadron());
  VERIFY(delta_bar.is_hadron());
  VERIFY(!deuteron.is_hadron());
  VERIFY(!antideutron.is_hadron());
}

TEST(nucleus) {
  VERIFY(!electron.is_nucleus());
  VERIFY(!pi_plus.is_nucleus());
  VERIFY(!proton.is_nucleus());
  VERIFY(!delta_bar.is_nucleus());
  VERIFY(deuteron.is_nucleus());
  VERIFY(antideutron.is_nucleus());
  VERIFY(He4.is_nucleus());
  VERIFY(H3L.is_nucleus());
}

TEST(lepton) {
  VERIFY(electron.is_lepton());
  VERIFY(antimu.is_lepton());

  VERIFY(!photon.is_lepton());
  VERIFY(!pi_plus.is_lepton());
  VERIFY(!proton.is_lepton());
  VERIFY(!deuteron.is_lepton());
}
TEST(is_meson) {
  VERIFY(pi_plus.is_meson());
  VERIFY(!deuteron.is_meson());
  VERIFY(!antideutron.is_meson());
  VERIFY(!proton.is_meson());
  VERIFY(!proton_bar.is_meson());
  VERIFY(!antimu.is_meson());
}

TEST(is_baryon) {
  VERIFY(!pi_plus.is_baryon());
  VERIFY(!deuteron.is_baryon());
  VERIFY(!antideutron.is_baryon());
  VERIFY(proton.is_baryon());
  VERIFY(proton_bar.is_baryon());
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
  COMPARE(pi_plus.baryon_number(), 0);
  COMPARE(k0.baryon_number(), 0);
  COMPARE(k0_long.baryon_number(), 0);
  COMPARE(k0_short.baryon_number(), 0);
  COMPARE(k_minus.baryon_number(), 0);
  COMPARE(d_minus.baryon_number(), 0);
  COMPARE(b0.baryon_number(), 0);
  COMPARE(bc_minus.baryon_number(), 0);
  COMPARE(eta_prim.baryon_number(), 0);
  COMPARE(j_psi.baryon_number(), 0);
  COMPARE(proton.baryon_number(), 1);
  COMPARE(delta_bar.baryon_number(), -1);
  COMPARE(sigma.baryon_number(), 1);
  COMPARE(lambda.baryon_number(), 1);
  COMPARE(xi_bar.baryon_number(), -1);
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
  COMPARE(pi_plus.isospin3(), +2);
  COMPARE(k0.isospin3(), -1);
  COMPARE(k_minus.isospin3(), -1);
  COMPARE(d_minus.isospin3(), -1);
  COMPARE(bs0.isospin3(), 0);
  COMPARE(bc_minus.isospin3(), 0);
  COMPARE(eta_prim.isospin3(), 0);
  COMPARE(j_psi.isospin3(), 0);
  COMPARE(proton.isospin3(), 1);
  COMPARE(delta_bar.isospin3(), -3);
  COMPARE(sigma_plus.isospin3(), +2);
  COMPARE(lambda.isospin3(), 0);
  COMPARE(xi_bar.isospin3(), +1);
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
  COMPARE(pi_plus.strangeness(), 0);
  COMPARE(k0.strangeness(), 1);
  COMPARE(k_minus.strangeness(), -1);
  COMPARE(d_minus.strangeness(), 0);
  COMPARE(bs0.strangeness(), -1);
  COMPARE(bc_minus.strangeness(), 0);
  COMPARE(eta_prim.strangeness(), 0);
  COMPARE(j_psi.strangeness(), 0);
  COMPARE(proton.strangeness(), 0);
  COMPARE(delta_bar.strangeness(), 0);
  COMPARE(sigma.strangeness(), -1);
  COMPARE(lambda.strangeness(), -1);
  COMPARE(xi_bar.strangeness(), +2);
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
  COMPARE(pi_plus.charmness(), 0);
  COMPARE(k0.charmness(), 0);
  COMPARE(k_minus.charmness(), 0);
  COMPARE(d_minus.charmness(), -1);
  COMPARE(b0.charmness(), 0);
  COMPARE(bc_minus.charmness(), -1);
  COMPARE(eta_prim.charmness(), 0);
  COMPARE(j_psi.charmness(), 0);
  COMPARE(proton.charmness(), 0);
  COMPARE(delta_bar.charmness(), 0);
  COMPARE(sigma.charmness(), 0);
  COMPARE(lambda.charmness(), 0);
  COMPARE(xi_bar.charmness(), 0);
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
  COMPARE(pi_plus.bottomness(), 0);
  COMPARE(k0.bottomness(), 0);
  COMPARE(k_minus.bottomness(), 0);
  COMPARE(d_minus.bottomness(), 0);
  COMPARE(b0.bottomness(), 1);
  COMPARE(bc_minus.bottomness(), -1);
  COMPARE(eta_prim.bottomness(), 0);
  COMPARE(j_psi.bottomness(), 0);
  COMPARE(proton.bottomness(), 0);
  COMPARE(delta_bar.bottomness(), 0);
  COMPARE(sigma.bottomness(), 0);
  COMPARE(lambda.bottomness(), 0);
  COMPARE(xi_bar.bottomness(), 0);
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
  COMPARE(pi_plus.frac_strange(), 0);
  COMPARE(k0.frac_strange(), 1. / 2);
  COMPARE(k_minus.frac_strange(), 1. / 2);
  COMPARE(d_minus.frac_strange(), 0);
  COMPARE(bs0.frac_strange(), 1. / 2);
  COMPARE(bc_minus.frac_strange(), 0);
  COMPARE(eta_prim.frac_strange(), 1);
  COMPARE(j_psi.frac_strange(), 0);
  COMPARE(proton.frac_strange(), 0);
  COMPARE(delta_bar.frac_strange(), 0);
  COMPARE(sigma.frac_strange(), 1. / 3);
  COMPARE(lambda.frac_strange(), 1. / 3);
  COMPARE(xi_bar.frac_strange(), 2. / 3);
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
  COMPARE(pi_plus.frac_bottom(), 0);
  COMPARE(k0.frac_bottom(), 0);
  COMPARE(k_minus.frac_bottom(), 0);
  COMPARE(d_minus.frac_bottom(), 0);
  COMPARE(b0.frac_bottom(), 1. / 2);
  COMPARE(bc_minus.frac_bottom(), 1. / 2);
  COMPARE(eta_prim.frac_bottom(), 0);
  COMPARE(j_psi.frac_bottom(), 0);
  COMPARE(proton.frac_bottom(), 0);
  COMPARE(delta_bar.frac_bottom(), 0);
  COMPARE(sigma.frac_bottom(), 0);
  COMPARE(lambda.frac_bottom(), 0);
  COMPARE(xi_bar.frac_bottom(), 0);
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
  COMPARE(pi_plus.frac_charm(), 0);
  COMPARE(k0.frac_charm(), 0);
  COMPARE(k_minus.frac_charm(), 0);
  COMPARE(d_minus.frac_charm(), 1. / 2);
  COMPARE(b0.frac_charm(), 0);
  COMPARE(bc_minus.frac_charm(), 1. / 2);
  COMPARE(eta_prim.frac_charm(), 0);
  COMPARE(j_psi.frac_charm(), 1);
  COMPARE(proton.frac_charm(), 0);
  COMPARE(delta_bar.frac_charm(), 0);
  COMPARE(sigma.frac_charm(), 0);
  COMPARE(lambda.frac_charm(), 0);
  COMPARE(xi_bar.frac_charm(), 0);
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
  COMPARE(pi_plus.is_heavy_flavor(), 0);
  COMPARE(k0.is_heavy_flavor(), 0);
  COMPARE(k_minus.is_heavy_flavor(), 0);
  COMPARE(d_minus.is_heavy_flavor(), 1);
  COMPARE(b0.is_heavy_flavor(), 1);
  COMPARE(bc_minus.is_heavy_flavor(), 1);
  COMPARE(eta_prim.is_heavy_flavor(), 0);
  COMPARE(j_psi.is_heavy_flavor(), 1);
  COMPARE(proton.is_heavy_flavor(), 0);
  COMPARE(delta_bar.is_heavy_flavor(), 0);
  COMPARE(sigma.is_heavy_flavor(), 0);
  COMPARE(lambda.is_heavy_flavor(), 0);
  COMPARE(xi_bar.is_heavy_flavor(), 0);
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
  COMPARE(pi_plus.charge(), +1);
  COMPARE(k0.charge(), 0);
  COMPARE(k0_long.charge(), 0);
  COMPARE(k0_short.charge(), 0);
  COMPARE(k_minus.charge(), -1);
  COMPARE(d_minus.charge(), -1);
  COMPARE(b0.charge(), 0);
  COMPARE(bc_minus.charge(), -1);
  COMPARE(eta_prim.charge(), 0);
  COMPARE(j_psi.charge(), 0);
  COMPARE(proton.charge(), +1);
  COMPARE(delta_bar.charge(), -2);
  COMPARE(sigma_plus.charge(), +1);
  COMPARE(lambda.charge(), 0);
  COMPARE(xi_bar.charge(), +1);
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
  COMPARE(pi_plus.quarks(), 0x021);
  COMPARE(k0.quarks(), 0x031);
  COMPARE(k0_long.quarks(), 0x013);
  COMPARE(k0_short.quarks(), 0x031);
  COMPARE(k_minus.quarks(), 0x032);
  COMPARE(d_minus.quarks(), 0x041);
  COMPARE(bs0.quarks(), 0x053);
  COMPARE(bc_minus.quarks(), 0x054);
  COMPARE(eta_prim.quarks(), 0x033);
  COMPARE(j_psi.quarks(), 0x044);
  COMPARE(proton.quarks(), 0x221);
  COMPARE(delta_bar.quarks(), 0x222);
  COMPARE(sigma_plus.quarks(), 0x322);
  COMPARE(lambda.quarks(), 0x312);
  COMPARE(xi_bar.quarks(), 0x331);
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
  COMPARE(pi_plus.spin(), 0u);
  COMPARE(k0.spin(), 0u);
  COMPARE(k0_long.spin(), 0u);
  COMPARE(k0_short.spin(), 0u);
  COMPARE(k_minus.spin(), 0u);
  COMPARE(d_minus.spin(), 0u);
  COMPARE(b0.spin(), 0u);
  COMPARE(bc_minus.spin(), 0u);
  COMPARE(eta_prim.spin(), 0u);
  COMPARE(j_psi.spin(), 2u);
  COMPARE(proton.spin(), 1u);
  COMPARE(delta_bar.spin(), 3u);
  COMPARE(sigma.spin(), 1u);
  COMPARE(lambda.spin(), 1u);
  COMPARE(xi_bar.spin(), 1u);
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
  COMPARE(pi_plus.spin_degeneracy(), 1u);
  COMPARE(k0.spin_degeneracy(), 1u);
  COMPARE(k0_long.spin_degeneracy(), 1u);
  COMPARE(k0_short.spin_degeneracy(), 1u);
  COMPARE(k_minus.spin_degeneracy(), 1u);
  COMPARE(d_minus.spin_degeneracy(), 1u);
  COMPARE(b0.spin_degeneracy(), 1u);
  COMPARE(bc_minus.spin_degeneracy(), 1u);
  COMPARE(eta_prim.spin_degeneracy(), 1u);
  COMPARE(j_psi.spin_degeneracy(), 3u);
  COMPARE(proton.spin_degeneracy(), 2u);
  COMPARE(delta_bar.spin_degeneracy(), 4u);
  COMPARE(sigma.spin_degeneracy(), 2u);
  COMPARE(lambda.spin_degeneracy(), 2u);
  COMPARE(xi_bar.spin_degeneracy(), 2u);
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
  VERIFY(pi_plus != eta_prim);
  PdgCode pi_plus2(0x211);
  VERIFY(pi_plus == pi_plus2);
  VERIFY(pi_plus2 < omega_bc);
}
TEST(antiparticle) { VERIFY(pi_plus.is_antiparticle_of(pi_minus)); }
TEST(get_antiparticle) { COMPARE(pi_plus.get_antiparticle(), pi_minus); }

TEST(from_decimal) {
  COMPARE(pi_plus, PdgCode::from_decimal(211));
  COMPARE(k_minus, PdgCode::from_decimal(-321));
  COMPARE(xi_bar, PdgCode::from_decimal(-103312));
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
  COMPARE(pi_plus.has_antiparticle(), true);
  COMPARE(pi0.has_antiparticle(), false);
  COMPARE(k0.has_antiparticle(), true);
  COMPARE(proton.has_antiparticle(), true);
  COMPARE(delta_pp.has_antiparticle(), true);
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
  std::array<int, 3> q;
  q = pi_plus.quark_content();
  COMPARE(q[0], 0);
  COMPARE(q[1], 2);
  COMPARE(q[2], -1);

  q = pi_minus.quark_content();
  COMPARE(q[0], 0);
  COMPARE(q[1], -2);
  COMPARE(q[2], 1);

  q = pi0.quark_content();
  COMPARE(q[0], 0);
  COMPARE(q[1], 1);
  COMPARE(q[2], -1);

  q = neutron.quark_content();
  COMPARE(q[0], 2);
  COMPARE(q[1], 1);
  COMPARE(q[2], 1);

  q = proton_bar.quark_content();
  COMPARE(q[0], -2);
  COMPARE(q[1], -2);
  COMPARE(q[2], -1);

  q = neutron_bar.quark_content();
  COMPARE(q[0], -2);
  COMPARE(q[1], -1);
  COMPARE(q[2], -1);

  q = electron.quark_content();
  COMPARE(q[0], 0);
  COMPARE(q[1], 0);
  COMPARE(q[2], 0);
}

TEST(net_quark_number) {
  // pi_plus+ has one u and one dbar
  VERIFY(pi_plus.net_quark_number(1) == -1);
  VERIFY(pi_plus.net_quark_number(2) == 1);

  // pi_plus- has one d and one ubar
  VERIFY(pi_minus.net_quark_number(1) == 1);
  VERIFY(pi_minus.net_quark_number(2) == -1);

  // pi_plus0 has one vanishing net u and d quark numbers
  VERIFY(pi0.net_quark_number(1) == 0);
  VERIFY(pi0.net_quark_number(2) == 0);

  // proton has two u and one d
  VERIFY(proton.net_quark_number(1) == 1);
  VERIFY(proton.net_quark_number(2) == 2);

  // neutron has one u and two d
  VERIFY(neutron.net_quark_number(1) == 2);
  VERIFY(neutron.net_quark_number(2) == 1);

  // proton_bar has two ubar and one dbar
  VERIFY(proton_bar.net_quark_number(1) == -1);
  VERIFY(proton_bar.net_quark_number(2) == -2);

  // antineutron has one ubar and two dbar
  VERIFY(neutron_bar.net_quark_number(1) == -2);
  VERIFY(neutron_bar.net_quark_number(2) == -1);
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
  // pi_plus0, meson, baryon number of 0
  VERIFY(pi0.contains_enough_valence_quarks(1));
  VERIFY(pi0.contains_enough_valence_quarks(-1));

  VERIFY(neutron_bar.contains_enough_valence_quarks(-2));
  VERIFY(neutron_bar.contains_enough_valence_quarks(-1));

  // proton, baryon, baryon number of 1
  VERIFY(proton.contains_enough_valence_quarks(1));
  VERIFY(proton.contains_enough_valence_quarks(2));

  // pi_plus-, meson, baryon number of 0
  VERIFY(!(pi_minus.contains_enough_valence_quarks(2)));
  VERIFY(!(pi_minus.contains_enough_valence_quarks(-2)));
}

TEST(nucleus_components) {
  VERIFY(H3L.nucleus_p() == 1);
  VERIFY(H3L.nucleus_n() == 1);
  VERIFY(H3L.nucleus_La() == 1);
  VERIFY(H3L.nucleus_ap() == 0);
  VERIFY(H3L.nucleus_an() == 0);
  VERIFY(H3L.nucleus_aLa() == 0);
}
namespace {

static void verify_quark_content(const PdgCode& hadron,
                                 std::initializer_list<int> present,
                                 std::initializer_list<int> absent) {
  for (const int q : present) {
    VERIFY(hadron.contains_quark(q));
  }
  for (const int q : absent) {
    VERIFY(!hadron.contains_quark(q));
  }
}

}  // namespace
TEST(contains_quark) {
  verify_quark_content(pi_minus, {1, -2}, {-1, 2, 3, -3});
  verify_quark_content(pi_plus, {-1, 2}, {1, -2, 3, -3});

  verify_quark_content(k_plus, {2, -3}, {-2, 3, 1, -1});
  verify_quark_content(k_minus, {-2, 3}, {2, -3, 1, -1});

  verify_quark_content(proton_bar, {-1, -2}, {1, 2, 3, -3});

  verify_quark_content(neutron, {1, 2}, {-1, -2, 3, -3});
  verify_quark_content(neutron_bar, {-1, -2}, {1, 2, 3, -3});

  verify_quark_content(lambda, {1, 2, 3}, {-1, -2, -3});
  verify_quark_content(lambda_bar, {-1, -2, -3}, {1, 2, 3});

  verify_quark_content(sigma_plus, {2, 3}, {1, -2, -3});
  verify_quark_content(sigma_minus, {1, 3}, {2, -1, -3});

  verify_quark_content(delta_pp, {2}, {1, 3, -2});
  verify_quark_content(delta_minus, {1}, {2, 3, -1});

  verify_quark_content(d_plus, {4, -1}, {-4, 1});
  verify_quark_content(d_minus, {-4, 1}, {4, -1});

  verify_quark_content(pi0, {1, -1, 2, -2}, {3, -3});
  verify_quark_content(rho0, {1, -1, 2, -2}, {3, -3});
  verify_quark_content(omega, {1, -1, 2, -2}, {3, -3});

  verify_quark_content(eta, {1, -1, 2, -2, 3, -3}, {});
  verify_quark_content(eta_prim, {1, -1, 2, -2, 3, -3}, {});

  verify_quark_content(phi, {3, -3}, {1, -1, 2, -2});
  verify_quark_content(j_psi, {4, -4}, {1, -1, 2, -2, 3, -3, 5, -5, 6, -6});
  verify_quark_content(upsilon, {5, -5}, {3, -3, 4, -4});

  verify_quark_content(k0, {1, -3}, {-1, 3});
  verify_quark_content(k0_bar, {-1, 3}, {1, -3});

  verify_quark_content(d0, {4, -2}, {-4, 2});
  verify_quark_content(d0_bar, {-4, 2}, {4, -2});

  verify_quark_content(b0, {1, -5}, {-1, 5});
  verify_quark_content(b0_bar, {-1, 5}, {1, -5});

  verify_quark_content(bs0, {3, -5}, {-3, 5});
  verify_quark_content(bs0_bar, {-3, 5}, {3, -5});
}
