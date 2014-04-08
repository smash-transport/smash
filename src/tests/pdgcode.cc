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

constexpr double accuracy = 1e-10;

TEST(set_code) {
  printf ("%08x\n", (0xf << 8) | 0xf0);


  PDGCode pion(0x211);
  printf ("pion: Anti: %08x\n", pion.antiparticle_);
  printf ("pion: n   : %08x\n", pion.n_);
  printf ("pion: nR  : %08x\n", pion.n_R_);
  printf ("pion: nL  : %08x\n", pion.n_L_);
  printf ("pion: nq1 : %08x\n", pion.n_q1_);
  printf ("pion: nq2 : %08x\n", pion.n_q2_);
  printf ("pion: nq3 : %08x\n", pion.n_q3_);
  printf ("pion: nJ  : %08x\n", pion.n_J_);
  printf ("pion: %08x\n", pion.code());

  FUZZY_COMPARE(pion.code(), 0x211) << "Pion";
  PDGCode proton(0x2212);
  FUZZY_COMPARE(proton.code(), 0x2212) << "Proton";
  PDGCode anti_proton(-0x2212);
  FUZZY_COMPARE(anti_proton.code(), ((2<<31)+0x2212)) << "Antiproton";
}

TEST_CATCH(set_invalid_code, PDGCode::InvalidPDGCode) {
  PDGCode code(211);
}
