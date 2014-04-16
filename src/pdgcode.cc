/*
 * Copyright (c) 2013
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include<istream>
#include "include/pdgcode.h"

namespace Smash {

std::istream& operator>>(std::istream& is, PdgCode& code) {
  std::string codestring("");
  // discard any whitespace at beginning:
  while (is.peek() == ' ' || is.peek() == '\t') {
    is.get();
  }
  // read sign if there is one:
  if (is.peek() == '+' || is.peek() == '-') {
    codestring += is.get();
  }
  // read a maximum of 7 characters
  for (int c = 0; c < 7; c++) {
    // look into the string. Is it a valid character?
    try {
      code.get_digit_from_char(is.peek());
    } catch (PdgCode::InvalidPdgCode) {
      // if not, end the loop.
      break;
    }
    // read one character from is:
    // char * s = new char[1];
    // is.read(s, 1);
    codestring += is.get();
  }
  try {
    // set the fields from the string:
    code.set_from_string(codestring);
  } catch (PdgCode::InvalidPdgCode) {
    is.setstate(std::ios::failbit);
  }
  // get as much whitespace as possible:
  while (is.peek() == ' ' || is.peek() == '\t') {
    is.get();
  }
  return is;
}

unsigned int PdgCode::isospin_total() const {
  // non-hadrons and η mesons (and ω and stuff):
  if (!is_hadron() || quarks() == 0x220) {
    return 0;
  }
  int number_of_u_or_d_quarks = 0;
  if (n_q3_ == 2 || n_q3_ == 1) { number_of_u_or_d_quarks++; }
  if (n_q2_ == 2 || n_q2_ == 1) { number_of_u_or_d_quarks++; }
  if (n_q1_ == 2 || n_q1_ == 1) { number_of_u_or_d_quarks++; }
  // Δ and N distinction. I don't know any smart algorithm for this; I am
  // confident that special casing is the only way to do that.
  if (number_of_u_or_d_quarks == 3) {
    std::int32_t multi = std::abs(multiplett());
    // first the most common ones: p/n and Δ
    if (multi == 0x1002) { return 1; }
    if (multi == 0x1004) { return 3; }
    // now the resonances:
    if (multi == 0x101002
     || multi == 0x121002
     || multi == 0x201002
     || multi == 0x211002
     || multi == 0x101004
     || multi == 0x111004
     || multi == 0x201004
     || multi == 0x211004
     || multi == 0x101006
     || multi == 0x201006
       ) { return 1; }
    if (multi == 0x111002
     || multi == 0x221002
     || multi == 0x121004
     || multi == 0x221004
     || multi == 0x211006
     || multi == 0x201008
       ) { return 3; }
    throw InvalidPdgCode(
          "Unknown Nucleon or Delta resonance, cannot determine isospin: "
          + string());
  }
  // special case: Λ. We know already that not three quarks are up or down,
  // so we need to find where the third is:
  // 312, 412, 512 are Λ, as is 213. And, well, while this is not yet part of
  // the standard, I'll throw in 214 and 215 as well.
  if ((quarks() & 0x0fff) == 0x0120 || (quarks() & 0xff0f) == 0x2100) {
    return 0;
  }
  return number_of_u_or_d_quarks;
}

int PdgCode::heavy_quarkness(const int quark) const {
  // input sanitization: Only quark numbers 1 through 8 are allowed.
  if (quark < 1 || quark > 8) {
    throw std::invalid_argument(std::string("PdgCode::heavy_quarkness(): ")
                 + std::string("Quark number must be in [1..8], received ")
                 + std::to_string(quark));
  }
  // non-hadrons and those that have none of this quark type: 0.
  if (! is_hadron() ||
        (n_q1_ != quark && n_q2_ != quark && n_q3_ != quark)) {
    return 0;
  }
  // baryons: count quarks.
  if (baryon_number() != 0) {
    // u,s,b quarks get negative *ness; d,c,t have positive *ness:
    int sign = (quark%2 == 0) ? +1 : -1;
    // and for anti-baryons, the sign changes:
    sign *= antiparticle_sign();
    return sign*((n_q1_ == quark) + (n_q2_ == quark) + (n_q3_ == quark));
  }
  // mesons.
  // quarkonium state? Not open heavy-quarkness.
  if (n_q3_ == quark && n_q2_ == quark) {
    return 0;
  }
  // this has covered all the easy stuff
  // get the "other" quark. (We know this must exist, since they are
  // not both the right one and one of them is the right one).
  int otherquark = (n_q2_ == quark) ? n_q3_ : n_q2_;
  // "our" quark is the heavier one: 1 for particles
  if (otherquark < quark) {
    return antiparticle_sign();
  }
  // ours is the lighter: If the heavier quark has the same SIGN, we
  // get a -1, else +1.
  if (otherquark%2 == quark%2) {
    return -1*antiparticle_sign();
  }
  return antiparticle_sign();
}

}  // namespace Smash
