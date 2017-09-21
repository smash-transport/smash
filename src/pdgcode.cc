/*
 * Copyright (c) 2013-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include <istream>

#include "include/pdgcode.h"

namespace Smash {

std::istream& operator>>(std::istream& is, PdgCode& code) {
  std::string codestring("");
  is >> codestring;
  if (!is) {
    code = PdgCode::invalid();
    return is;
  }
  try {
    // set the fields from the string:
    code.set_from_string(codestring);
  } catch (PdgCode::InvalidPdgCode) {
    is.setstate(std::ios::failbit);
    code = PdgCode::invalid();
  }
  return is;
}

int PdgCode::net_quark_number(const int quark) const {
  // input sanitization: Only quark numbers 1 through 8 are allowed.
  if (quark < 1 || quark > 8) {
    throw std::invalid_argument(
        std::string("PdgCode::net_quark_number(): ") +
        std::string("Quark number must be in [1..8], received ") +
        std::to_string(quark));
  }
  // non-hadrons and those that have none of this quark type: 0.
  if (!is_hadron() || (digits_.n_q1_ != quark && digits_.n_q2_ != quark &&
                       digits_.n_q3_ != quark)) {
    return 0;
  }
  // baryons: count quarks.
  if (baryon_number() != 0) {
    // for anti-baryons, the sign changes:
    return antiparticle_sign() *
           ((digits_.n_q1_ == quark) + (digits_.n_q2_ == quark) +
            (digits_.n_q3_ == quark));
  }
  // mesons.
  // quarkonium state? Not open net_quark_number.
  if (digits_.n_q3_ == quark && digits_.n_q2_ == quark) {
    return 0;
  }
  // this has covered all the easy stuff
  // get the "other" quark. (We know this must exist, since they are
  // not both the right one and one of them is the right one).
  int otherquark = (digits_.n_q2_ == quark) ? digits_.n_q3_ : digits_.n_q2_;
  // "our" quark is the heavier one: 1 for u,c,t; -1 for d,s,b (and of
  // course the antiparticle sign)
  if (quark > otherquark) {
    return ((quark % 2 == 0) ? 1 : -1) * antiparticle_sign();
  }
  // ours is the lighter: If the heavier particle is u,c,t, the lighter
  // one (ours) is an antiquark.
  return ((otherquark % 2 == 0) ? -1 : 1) * antiparticle_sign();
}

std::ostream& operator<<(std::ostream& s, const PdgCode& code) {
  return s << code.string();
}

}  // namespace Smash
