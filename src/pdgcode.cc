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

}  // namespace Smash
