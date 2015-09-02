/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteractionnucleonkaon.h"

#include "include/cxx14compat.h"
#include "include/parametrizations.h"

namespace Smash {

float ScatterActionNucleonKaon::elastic_parametrization() {
  const PdgCode &pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode &pdg_b = incoming_particles_[1].type().pdgcode();

  const PdgCode &nucleon = pdg_a.is_nucleon() ? pdg_a : pdg_b;
  const PdgCode &kaon = pdg_a.is_nucleon() ? pdg_b : pdg_a;
  assert(kaon != nucleon);

  const double s = mandelstam_s();

  float sig_el;
  switch (nucleon.code()) {
    case 0x2212:  // p
      switch (kaon.code()) {
        case 0x321:  // K+
          sig_el = kplusp_elastic(s);
          break;
        case -0x321:  // K-
          sig_el = kminusp_elastic(s);
          break;
        case 0x311:  // K0
          sig_el = k0p_elastic(s);
          break;
        case -0x311:  // Kbar0
          sig_el = kbar0p_elastic(s);
          break;
      }
      break;
    case 0x2112:  // n
      switch (kaon.code()) {
        case 0x321:  // K+
          sig_el = kplusn_elastic(s);
          break;
        case -0x321:  // K-
          sig_el = kminusn_elastic(s);
          break;
        case 0x311:  // K0
          sig_el = k0n_elastic(s);
          break;
        case -0x311:  // Kbar0
          sig_el = kbar0n_elastic(s);
          break;
      }
      break;
    default:
      throw std::runtime_error("elastic cross section for antinucleon-kaon not implemented");
  }

  if (sig_el > 0) {
    return sig_el;
  } else {
    std::stringstream ss;
    const auto name_a = incoming_particles_[0].type().name();
    const auto name_b = incoming_particles_[1].type().name();
    ss << "problem in CrossSections::elastic: a=" << name_a
       << " b=" << name_b << " j_a=" << pdg_a.spin() << " j_b="
       << pdg_b.spin() << " sigma=" << sig_el << " s=" << s;
    throw std::runtime_error(ss.str());
  }
}

void ScatterActionNucleonKaon::format_debug_output(std::ostream &out) const {
  out << "Nucleon-Kaon  ";
  ScatterAction::format_debug_output(out);
}


}  // namespace Smash
