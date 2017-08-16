/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteractionnucleonpion.h"

#include "include/clebschgordan.h"
#include "include/cxx14compat.h"
#include "include/parametrizations.h"
#include "include/pdgcode_constants.h"

namespace Smash {

double ScatterActionNucleonPion::elastic_parametrization() {
  const PdgCode &pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode &pdg_b = incoming_particles_[1].type().pdgcode();

  const PdgCode &nucleon = pdg_a.is_nucleon() ? pdg_a : pdg_b;
  const PdgCode &pion = pdg_a.is_nucleon() ? pdg_b : pdg_a;
  assert(pion != nucleon);

  const double s = mandelstam_s();

  double sig_el = 0.;
  switch (nucleon.code()) {
    case pdg::p:
      switch (pion.code()) {
        case pdg::pi_p:
          sig_el = piplusp_elastic(s);
          break;
        case pdg::pi_m:
          sig_el = piminusp_elastic(s);
          break;
        case pdg::pi_z:
          sig_el = 0.5 * (piplusp_elastic(s) + piminusp_elastic(s));
          break;
      }
      break;
    case pdg::n:
      switch (pion.code()) {
        case pdg::pi_p:
          sig_el = piminusp_elastic(s);
          break;
        case pdg::pi_m:
          sig_el = piplusp_elastic(s);
          break;
        case pdg::pi_z:
          sig_el = 0.5 * (piplusp_elastic(s) + piminusp_elastic(s));
          break;
      }
      break;
    case -pdg::p:
      switch (pion.code()) {
        case pdg::pi_p:
          sig_el = piminusp_elastic(s);
          break;
        case pdg::pi_m:
          sig_el = piplusp_elastic(s);
          break;
        case pdg::pi_z:
          sig_el = 0.5 * (piplusp_elastic(s) + piminusp_elastic(s));
          break;
      }
      break;
    case -pdg::n:
      switch (pion.code()) {
        case pdg::pi_p:
          sig_el = piplusp_elastic(s);
          break;
        case pdg::pi_m:
          sig_el = piminusp_elastic(s);
          break;
        case pdg::pi_z:
          sig_el = 0.5 * (piplusp_elastic(s) + piminusp_elastic(s));
          break;
      }
      break;
    default:
      throw std::runtime_error(
          "only the elastic cross section for proton-pion "
          "is implemented");
  }

  if (sig_el > 0) {
    return sig_el;
  } else {
    std::stringstream ss;
    const auto name_a = incoming_particles_[0].type().name();
    const auto name_b = incoming_particles_[1].type().name();
    ss << "problem in CrossSections::elastic: a=" << name_a << " b=" << name_b
       << " j_a=" << pdg_a.spin() << " j_b=" << pdg_b.spin()
       << " sigma=" << sig_el << " s=" << s;
    throw std::runtime_error(ss.str());
  }
}

void ScatterActionNucleonPion::format_debug_output(std::ostream &out) const {
  out << "Nucleon-Pion  ";
  ScatterAction::format_debug_output(out);
}

}  // namespace Smash
