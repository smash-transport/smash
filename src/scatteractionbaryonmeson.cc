/*
 *
 *    Copyright (c) 2015-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/clebschgordan.h"
#include "include/scatteractionbaryonmeson.h"
#include "include/scatteractionbaryonbaryon.h"
#include "include/parametrizations.h"
#include "include/pow.h"

namespace smash {

void ScatterActionBaryonMeson::format_debug_output(std::ostream &out) const {
  out << "Baryon-Meson  ";
  ScatterAction::format_debug_output(out);
}

double ScatterActionBaryonMeson::high_energy_cross_section() const {
  const PdgCode &pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode &pdg_b = incoming_particles_[1].type().pdgcode();
  const double s = mandelstam_s();

  /* Currently only include pion nucleon interaction. */
  if ((pdg_a == pdg::pi_p && pdg_b == pdg::p) ||
      (pdg_b == pdg::pi_p && pdg_a == pdg::p) ||
      (pdg_a == pdg::pi_m && pdg_b == pdg::n) ||
      (pdg_b == pdg::pi_m && pdg_a == pdg::n)) {
    return piplusp_high_energy(s);  // pi+ p, pi- n
  } else if ((pdg_a == pdg::pi_m && pdg_b == pdg::p) ||
             (pdg_b == pdg::pi_m && pdg_a == pdg::p) ||
             (pdg_a == pdg::pi_p && pdg_b == pdg::n) ||
             (pdg_b == pdg::pi_p && pdg_a == pdg::n)) {
    return piminusp_high_energy(s);  // pi- p, pi+ n
  } else {
    return 0;
  }
}

double ScatterActionBaryonMeson::string_hard_cross_section() const {
  // const PdgCode &pdg_a = incoming_particles_[0].type().pdgcode();
  // const PdgCode &pdg_b = incoming_particles_[1].type().pdgcode();
  const double s = mandelstam_s();

  /**
   * Currently nucleon-pion cross section is used for all case.
   * This will be changed later by applying additive quark model.
   */
  return Npi_string_hard(s);
}

CollisionBranchList ScatterActionBaryonMeson::two_to_two_cross_sections(
  ReactionsBitSet /*included_2to2*/) {
  CollisionBranchList process_list;
  const ParticleType &type_a = incoming_particles_[0].type();
  const ParticleType &type_b = incoming_particles_[1].type();

  // pi d -> N N
  if ((type_a.is_deuteron() && type_b.pdgcode().is_pion()) ||
      (type_b.is_deuteron() && type_a.pdgcode().is_pion())) {
    const int baryon_number = type_a.baryon_number() + type_b.baryon_number();
    ParticleTypePtrList nuc = (baryon_number > 0) ?
                                  ParticleType::list_nucleons() :
                                  ParticleType::list_anti_nucleons();
    const double s = mandelstam_s();
    const double sqrts = std::sqrt(s);
    for (ParticleTypePtr nuc_a : nuc) {
      for (ParticleTypePtr nuc_b : nuc) {
        if (type_a.charge() + type_b.charge() !=
             nuc_a->charge() + nuc_b->charge()) {
          continue;
        }
        // loop over total isospin
        for (const int twoI : I_tot_range(*nuc_a, *nuc_b)) {
          const double isospin_factor = isospin_clebsch_gordan_sqr_2to2(
              type_a, type_b, *nuc_a, *nuc_b, twoI);
          /* If Clebsch-Gordan coefficient is zero, don't bother with the rest */
          if (std::abs(isospin_factor) < really_small) {
            continue;
          }

          /* Calculate matrix element for inverse process. */
          const double matrix_element =
              ScatterActionBaryonBaryon::nn_to_resonance_matrix_element(
                  sqrts, type_a, type_b, twoI);
          if (matrix_element <= 0.) {
            continue;
          }

          const double spin_factor = (nuc_a->spin() + 1) * (nuc_b->spin() + 1);
          const int sym_fac_in =
              (type_a.iso_multiplet() == type_b.iso_multiplet()) ? 2 : 1;
          const int sym_fac_out =
              (nuc_a->iso_multiplet() == nuc_b->iso_multiplet()) ? 2 : 1;
          double p_cm_final = pCM_from_s(s, nuc_a->mass(), nuc_b->mass());
          const double xsection = isospin_factor * spin_factor * sym_fac_in /
                                  sym_fac_out * p_cm_final * matrix_element /
                                  (s * cm_momentum());

          if (xsection > really_small) {
            process_list.push_back(make_unique<CollisionBranch>(
                *nuc_a, *nuc_b, xsection, ProcessType::TwoToTwo));
            const auto &log = logger<LogArea::ScatterAction>();
            log.debug(type_a.name(), type_b.name(), "->",
                     nuc_a->name(), nuc_b->name(), " at sqrts [GeV] = ",
                     sqrts, " with cs[mb] = ", xsection);
          }
        }
      }
    }
  }

  // pi d -> pi d' (effectively pi d -> pi p n)  AND reverse, pi d' -> pi d
  if (( (type_a.is_deuteron() || type_a.is_dprime()) &&
         type_b.pdgcode().is_pion() ) ||
      ( (type_b.is_deuteron() || type_b.is_dprime()) &&
         type_a.pdgcode().is_pion()) ) {
    const ParticleType &type_pi = type_a.pdgcode().is_pion() ? type_a : type_b;
    const ParticleType &type_nucleus = type_a.is_nucleus() ? type_a : type_b;
    ParticleTypePtrList nuclei = ParticleType::list_light_nuclei();
    const double s = mandelstam_s();
    const double sqrts = std::sqrt(s);
    for (ParticleTypePtr produced_nucleus : nuclei) {
      // No elastic collisions for now
      if (produced_nucleus == &type_nucleus ||
          produced_nucleus->charge() != type_nucleus.charge() ||
          produced_nucleus->baryon_number() != type_nucleus.baryon_number()) {
        continue;
      }
      // same matrix element for πd and πd̅
      const double tmp = sqrts - type_a.min_mass_kinematic() -
                                 type_b.min_mass_kinematic();
      // Purely experimental-data-based, ad hoc parametrization
      const double matrix_element = 295.5 +
          2.862 / (0.00283735 + pow_int(sqrts - 2.181, 2)) +
          0.0672 / pow_int(tmp, 2) - 6.61753 / tmp;
      const double spin_factor = (produced_nucleus->spin() + 1) *
                                 (type_pi.spin() + 1);
      // Isospin factor is always the same, so it is included into matrix element
      // Symmetry factor is always 1 here
      // Absorb (hbarc)^2/16 pi factor into matrix element
      double xsection = matrix_element * spin_factor / (s * cm_momentum());
      if (produced_nucleus->is_stable()) {
        assert(!type_nucleus.stable());
        xsection *= pCM_from_s(s, type_pi.mass(), produced_nucleus->mass());
      } else {
        assert(type_nucleus.stable());
        const double resonance_integral =
            produced_nucleus->iso_multiplet()->get_integral_piR(sqrts);
        xsection *= resonance_integral;
        //std::cout << "Resonance integral " << resonance_integral
        //          << ", matrix element: " << matrix_element
        //          << ", cm_momentum: " << cm_momentum() << std::endl;
      }
      process_list.push_back(make_unique<CollisionBranch>(
            type_pi, *produced_nucleus, xsection, ProcessType::TwoToTwo));
      const auto &log = logger<LogArea::ScatterAction>();
      log.debug(type_pi.name(), type_nucleus.name(), "→ ",
               type_pi.name(), produced_nucleus->name(), " at ", sqrts, " GeV, xs[mb] = ",
               xsection);
    }
  }
  return process_list;
}

}  // namespace smash
