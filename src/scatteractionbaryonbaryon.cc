
/*
 *
 *    Copyright (c) 2015-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteractionbaryonbaryon.h"

#include "include/clebschgordan.h"
#include "include/constants.h"
#include "include/cxx14compat.h"
#include "include/parametrizations.h"
#include "include/pow.h"

namespace smash {

double ScatterActionBaryonBaryon::high_energy_cross_section() const {
  const PdgCode &pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode &pdg_b = incoming_particles_[1].type().pdgcode();
  const double s = mandelstam_s();

  /* Currently all BB collisions use the nucleon-nucleon parametrizations. */
  if (pdg_a == pdg_b) {
    return pp_high_energy(s);  // pp, nn
  } else if (pdg_a.is_antiparticle_of(pdg_b)) {
    return ppbar_high_energy(s);  // ppbar, nnbar
  } else if (pdg_a.antiparticle_sign() * pdg_b.antiparticle_sign() == 1) {
    return np_high_energy(s);  // np, nbarpbar
  } else {
    return npbar_high_energy(s);  // npbar, nbarp
  }
}

double ScatterActionBaryonBaryon::string_hard_cross_section() const {
  // const PdgCode &pdg_a = incoming_particles_[0].type().pdgcode();
  // const PdgCode &pdg_b = incoming_particles_[1].type().pdgcode();
  const double s = mandelstam_s();

  /**
   * Currently nucleon-nucleon cross section is used for all case.
   * This will be changed later by applying additive quark model.
   */
  return NN_string_hard(s);
}

CollisionBranchList ScatterActionBaryonBaryon::two_to_two_cross_sections() {
  CollisionBranchList process_list;
  const ParticleType &type_a = incoming_particles_[0].type();
  const ParticleType &type_b = incoming_particles_[1].type();

  if (type_a.is_nucleon() || type_a.is_Delta() || type_b.is_nucleon() ||
      type_b.is_Delta()) {
    if (type_a.antiparticle_sign() == 1 && type_b.antiparticle_sign() == 1) {
      /* N R → N N, Δ R → N N */
      process_list = bar_bar_to_nuc_nuc(false);
    } else if (type_a.antiparticle_sign() == -1 &&
               type_b.antiparticle_sign() == -1) {
      /* N̅ R → N̅ N̅, Δ̅ R → N̅ N̅ */
      process_list = bar_bar_to_nuc_nuc(true);
    }
  }

  if ( (type_a.is_nucleon() && type_b.is_nucleus()) ||
       (type_b.is_nucleon() && type_a.is_nucleus()) ) {
    // Nd → Nd', N̅d →  N̅d', N̅d̅→ N̅d̅', Nd̅→ Nd̅'
    // and reverse (for example Nd'→ Nd).
    process_list = n_nucleus_to_n_nucleus();
  }

  return process_list;
}

CollisionBranchList ScatterActionBaryonBaryon::n_nucleus_to_n_nucleus() {
  const ParticleType &type_a = incoming_particles_[0].type();
  const ParticleType &type_b = incoming_particles_[1].type();
  const ParticleType &type_N = type_a.is_nucleon() ? type_a : type_b;
  const ParticleType &type_nucleus = type_a.is_nucleus() ? type_a : type_b;
  CollisionBranchList process_list;
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
    double matrix_element = 0.0;
    if (std::signbit(type_N.baryon_number()) ==
        std::signbit(type_nucleus.baryon_number())) {
      // Nd → Nd', N̅d̅→ N̅d̅' and reverse
      const double tmp = sqrts - type_N.min_mass_kinematic() -
                                 type_nucleus.min_mass_kinematic();
      assert(tmp >= 0.0);
      matrix_element = 79.0474 / std::pow(tmp, 0.7897) + 654.596 * tmp;
    } else {
      // N̅d →  N̅d', Nd̅→ Nd̅' and reverse
      matrix_element = 681.4;
    }
    const double spin_factor = (produced_nucleus->spin() + 1) *
                               (type_N.spin() + 1);
    // Isospin factor is always the same, so it is included into matrix element
    // Symmetry factor is always 1 here
    // Absorb (hbarc)^2/16 pi factor into matrix element
    double xsection = matrix_element * spin_factor / (s * cm_momentum());
    if (produced_nucleus->is_stable()) {
      assert(!type_nucleus.stable());
      xsection *= pCM_from_s(s, type_N.mass(), produced_nucleus->mass());
    } else {
      assert(type_nucleus.stable());
      const double resonance_integral =
          produced_nucleus->iso_multiplet()->get_integral_NR(sqrts);
      xsection *= resonance_integral;
      //std::cout << "Resonance integral " << resonance_integral
      //          << ", matrix element: " << matrix_element
      //          << ", cm_momentum: " << cm_momentum() << std::endl;
    }
    process_list.push_back(make_unique<CollisionBranch>(
          type_N, *produced_nucleus, xsection, ProcessType::TwoToTwo));
    const auto &log = logger<LogArea::ScatterAction>();
    log.debug(type_N.name(), type_nucleus.name(), "→ ",
             type_N.name(), produced_nucleus->name(),
             " at ", sqrts, " GeV, xs[mb] = ", xsection);
  }
  return process_list;
}

CollisionBranchList ScatterActionBaryonBaryon::bar_bar_to_nuc_nuc(
    const bool is_anti_particles) {
  const ParticleType &type_a = incoming_particles_[0].type();
  const ParticleType &type_b = incoming_particles_[1].type();
  CollisionBranchList process_list;

  const double s = mandelstam_s();
  const double sqrts = sqrt_s();
  /* CM momentum in final state */
  double p_cm_final = std::sqrt(s - 4. * nucleon_mass * nucleon_mass) / 2.;

  ParticleTypePtrList nuc_or_anti_nuc;
  if (is_anti_particles) {
    nuc_or_anti_nuc = ParticleType::list_anti_nucleons();
  } else {
    nuc_or_anti_nuc = ParticleType::list_nucleons();
  }

  /* Loop over all nucleon or anti-nucleon charge states. */
  for (ParticleTypePtr nuc_a : nuc_or_anti_nuc) {
    for (ParticleTypePtr nuc_b : nuc_or_anti_nuc) {
      /* Check for charge conservation. */
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
            nn_to_resonance_matrix_element(sqrts, type_a, type_b, twoI);
        if (matrix_element <= 0.) {
          continue;
        }

        /** Cross section for 2->2 resonance absorption, obtained via detailed
         * balance from the inverse reaction.
         * See eqs. (B.6), (B.9) and (181) in \iref{Buss:2011mx}.
         * There are factors for spin, isospin and symmetry involved. */
        const double spin_factor = (nuc_a->spin() + 1) * (nuc_b->spin() + 1);
        const int sym_fac_in =
            (type_a.iso_multiplet() == type_b.iso_multiplet()) ? 2 : 1;
        const int sym_fac_out =
            (nuc_a->iso_multiplet() == nuc_b->iso_multiplet()) ? 2 : 1;
        const double xsection = isospin_factor * spin_factor * sym_fac_in /
                                sym_fac_out * p_cm_final * matrix_element /
                                (s * cm_momentum());

        if (xsection > really_small) {
          process_list.push_back(make_unique<CollisionBranch>(
              *nuc_a, *nuc_b, xsection, ProcessType::TwoToTwo));
          const auto &log = logger<LogArea::ScatterAction>();
          log.debug("2->2 absorption with original particles: ", type_a,
                    type_b);
        }
      }
    }
  }
  return process_list;
}

double ScatterActionBaryonBaryon::nn_to_resonance_matrix_element(
    double sqrts, const ParticleType &type_a, const ParticleType &type_b,
    const int twoI) {
  const double m_a = type_a.mass();
  const double m_b = type_b.mass();
  const double msqr = 2. * (m_a * m_a + m_b * m_b);
  /* If the c.m. energy is larger than the sum of the pole masses of the
   * outgoing particles plus three times of the sum of the widths plus 3 GeV,
   * the collision will be neglected.*/
  const double w_a = type_a.width_at_pole();
  const double w_b = type_b.width_at_pole();
  const double uplmt = m_a + m_b + 3.0 * (w_a + w_b) + 3.0;
  if (sqrts > uplmt) {
    return 0.;
  }
  /** NN → NΔ: fit sqrt(s)-dependence to OBE model [\iref{Dmitriev:1986st}] */
  if (((type_a.is_Delta() && type_b.is_nucleon()) ||
       (type_b.is_Delta() && type_a.is_nucleon())) &&
      (type_a.antiparticle_sign() == type_b.antiparticle_sign())) {
    return 68. / std::pow(sqrts - 1.104, 1.951);
    /** All other processes use a constant matrix element,
     *  similar to \iref{Bass:1998ca}, equ. (3.35). */
  } else if (((type_a.is_Nstar() && type_b.is_nucleon()) ||
              (type_b.is_Nstar() && type_a.is_nucleon())) &&
             type_a.antiparticle_sign() == type_b.antiparticle_sign()) {
    // NN → NN*
    if (twoI == 2) {
      return 7. / msqr;
    } else if (twoI == 0) {
      const double parametrization = 14. / msqr;
      /* pn → pnη cross section is known to be larger than the corresponding
       * pp → ppη cross section by a factor of 6.5 [\iref{Calen:1998vh}].
       * Since the eta is mainly produced by an intermediate N*(1535) we
       * introduce an explicit isospin asymmetry for the production of N*(1535)
       * produced in pn vs. pp similar to [\iref{Teis:1996kx}], eq. 29. */
      if (type_a.is_Nstar1535() || type_b.is_Nstar1535()) {
        return 6.5 * parametrization;
      } else {
        return parametrization;
      }
    }
  } else if (((type_a.is_Deltastar() && type_b.is_nucleon()) ||
              (type_b.is_Deltastar() && type_a.is_nucleon())) &&
             type_a.antiparticle_sign() == type_b.antiparticle_sign()) {
    // NN → NΔ*
    return 15. / msqr;
  } else if ((type_a.is_Delta() && type_b.is_Delta()) &&
             (type_a.antiparticle_sign() == type_b.antiparticle_sign())) {
    // NN → ΔΔ
    if (twoI == 2) {
      return 45. / msqr;
    } else if (twoI == 0) {
      return 120. / msqr;
    }
  } else if (((type_a.is_Nstar() && type_b.is_Delta()) ||
              (type_b.is_Nstar() && type_a.is_Delta())) &&
             type_a.antiparticle_sign() == type_b.antiparticle_sign()) {
    // NN → ΔN*
    return 7. / msqr;
  } else if (((type_a.is_Deltastar() && type_b.is_Delta()) ||
              (type_b.is_Deltastar() && type_a.is_Delta())) &&
             type_a.antiparticle_sign() == type_b.antiparticle_sign()) {
    // NN → ΔΔ*
    if (twoI == 2) {
      return 15. / msqr;
    } else if (twoI == 0) {
      return 25. / msqr;
    }
  } else if ((type_a.is_nucleus() && type_b.pdgcode().is_pion()) ||
             (type_b.is_nucleus() && type_a.pdgcode().is_pion())) {
    // This parametrization is the result of fitting d+pi->NN cross-section.
    // Already Breit-Wigner-like part provides a good fit, exponential fixes
    // behaviour around the treshold.
    return 0.055 / (pow_int(sqrts - 2.145, 2) + pow_int(0.065, 2)) *
           (1.0 - std::exp(-(sqrts - 2.0) * 20.0));
  }
  // all cases not listed: zero!
  return 0.;
}

void ScatterActionBaryonBaryon::format_debug_output(std::ostream &out) const {
  out << "Baryon-Baryon ";
  ScatterAction::format_debug_output(out);
}

}  // namespace smash
