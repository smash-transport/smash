/*
 *
 *    Copyright (c) 2015
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

namespace Smash {


float ScatterActionBaryonBaryon::total_cross_section() const {
  const PdgCode &pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode &pdg_b = incoming_particles_[1].type().pdgcode();
  const double s = mandelstam_s();

  /* Currently all BB collisions use the nucleon-nucleon parametrizations. */
  if (pdg_a == pdg_b) {
    return pp_total(s);     // pp, nn
  } else if (pdg_a.is_antiparticle_of(pdg_b)) {
    return ppbar_total(s);  // NNbar
  } else {
    return np_total(s);     // np
  }
}


CollisionBranchList ScatterActionBaryonBaryon::two_to_two_cross_sections() {
  CollisionBranchList process_list;
  const ParticleType &type_a = incoming_particles_[0].type();
  const ParticleType &type_b = incoming_particles_[1].type();

  if (type_a.pdgcode().is_nucleon() || type_a.pdgcode().is_Delta() ||
      type_b.pdgcode().is_nucleon() || type_b.pdgcode().is_Delta()) {
    /* N R -> N N, Delta R -> N N */
    process_list = bar_bar_to_nuc_nuc(type_a, type_b);
  }

  return process_list;
}


CollisionBranchList ScatterActionBaryonBaryon::bar_bar_to_nuc_nuc(
                            const ParticleType &type_a,
                            const ParticleType &type_b) {
  CollisionBranchList process_list;

  const double s = mandelstam_s();
  const double sqrts = sqrt_s();
  /* CM momentum in final state */
  double p_cm_final = std::sqrt(s - 4.*nucleon_mass*nucleon_mass)/2.;

  /* Loop over all nucleon charge states. */
  for (ParticleTypePtr nuc_a : ParticleType::list_nucleons()) {
    for (ParticleTypePtr nuc_b : ParticleType::list_nucleons()) {
      /* Check for charge conservation. */
      if (type_a.charge() + type_b.charge()
          != nuc_a->charge() + nuc_b->charge()) {
        continue;
      }

      double isospin_factor = isospin_clebsch_gordan(type_a, type_b,
                                                     *nuc_a, *nuc_b);
      /* If Clebsch-Gordan coefficient is zero, don't bother with the rest */
      if (std::abs(isospin_factor) < really_small) {
        continue;
      }

      /* Calculate matrix element for inverse process. */
      const float matrix_element =
          nn_to_resonance_matrix_element(sqrts, type_a, type_b);
      if (matrix_element <= 0.) {
        continue;
      }

      /** Cross section for 2->2 resonance absorption, obtained via detailed
       * balance from the inverse reaction.
       * See eqs. (B.6), (B.9) and (181) in \iref{Buss:2011mx}.
       * There is a symmetry factor 1/2 and a spin factor 4/(2S_a+1)/(2S_b+1)
       * involved. */
      float xsection = isospin_factor * isospin_factor
                     * p_cm_final * matrix_element * 2.
                     / ((type_a.spin()+1)*(type_b.spin()+1) * s*cm_momentum());

      if (xsection > really_small) {
        process_list.push_back(make_unique<CollisionBranch>
                               (*nuc_a, *nuc_b, xsection,
                                ProcessType::TwoToTwo));
        const auto &log = logger<LogArea::ScatterAction>();
        log.debug("2->2 absorption with original particles: ",
                  type_a, type_b);
      }
    }
  }
  return process_list;
}


float ScatterActionBaryonBaryon::nn_to_resonance_matrix_element(double sqrts,
      const ParticleType &type_a, const ParticleType &type_b) const {
  const float spin_factor = (type_a.spin()+1) * (type_b.spin()+1);
  const float m_plus = type_a.mass() + type_b.mass();
  const float m_minus = type_a.mass() - type_b.mass();

  const PdgCode pdg_a = type_a.pdgcode();
  const PdgCode pdg_b = type_b.pdgcode();

  if (pdg_a.is_Delta() && pdg_b.is_nucleon()) {
    /** \f$ NN \rightarrow N\Delta \f$:
      * fit sqrt(s)-dependence to OBE model [\iref{Dmitriev:1986st}] */
    return 57.375 * spin_factor / std::pow(sqrts - 1.104, 1.951);
  } else if (pdg_a.is_Nstar() && pdg_b.is_nucleon()) {
    /** \f$ NN \rightarrow NN^* \f$:
      * constant matrix element, cf. \iref{Bass:1998ca}, equ. (3.35). */
    return 15. * spin_factor / (m_plus * m_plus + m_minus * m_minus);
  } else if (pdg_a.is_Deltastar() && pdg_b.is_nucleon()) {
    /** \f$ NN \rightarrow N\Delta^* \f$:
      * constant matrix element, cf. \iref{Bass:1998ca}, equ. (3.35). */
    return 25. * spin_factor / (m_plus * m_plus + m_minus * m_minus);
  } else if (pdg_a.is_Delta() && pdg_b.is_Delta()) {
    /** \f$ NN \rightarrow \Delta\Delta \f$:
      * constant matrix element, cf. \iref{Bass:1998ca}, equ. (3.35). */
    return 20. * spin_factor / (m_plus * m_plus + m_minus * m_minus);
  } else if (pdg_a.is_Nstar() && pdg_b.is_Delta()) {
    /** \f$ NN \rightarrow \Delta N^* \f$:
      * constant matrix element, cf. \iref{Bass:1998ca}, equ. (3.35). */
    return 10. * spin_factor / (m_plus * m_plus + m_minus * m_minus);
  } else if (pdg_a.is_Deltastar() && pdg_b.is_Delta()) {
    /** \f$ NN \rightarrow \Delta\Delta^* \f$:
      * constant matrix element, cf. \iref{Bass:1998ca}, equ. (3.35). */
    return 15. * spin_factor / (m_plus * m_plus + m_minus * m_minus);
  } else {
    return 0.0;
  }
}


void ScatterActionBaryonBaryon::format_debug_output(std::ostream &out) const {
  out << "Baryon-Baryon ";
  ScatterAction::format_debug_output(out);
}


}  // namespace Smash
