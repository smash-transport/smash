/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteractionbaryonbaryon.h"

#include "include/cxx14compat.h"
#include "include/parametrizations.h"
#include "include/resonances.h"

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
  const ParticleType &type_particle_a = incoming_particles_[0].type();
  const ParticleType &type_particle_b = incoming_particles_[1].type();

  if (type_particle_a.pdgcode().iso_multiplet() == 0x1112 ||
      type_particle_b.pdgcode().iso_multiplet() == 0x1112) {
    /* Nucleon+Resonance: absorption */
    process_list = nuc_res_to_nuc_nuc(type_particle_a, type_particle_b);
  }

  return process_list;
}


CollisionBranchList ScatterActionBaryonBaryon::nuc_res_to_nuc_nuc(
                            const ParticleType &type_particle_a,
                            const ParticleType &type_particle_b) {
  ParticleTypePtr type_resonance, type_nucleon;
  CollisionBranchList process_list;

  if (type_particle_a.pdgcode().iso_multiplet() == 0x1112) {
    type_nucleon = &type_particle_a;
    type_resonance = &type_particle_b;
  } else if (type_particle_b.pdgcode().iso_multiplet() == 0x1112) {
    type_nucleon = &type_particle_b;
    type_resonance = &type_particle_a;
  } else {
    throw std::runtime_error("Error: no nucleon found in nuc_res_to_nuc_nuc!");
  }

  const double s = mandelstam_s();
  const double srts = sqrt_s();
  /* CM momentum in final state */
  double p_cm_final = std::sqrt(s - 4.*type_nucleon->mass_sqr())/2.;

  /* Loop over all nucleon charge states. */
  for (ParticleTypePtr nuc_a : ParticleType::list_nucleons()) {
    for (ParticleTypePtr nuc_b : ParticleType::list_nucleons()) {
      /* Check for charge conservation. */
      if (type_resonance->charge() + type_nucleon->charge()
          != nuc_a->charge() + nuc_b->charge()) {
        continue;
      }

      int I_z = type_resonance->isospin3() + type_nucleon->isospin3();

      /* Compute total isospin range with given initial and final particles. */
      int I_max = std::min(type_resonance->isospin() + type_nucleon->isospin(),
                          nuc_a->isospin() + nuc_b->isospin());
      int I_min = std::max(
          std::abs(type_resonance->isospin() - type_nucleon->isospin()),
          std::abs(nuc_a->isospin() - nuc_b->isospin()));
      I_min = std::max(I_min, std::abs(I_z));

      /* Loop over total isospin in allowed range.
      * Use decrement of 2, since isospin is multiplied by 2. */
      double isospin_factor = 0.;
      for (int I_tot = I_max; I_tot >= I_min; I_tot -= 2) {
        isospin_factor = isospin_factor +
            isospin_clebsch_gordan(*nuc_a, *nuc_b, I_tot, I_z)
          * isospin_clebsch_gordan(*type_resonance, *type_nucleon, I_tot, I_z);
      }

      /* If Clebsch-Gordan coefficient is zero, don't bother with the rest */
      if (std::abs(isospin_factor) < really_small) {
        continue;
      }

      /* Calculate matrix element for inverse process. */
      const float matrix_element =
          nn_to_resonance_matrix_element(srts, *type_resonance, *type_nucleon);
      if (matrix_element <= 0.) {
        continue;
      }

      /** Cross section for 2->2 resonance absorption, obtained via detailed
       * balance from the inverse reaction.
       * See eqs. (B.6), (B.9) and (181) in \iref{Buss:2011mx}.
       * There is a symmetry factor 1/2 and a spin factor 2/(2S+1) involved,
       * which combine to 1/(2S+1). */
      float xsection = isospin_factor * isospin_factor
                     * p_cm_final * matrix_element
                     / ((type_resonance->spin()+1) * s * cm_momentum());

      if (xsection > really_small) {
        process_list.push_back(make_unique<CollisionBranch>
                               (*nuc_a, *nuc_b, xsection,
                                ProcessType::TwoToTwo));
        const auto &log = logger<LogArea::ScatterAction>();
        log.debug("Found 2->2 absoption process for resonance ",
                  *type_resonance);
        log.debug("2->2 with original particles: ",
                  type_particle_a, type_particle_b);
      }
    }
  }
  return process_list;
}


float ScatterActionBaryonBaryon::nn_to_resonance_matrix_element(
      const double srts,
      const ParticleType &type_a, const ParticleType &type_b) const {
  if (type_a.pdgcode().iso_multiplet() == type_b.pdgcode().iso_multiplet()) {
    return 0.;
  }

  int delta = PdgCode("2224").iso_multiplet();
  float spin_factor = (type_a.spin()+1) * (type_b.spin()+1);
  float m_plus = type_a.mass() + type_b.mass();
  float m_minus = type_a.mass() - type_b.mass();

  if (type_a.pdgcode().iso_multiplet() == delta
      || type_b.pdgcode().iso_multiplet() == delta) {
    /** \f$ NN \rightarrow N\Delta \f$:
      * fit sqrt(s)-dependence to OBE model [\iref{Dmitriev:1986st}] */
    return 57.375 * spin_factor / std::pow(srts - 1.104, 1.951);
  } else if (type_a.isospin() == 1 && type_b.isospin() == 1) {
    /** \f$ NN \rightarrow NN^* \f$:
      * constant matrix element, cf. \iref{Bass:1998ca}, equ. (3.35). */
    return 30. * spin_factor / (m_plus * m_plus + m_minus * m_minus);
  } else if (type_a.isospin() == 3 || type_b.isospin() == 3) {
    /** \f$ NN \rightarrow N\Delta^* \f$:
      * constant matrix element, cf. [\iref{Bass:1998ca}, equ. (3.35)]. */
    return 30. * spin_factor / (m_plus * m_plus + m_minus * m_minus);
  } else {
    return 0.0;
  }
}


void ScatterActionBaryonBaryon::format_debug_output(std::ostream &out) const {
  out << "Baryon-Baryon ";
  ScatterAction::format_debug_output(out);
}


}  // namespace Smash
