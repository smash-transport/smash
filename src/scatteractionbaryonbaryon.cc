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
  /* CM momentum in final state */
  double p_cm_final = sqrt(s - 4.*type_nucleon->mass_sqr())/2.;

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
      int I_min = std::max(abs(type_resonance->isospin()
                               - type_nucleon->isospin()),
                          abs(nuc_a->isospin() - nuc_b->isospin()));
      I_min = std::max(I_min, abs(I_z));

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

      /* Calculate matrix element. */
      const float matrix_element = nn_to_resonance_matrix_element(s,
                                                *type_resonance, *type_nucleon);
      if (matrix_element <= 0.) {
        continue;
      }

      /* Cross section for 2->2 resonance absorption, obtained via detailed
       * balance from the inverse reaction.
       * See eqs. (B.6), (B.9) and (181) in the GiBUU review paper.
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


void ScatterActionBaryonBaryon::format_debug_output(std::ostream &out) const {
  out << "Baryon-Baryon ";
  ScatterAction::format_debug_output(out);
}


}  // namespace Smash
