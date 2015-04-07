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


CollisionBranch* ScatterActionBaryonBaryon::elastic_cross_section(
                                                              float elast_par) {
  const PdgCode &pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode &pdg_b = incoming_particles_[1].type().pdgcode();

  const double s = mandelstam_s();

  if ((pdg_a.iso_multiplet() == 0x1112) &&
      (pdg_b.iso_multiplet() == 0x1112)) {
    /* Nucleon-Nucleon scattering: use parametrized cross sections. */
    float sig_el;
    if (pdg_a == pdg_b) {                          /* pp */
      sig_el = pp_elastic(s);
    } else if (pdg_a.is_antiparticle_of(pdg_b)) {  /* ppbar */
      sig_el = ppbar_elastic(s);
    } else {                                     /* np */
      sig_el = np_elastic(s);
    }
    if (sig_el > 0.) {
      return new CollisionBranch(incoming_particles_[0].type(),
                                 incoming_particles_[1].type(),
                                 sig_el, ProcessBranch::Elastic);
    } else {
      std::stringstream ss;
      ss << "problem in CrossSections::elastic: " << pdg_a.string().c_str()
        << " " << pdg_b.string().c_str() << " " << pdg_a.spin() << " "
        << pdg_b.spin() << " " << sig_el << " " << s;
      throw std::runtime_error(ss.str());
    }
  } else {
    /* Default: Fall back to parent routine. */
    return ScatterAction::elastic_cross_section(elast_par);
  }
}


ProcessBranchList ScatterActionBaryonBaryon::two_to_two_cross_sections() {
  ProcessBranchList process_list;
  const ParticleType &type_particle_a = incoming_particles_[0].type();
  const ParticleType &type_particle_b = incoming_particles_[1].type();

  if (type_particle_a.pdgcode().iso_multiplet() == 0x1112 &&
      type_particle_b.pdgcode().iso_multiplet() == 0x1112) {
    /* Nucleon+Nucleon: find all resonance production channels */
      process_list = nuc_nuc_to_nuc_res(type_particle_a, type_particle_b);
  } else if (type_particle_a.pdgcode().iso_multiplet() == 0x1112 ||
             type_particle_b.pdgcode().iso_multiplet() == 0x1112) {
    /* Nucleon+Resonance: absorption */
    process_list = nuc_res_to_nuc_nuc(type_particle_a, type_particle_b);
  }

  return process_list;
}


/**
 * Scattering matrix amplitude squared for \f$NN \rightarrow NR\f$ processes,
 * where R is a baryon resonance (Delta, N*, Delta*).
 *
 * \param[in] mandelstam_s Mandelstam-s, i.e. collision CMS energy squared.
 * \param[in] type_final_a Type information for the first final state particle.
 * \param[in] type_final_b Type information for the second final state particle.
 *
 * \return Matrix amplitude squared \f$|\mathcal{M}(\sqrt{s})|^2/16\pi\f$.
 */
static float nn_to_resonance_matrix_element(const double mandelstam_s,
  const ParticleType &type_final_a, const ParticleType &type_final_b) {
  PdgCode delta = PdgCode("2224");
  if (type_final_a.pdgcode().iso_multiplet()
      != type_final_b.pdgcode().iso_multiplet()) {
    /* N + N -> N + Delta: fit to Dmitriev OBE model,
     * Nucl. Phys. A 459, 503 (1986) */
    if (type_final_a.pdgcode().iso_multiplet() == delta.iso_multiplet()
        || type_final_b.pdgcode().iso_multiplet() == delta.iso_multiplet()) {
      return 459. / std::pow(std::sqrt(mandelstam_s) - 1.104, 1.951);
    } else {
      return 0.0;
    }
  } else {
    return 0.0;
  }
}


ProcessBranchList ScatterActionBaryonBaryon::nuc_nuc_to_nuc_res(
                            const ParticleType &type_particle_a,
                            const ParticleType &type_particle_b) {
  const auto &log = logger<LogArea::ScatterAction>();
  ProcessBranchList process_list;
  const double s = mandelstam_s();

  /* Loop over all baryon resonances. */
  for (ParticleTypePtr type_resonance :
       ParticleType::list_baryon_resonances()) {
    /* Loop over second particle (nucleon). */
    for (ParticleTypePtr second_type : ParticleType::list_nucleons()) {
      /* Check for charge conservation. */
      if (type_resonance->charge() + second_type->charge() !=
          type_particle_a.charge() + type_particle_b.charge()) {
        continue;
      }

      int I_z = type_resonance->isospin3() + second_type->isospin3();

      /* Compute total isospin range with given initial and final particles. */
      int I_max =
          std::min(type_resonance->isospin() + second_type->isospin(),
                   type_particle_a.isospin() + type_particle_b.isospin());
      int I_min =
          std::max(abs(type_resonance->isospin() - second_type->isospin()),
                   abs(type_particle_a.isospin() - type_particle_b.isospin()));
      I_min = std::max(I_min, abs(I_z));

      /* Loop over total isospin in allowed range.
      * Use decrement of 2, since isospin is multiplied by 2. */
      double isospin_factor = 0.;
      for (int I_tot = I_max; I_tot >= I_min; I_tot -= 2) {
        isospin_factor = isospin_factor +
                         isospin_clebsch_gordan(type_particle_a,
                                                type_particle_b, I_tot, I_z) *
                             isospin_clebsch_gordan(*type_resonance,
                                                    *second_type, I_tot, I_z);
      }

      /* If Clebsch-Gordan coefficient is zero, don't bother with the rest. */
      if (std::abs(isospin_factor) < really_small) {
        continue;
      }

      /* Integration limits. */
      double lower_limit = type_resonance->minimum_mass();
      double upper_limit = std::sqrt(s) - second_type->mass();
      /* Check the available energy (requiring it to be a little above the
      * threshold, because the integration will not work if it's too close). */
      if (upper_limit - lower_limit < 1E-3) {
        continue;
      }

      /* Calculate matrix element. */
      const float matrix_element =
          nn_to_resonance_matrix_element(s, *type_resonance, *second_type);
      if (matrix_element <= 0.) {
        continue;
      }

      /* Calculate resonance production cross section
       * using the Breit-Wigner distribution as probability amplitude.
       * Integrate over the allowed resonance mass range. */
      IntegrandParameters params = {type_resonance, second_type->mass(), s};
      log.debug("Process: ", type_particle_a, type_particle_b, " -> ",
                *second_type, *type_resonance);
      log.debug("Limits: ", lower_limit, " ", upper_limit);
      double resonance_integral, integral_error;
      quadrature_1d(&spectral_function_integrand, &params,
                    lower_limit, upper_limit,
                    &resonance_integral, &integral_error);
      log.debug("Integral value: ", resonance_integral,
                " Error: ", integral_error);

      /* Cross section for 2->2 process with one resonance in final state.
       * Based on Eq. (46) in PhD thesis of J. Weil
       * (https://gibuu.hepforge.org/trac/chrome/site/files/phd/weil.pdf) */
      float xsection = isospin_factor * isospin_factor * matrix_element
                     * resonance_integral / (s * cm_momentum());

      if (xsection > really_small) {
        process_list.push_back(make_unique<CollisionBranch>
                               (*type_resonance, *second_type, xsection,
                                ProcessBranch::TwoToTwo));
        log.debug("Found 2->2 creation process for resonance ",
                  *type_resonance);
        log.debug("2->2 with original particles: ",
                  type_particle_a, type_particle_b);
      }
    }
  }
  return process_list;
}


ProcessBranchList ScatterActionBaryonBaryon::nuc_res_to_nuc_nuc(
                            const ParticleType &type_particle_a,
                            const ParticleType &type_particle_b) {
  ParticleTypePtr type_resonance, type_nucleon;
  ProcessBranchList process_list;

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
                                ProcessBranch::TwoToTwo));
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
