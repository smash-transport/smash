/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <memory>

#include "include/action.h"
#include "include/angles.h"
#include "include/constants.h"
#include "include/cxx14compat.h"
#include "include/distributions.h"
#include "include/kinematics.h"
#include "include/logging.h"
#include "include/parametrizations.h"
#include "include/pdgcode.h"
#include "include/random.h"
#include "include/resonances.h"

namespace Smash {

ScatterAction::ScatterAction(const ParticleData &in_part_a,
                             const ParticleData &in_part_b,
                             float time_of_execution)
    : Action({in_part_a, in_part_b}, time_of_execution),
      total_cross_section_(0.) {}

void ScatterAction::add_collision(CollisionBranchPtr p) {
  add_process<CollisionBranch>(p, collision_channels_, total_cross_section_);
}

void ScatterAction::add_collisions(CollisionBranchList pv) {
  add_processes<CollisionBranch>(std::move(pv), collision_channels_,
      total_cross_section_);
}

void ScatterAction::generate_final_state() {
  const auto &log = logger<LogArea::ScatterAction>();

  log.debug("Incoming particles: ", incoming_particles_);

  /* Decide for a particular final state. */
  const CollisionBranch* proc = choose_channel <CollisionBranch>(
      collision_channels_, total_cross_section_);
  process_type_ = proc->get_type();
  outgoing_particles_ = proc->particle_list();

  log.debug("Chosen channel: ", process_type_, outgoing_particles_);

  /* The production point of the new particles.  */
  FourVector middle_point = get_interaction_point();

  switch (process_type_) {
    case ProcessType::Elastic:
      /* 2->2 elastic scattering */
      log.debug("Process: Elastic collision.", process_type_);
      momenta_exchange();
      break;
    case ProcessType::TwoToOne:
      /* resonance formation */
      log.debug("Process: Resonance formation.", process_type_);
      /* processes computed in the center of momenta */
      resonance_formation();
      break;
    case ProcessType::TwoToTwo:
      /* 2->2 inelastic scattering */
      log.debug("Process: Inelastic scattering.", process_type_);
      /* Sample the particle momenta in CM system. */
      sample_cms_momenta();
      break;
    case ProcessType::String:
      /* string excitation */
      log.debug("Process: String Excitation.");
      /// string_excitation(incoming_particles_, outgoing_particles_);
      break;
    case ProcessType::None:
      log.debug("ProcessType None should not have been selected");
      break;
    case ProcessType::Decay:
      log.debug("ProcessType Decay should have been handled as DecayAction");
      break;
    default:
      throw InvalidScatterAction(
        "ScatterAction::perform: Unknown Process Type. "
        "ProcessType " + std::to_string(static_cast<int>(process_type_)) +
        " was requested. (PDGcode1=" + incoming_particles_[0].pdgcode().string()
        + ", PDGcode2=" + incoming_particles_[1].pdgcode().string()
        + ")");
  }

  /* Set positions & boost to computational frame. */
  for (ParticleData &new_particle : outgoing_particles_) {
    new_particle.set_4position(middle_point);
    new_particle.boost_momentum(-beta_cm());
  }
}


float ScatterAction::raw_weight_value() const {
  return total_cross_section_;
}


ThreeVector ScatterAction::beta_cm() const {
  FourVector mom = incoming_particles_[0].momentum() +
                   incoming_particles_[1].momentum();
  return mom.threevec() / mom.x0();
}


double ScatterAction::mandelstam_s() const {
  return (incoming_particles_[0].momentum() +
          incoming_particles_[1].momentum()).sqr();
}

double ScatterAction::sqrt_s() const {
  return (incoming_particles_[0].momentum() +
          incoming_particles_[1].momentum()).abs();
}

double ScatterAction::cm_momentum() const {
  const double m1 = incoming_particles_[0].effective_mass();
  const double m2 = incoming_particles_[1].effective_mass();
  return pCM(sqrt_s(), m1, m2);
}

double ScatterAction::cm_momentum_squared() const {
  const double m1 = incoming_particles_[0].effective_mass();
  const double m2 = incoming_particles_[1].effective_mass();
  return pCM_sqr(sqrt_s(), m1, m2);
}


double ScatterAction::particle_distance() const {
  const auto &log = logger<LogArea::ScatterAction>();
  // local copy of particles (since we need to boost them)
  ParticleData p_a = incoming_particles_[0];
  ParticleData p_b = incoming_particles_[1];
  /* Boost particles to center-of-momentum frame. */
  ThreeVector velocity = beta_cm();
  p_a.boost(velocity);
  p_b.boost(velocity);
  ThreeVector pos_diff = p_a.position().threevec() - p_b.position().threevec();
  ThreeVector mom_diff = p_a.momentum().threevec() - p_b.momentum().threevec();
  log.debug("Particle ", incoming_particles_, " position difference [fm]: ",
            pos_diff, ", momentum difference [GeV]: ", mom_diff);
  /* Zero momentum leads to infite distance. */
  if (std::abs(mom_diff.sqr()) < really_small)
    return  pos_diff.sqr();

  /* UrQMD squared distance criteria:
   * arXiv:nucl-th/9803035 (3.27): in center of momemtum frame
   * position of particle a: x_a
   * position of particle b: x_b
   * velocity of particle a: v_a
   * velocity of particle b: v_b
   * d^2_{coll} = (x_a - x_b)^2 - ((x_a - x_a) . (v_a - v_b))^2 / (v_a - v_b)^2
   */
  return pos_diff.sqr()
         - (pos_diff * mom_diff) * (pos_diff * mom_diff) / mom_diff.sqr();
}


CollisionBranchPtr ScatterAction::elastic_cross_section(float elast_par) {
  return make_unique<CollisionBranch>(incoming_particles_[0].type(),
                                      incoming_particles_[1].type(),
                                      elast_par, ProcessType::Elastic);
}

CollisionBranchPtr ScatterAction::string_excitation_cross_section() {
  /* Calculate string-excitation cross section:
   * Parametrized total minus all other present channels. */
  /* TODO(weil): This is currently set to zero,
   * since Pythia is not yet implemented. */
  float sig_string = 0.f;
  // = std::max(0.f, total_cross_section() - total_weight_);

  return make_unique<CollisionBranch>(sig_string, ProcessType::String);
}


double ScatterAction::two_to_one_formation(const ParticleType &type_resonance,
                                           double s, double cm_momentum_sqr) {
  const ParticleType &type_particle_a = incoming_particles_[0].type();
  const ParticleType &type_particle_b = incoming_particles_[1].type();
  /* Check for charge conservation. */
  if (type_resonance.charge() != type_particle_a.charge()
                               + type_particle_b.charge()) {
    return 0.;
  }

  /* Check for baryon-number conservation. */
  if (type_resonance.baryon_number() != type_particle_a.baryon_number()
                                      + type_particle_b.baryon_number()) {
    return 0.;
  }

  /* Calculate partial in-width. */
  double srts = std::sqrt(s);
  float partial_width = type_resonance.get_partial_in_width(srts,
                                incoming_particles_[0], incoming_particles_[1]);
  if (partial_width <= 0.) {
    return 0.;
  }

  /* Calculate spin factor */
  const double spinfactor = (type_resonance.spin() + 1)
    / ((type_particle_a.spin() + 1) * (type_particle_b.spin() + 1));
  const int sym_factor = (type_particle_a.pdgcode() ==
                          type_particle_b.pdgcode()) ? 2 : 1;
  float resonance_width = type_resonance.total_width(srts);
  float resonance_mass = type_resonance.mass();
  /* Calculate resonance production cross section
   * using the Breit-Wigner distribution as probability amplitude.
   * See Eq. (176) in Buss et al., Physics Reports 512, 1 (2012). */
  return spinfactor * sym_factor * 4.0 * M_PI / cm_momentum_sqr
         * breit_wigner(s, resonance_mass, resonance_width)
         * partial_width/resonance_width
         * hbarc * hbarc / fm2_mb;
}


CollisionBranchList ScatterAction::resonance_cross_sections() {
  const auto &log = logger<LogArea::ScatterAction>();
  CollisionBranchList resonance_process_list;
  const ParticleType &type_particle_a = incoming_particles_[0].type();
  const ParticleType &type_particle_b = incoming_particles_[1].type();

  const double s = mandelstam_s();
  const double p_cm_sqr = cm_momentum_squared();

  /* Find all the possible resonances */
  for (const ParticleType &type_resonance : ParticleType::list_all()) {
    /* Not a resonance, go to next type of particle */
    if (type_resonance.is_stable()) {
      continue;
    }

    /* Same resonance as in the beginning, ignore */
    if ((!type_particle_a.is_stable()
         && type_resonance.pdgcode() == type_particle_a.pdgcode())
        || (!type_particle_b.is_stable()
            && type_resonance.pdgcode() == type_particle_b.pdgcode())) {
      continue;
    }

    float resonance_xsection = two_to_one_formation(type_resonance,
                                                    s, p_cm_sqr);

    /* If cross section is non-negligible, add resonance to the list */
    if (resonance_xsection > really_small) {
      resonance_process_list.push_back(make_unique<CollisionBranch>
                                       (type_resonance, resonance_xsection,
                                        ProcessType::TwoToOne));
      log.debug("Found resonance: ", type_resonance);
      log.debug("2->1 with original particles: ", type_particle_a,
                type_particle_b);
    }
  }
  return resonance_process_list;
}


void ScatterAction::momenta_exchange() {
  const auto &log = logger<LogArea::ScatterAction>();
  outgoing_particles_[0] = incoming_particles_[0];
  outgoing_particles_[1] = incoming_particles_[1];

  /* Determine absolute momentum in center-of-mass frame. */
  const double momentum_radial = cm_momentum();

  /* Particle exchange momenta and scatter to random direction.
   * TODO: Angles should be sampled from differential cross section
   * of this process. */
  Angles phitheta;
  phitheta.distribute_isotropically();
  log.debug("Random momentum: ", momentum_radial, " ", phitheta);

  /* Set 4-momentum: Masses stay the same, 3-momentum changes. */
  outgoing_particles_[0].set_4momentum(outgoing_particles_[0].effective_mass(),
                                       phitheta.threevec() * momentum_radial);
  outgoing_particles_[1].set_4momentum(outgoing_particles_[1].effective_mass(),
                                       -phitheta.threevec() * momentum_radial);

  /* debug output */
  log.debug("exchanged momenta a", outgoing_particles_[0].momentum());
  log.debug("exchanged momenta b", outgoing_particles_[1].momentum());
}


void ScatterAction::resonance_formation() {
  const auto &log = logger<LogArea::ScatterAction>();

  if (outgoing_particles_.size() != 1) {
    std::string s = "resonance_formation: "
                    "Incorrect number of particles in final state: ";
    s += std::to_string(outgoing_particles_.size()) + " (";
    s += incoming_particles_[0].pdgcode().string() + " + ";
    s += incoming_particles_[1].pdgcode().string() + ")";
    throw InvalidResonanceFormation(s);
  }

  /* 1 particle in final state: Center-of-momentum frame of initial particles
   * is the rest frame of the resonance.  */
  outgoing_particles_[0].set_4momentum(FourVector(sqrt_s(), 0., 0., 0.));

  log.debug("Momentum of the new particle: ",
            outgoing_particles_[0].momentum());
}


/***** ScatterActionBaryonBaryon **********************************************/

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


CollisionBranchPtr ScatterActionBaryonBaryon::elastic_cross_section(
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
      return make_unique<CollisionBranch>(incoming_particles_[0].type(),
                                          incoming_particles_[1].type(),
                                          sig_el, ProcessType::Elastic);
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


CollisionBranchList ScatterActionBaryonBaryon::two_to_two_cross_sections() {
  CollisionBranchList process_list;
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


CollisionBranchList ScatterActionBaryonBaryon::nuc_nuc_to_nuc_res(
                            const ParticleType &type_particle_a,
                            const ParticleType &type_particle_b) {
  const auto &log = logger<LogArea::ScatterAction>();
  CollisionBranchList process_list;
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
      int I_min = std::max(
          std::abs(type_resonance->isospin() - second_type->isospin()),
          std::abs(type_particle_a.isospin() - type_particle_b.isospin()));
      I_min = std::max(I_min, std::abs(I_z));

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
                                ProcessType::TwoToTwo));
        log.debug("Found 2->2 creation process for resonance ",
                  *type_resonance);
        log.debug("2->2 with original particles: ",
                  type_particle_a, type_particle_b);
      }
    }
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
  double p_cm_final = std::sqrt(s - 4. * type_nucleon->mass_sqr()) / 2.;

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


void ScatterAction::format_debug_output(std::ostream &out) const {
  out << "Scatter of " << incoming_particles_;
  if (outgoing_particles_.empty()) {
    out << " (not performed)";
  } else {
    out << " to " << outgoing_particles_;
  }
}
void ScatterActionMesonMeson::format_debug_output(std::ostream &out) const {
  out << " Meson-Meson  ";
  ScatterAction::format_debug_output(out);
}
void ScatterActionBaryonMeson::format_debug_output(std::ostream &out) const {
  out << "Baryon-Meson  ";
  ScatterAction::format_debug_output(out);
}
void ScatterActionBaryonBaryon::format_debug_output(std::ostream &out) const {
  out << "Baryon-Baryon ";
  ScatterAction::format_debug_output(out);
}


}  // namespace Smash
