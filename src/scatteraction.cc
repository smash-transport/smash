/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/action.h"

#include "include/constants.h"
#include "include/outputroutines.h"
#include "include/pdgcode.h"
#include "include/random.h"
#include "include/resonances.h"
#include "include/angles.h"
#include "include/parametrizations.h"

namespace Smash {

ScatterAction::ScatterAction(const ParticleData &in_part1,
                             const ParticleData &in_part2,
                             float time_of_execution)
    : Action({in_part1, in_part2}, time_of_execution) {}


void ScatterAction::perform(Particles *particles, size_t &id_process) {
  /* Relevant particle IDs for the collision. */
  int id1 = incoming_particles_[0].id();
  int id2 = incoming_particles_[1].id();

  printd("Process %zu particle %s<->%s colliding %d<->%d time %g\n",
         id_process, incoming_particles_[0].type().name().c_str(),
         incoming_particles_[1].type().name().c_str(), id1, id2,
         incoming_particles_[0].position().x0());
  printd_momenta("particle 1 momenta before", incoming_particles_[0]);
  printd_momenta("particle 2 momenta before", incoming_particles_[1]);

  /* Decide for a particular final state. */
  outgoing_particles_ = choose_channel();

  if (is_elastic()) {
    /* 2->2 elastic scattering */
    printd("Process: Elastic collision.\n");

    momenta_exchange();

    // store the process id in the Particle data
    outgoing_particles_[0].set_id_process(id_process);
    outgoing_particles_[1].set_id_process(id_process);

    particles->data(id1) = outgoing_particles_[0];
    particles->data(id2) = outgoing_particles_[1];
  } else {
    /* resonance formation */
    printd("Process: Resonance formation. ");

    /* The starting point of resonance is between the two initial particles:
     * x_middle = (x_a + x_b) / 2   */
    FourVector middle_point = (incoming_particles_[0].position()
                               + incoming_particles_[1].position()) / 2.;

    /* processes computed in the center of momenta */
    resonance_formation();

    /* Set positions & boost to computational frame. */
    for (ParticleData &new_particle : outgoing_particles_) {
      new_particle.set_4position(middle_point);

      new_particle.set_4momentum(
          new_particle.momentum().LorentzBoost(-beta_cm()));

      // store the process id in the Particle data
      new_particle.set_id_process(id_process);

      printd("Resonance %s with ID %i \n",
             new_particle.type().name().c_str(), new_particle.id());
      printd_momenta("momentum in comp frame", new_particle);
      printd_position("position in comp frame", new_particle);

      new_particle.set_id(particles->add_data(new_particle));
    }

    /* Remove the initial particles */
    particles->remove(id1);
    particles->remove(id2);

    printd("Particle map has now %zu elements. \n", particles->size());
  }

  check_conservation(id_process);

  id_process++;
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


double ScatterAction::cm_momentum_squared() const {
  return (incoming_particles_[0].momentum().Dot(incoming_particles_[1].momentum())
       * incoming_particles_[0].momentum().Dot(incoming_particles_[1].momentum())
       - incoming_particles_[0].type().mass_sqr()
       * incoming_particles_[1].type().mass_sqr()) / mandelstam_s();
}


double ScatterAction::particle_distance() const {
  // local copy of particles (since we need to boost them)
  ParticleData p1 = incoming_particles_[0];
  ParticleData p2 = incoming_particles_[1];
  /* Boost particles to center-of-momentum frame. */
  ThreeVector velocity = beta_cm();
  p1.boost(velocity);
  p2.boost(velocity);
  ThreeVector pos_diff = p1.position().threevec() - p2.position().threevec();
  printd("Particle %d<->%d position difference: %g %g %g %g [fm]\n",
         p1.id(), p2.id(), pos_diff.x1(), pos_diff.x2(), pos_diff.x3());
  ThreeVector mom_diff = p1.momentum().threevec() - p2.momentum().threevec();
  printd("Particle %d<->%d momentum difference: %g %g %g %g [fm]\n",
         p1.id(), p2.id(), mom_diff.x1(), mom_diff.x2(), mom_diff.x3());
  /* Zero momentum leads to infite distance. */
  if (fabs(mom_diff.sqr()) < really_small)
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


bool ScatterAction::is_elastic() const {
  return outgoing_particles_.size() == 2 &&
         outgoing_particles_[0].pdgcode() == incoming_particles_[0].pdgcode() &&
         outgoing_particles_[1].pdgcode() == incoming_particles_[1].pdgcode();
}


ProcessBranch ScatterAction::elastic_cross_section(float elast_par) {
  return ProcessBranch(incoming_particles_[0].type().pdgcode(),
                       incoming_particles_[1].type().pdgcode(), elast_par);
}


ProcessBranchList ScatterAction::resonance_cross_sections() {
  ProcessBranchList resonance_process_list;
  ParticleType type_particle1 = incoming_particles_[0].type(),
               type_particle2 = incoming_particles_[1].type();

  /* Isospin symmetry factor, by default 1 */
  int symmetryfactor = 1;
  /* The isospin symmetry factor is 2 if both particles are in the same
   * isospin multiplet. */
  if (type_particle1.pdgcode().iso_multiplet()
      == type_particle2.pdgcode().iso_multiplet()) {
    symmetryfactor = 2;
  }

  const double s = mandelstam_s();
  const double p_cm_sqr = cm_momentum_squared();

  /* Find all the possible resonances */
  for (const ParticleType &type_resonance : ParticleType::list_all()) {
    /* Not a resonance, go to next type of particle */
    if (type_resonance.is_stable()) {
      continue;
    }

    /* Same resonance as in the beginning, ignore */
    if ((!type_particle1.is_stable()
         && type_resonance.pdgcode() == type_particle1.pdgcode())
        || (!type_particle2.is_stable()
            && type_resonance.pdgcode() == type_particle2.pdgcode())) {
      continue;
    }

    float resonance_xsection
      = symmetryfactor * two_to_one_formation(type_particle1, type_particle2,
                                              type_resonance, s, p_cm_sqr);

    /* If cross section is non-negligible, add resonance to the list */
    if (resonance_xsection > really_small) {
      resonance_process_list.push_back(ProcessBranch(type_resonance.pdgcode(),
                                                     resonance_xsection));

      printd("Found resonance %s (%s) with mass %f and width %f.\n",
             type_resonance.pdgcode().string().c_str(),
             type_resonance.name().c_str(),
             type_resonance.mass(), type_resonance.width_at_pole());
      printd("2->1 with original particles: %s %s Charges: %i %i \n",
             type_particle1.name().c_str(), type_particle2.name().c_str(),
             type_particle1.charge(), type_particle2.charge());
    }
  }
  return resonance_process_list;
}


void ScatterAction::momenta_exchange() {
  outgoing_particles_[0] = incoming_particles_[0];
  outgoing_particles_[1] = incoming_particles_[1];

  ParticleData *p1 = &outgoing_particles_[0];
  ParticleData *p2 = &outgoing_particles_[1];

  /* Boost to CM frame. */
  ThreeVector velocity_CM = beta_cm();
  p1->boost(velocity_CM);
  p2->boost(velocity_CM);

  /* debug output */
  printd_momenta("center of momenta 1", *p1);
  printd_momenta("center of momenta 2", *p2);

  /* We are in the center of momentum,
     hence this is equal for both particles. */
  const double momentum_radial = p1->momentum().abs3();

  /* Particle exchange momenta and scatter to random direction.
   * XXX: Angles should be sampled from differential cross section
   * of this process. */
  Angles phitheta;
  phitheta.distribute_isotropically();
  printd("Random momentum: %g %g %g %g \n", momentum_radial, phitheta.phi(),
        phitheta.costheta(), phitheta.sintheta());

  /* Only direction of 3-momentum, not magnitude, changes in CM frame.
   * Thus particle energies remain the same (Lorentz boost will change them for
   * computational frame, however). */
  p1->set_3momentum(phitheta.threevec() * momentum_radial);
  p2->set_3momentum(-phitheta.threevec() * momentum_radial);

  /* debug output */
  printd_momenta("exchanged momenta 1", *p1);
  printd_momenta("exchanged momenta 2", *p2);

  /* Boost back. */
  p1->boost(-velocity_CM);
  p2->boost(-velocity_CM);
}


void ScatterAction::resonance_formation() {

  switch (outgoing_particles_.size()) {
  case 1:
    /* 1 particle in final state: Center-of-momentum frame of initial
     * particles is the rest frame of the resonance.
     */
    outgoing_particles_[0].set_4momentum(FourVector(sqrt_s(), 0., 0., 0.));

    printd("Momentum of the new particle: %g %g %g %g \n",
           outgoing_particles_[0].momentum().x0(),
           outgoing_particles_[0].momentum().x1(),
           outgoing_particles_[0].momentum().x2(),
           outgoing_particles_[0].momentum().x3());
    break;
  case 2:
    /* 2 particles in final state: Sample the particle momenta. */
    sample_cms_momenta();
    break;
  default:
    std::string s = "resonance_formation: "
                    "Incorrect number of particles in final state: ";
    s += std::to_string(outgoing_particles_.size()) + " (";
    s += incoming_particles_[0].pdgcode().string() + " + ";
    s += incoming_particles_[1].pdgcode().string() + ")";
    throw InvalidResonanceFormation(s);
  }
}


ProcessBranch ScatterActionBaryonBaryon::elastic_cross_section(float elast_par) {

  const PdgCode &pdg1 = incoming_particles_[0].type().pdgcode();
  const PdgCode &pdg2 = incoming_particles_[1].type().pdgcode();

  const double s = mandelstam_s();

  if ((pdg1.iso_multiplet() == 0x1112) &&
      (pdg2.iso_multiplet() == 0x1112)) {
    /* Nucleon-Nucleon scattering: use parametrized cross sections. */
    float sig_el;
    if (pdg1 == pdg2) {                          /* pp */
      sig_el = pp_elastic(s);
    } else if (pdg1.is_antiparticle_of(pdg2)) {  /* ppbar */
      sig_el = ppbar_elastic(s);
    } else {                                     /* np */
      sig_el = np_elastic(s);
    }
    if (sig_el>0.) {
      return ProcessBranch(pdg1, pdg2, sig_el);
    } else {
      std::stringstream ss;
      ss << "problem in CrossSections::elastic: " << pdg1.string().c_str()
        << " " << pdg2.string().c_str() << " " << pdg1.spin() << " "
        << pdg2.spin() << " " << sig_el << " " << s;
      throw std::runtime_error(ss.str());
    }
  } else {
    /* Default: Fall back to parent routine. */
    return ScatterAction::elastic_cross_section(elast_par);
  }
}

ProcessBranchList ScatterActionBaryonBaryon::two_to_two_cross_sections() {
  ProcessBranchList process_list;
  ParticleType type_particle1 = incoming_particles_[0].type(),
               type_particle2 = incoming_particles_[1].type();

  if (type_particle1.pdgcode().iso_multiplet() == 0x1112 &&
      type_particle2.pdgcode().iso_multiplet() == 0x1112) {
    /* Nucleon+Nucleon: find all resonance production channels */
      process_list = NucNuc_to_NucRes (type_particle1, type_particle2);
  } else if (type_particle1.pdgcode().iso_multiplet() == 0x1112 ||
             type_particle2.pdgcode().iso_multiplet() == 0x1112) {
    /* Nucleon+Resonance: absorption */
    process_list = NucRes_to_NucNuc (type_particle1, type_particle2);
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


ProcessBranchList ScatterActionBaryonBaryon::NucNuc_to_NucRes (
                            const ParticleType &type_particle1,
                            const ParticleType &type_particle2) {

  ProcessBranchList process_list;
  const double s = mandelstam_s();

  /* Loop over all baryon resonances. */
  for (const ParticleType &type_resonance :
       ParticleType::list_baryon_resonances()) {
    /* Loop over second particle (nucleon). */
    for (const ParticleType &second_type : ParticleType::list_nucleons()) {

      /* Check for charge conservation. */
      if (type_resonance.charge() + second_type.charge()
          != type_particle1.charge() + type_particle2.charge()) {
        continue;
      }

      int I_z = type_resonance.isospin3() + second_type.isospin3();

      /* Compute total isospin range with given initial and final particles. */
      int I_max = std::min(type_resonance.isospin() + second_type.isospin(),
                          type_particle1.isospin() + type_particle2.isospin());
      int I_min = std::max(abs(type_resonance.isospin()
                              - second_type.isospin()),
                          abs(type_particle1.isospin()
                              - type_particle2.isospin()));
      I_min = std::max(I_min, abs(I_z));

      /* Loop over total isospin in allowed range.
      * Use decrement of 2, since isospin is multiplied by 2. */
      double isospin_factor = 0.;
      for (int I_tot = I_max; I_tot >= I_min; I_tot -= 2) {
        isospin_factor = isospin_factor +
            isospin_clebsch_gordan(type_particle1, type_particle2, I_tot, I_z)
          * isospin_clebsch_gordan(type_resonance, second_type, I_tot, I_z);
      }

      /* If Clebsch-Gordan coefficient is zero, don't bother with the rest. */
      if (std::abs(isospin_factor) < really_small) {
        continue;
      }

      /* Integration limits. */
      double lower_limit = type_resonance.minimum_mass();
      double upper_limit = std::sqrt(s) - second_type.mass();
      /* Check the available energy (requiring it to be a little above the
      * threshold, because the integration will not work if it's too close). */
      if (upper_limit - lower_limit < 1E-3) {
        continue;
      }

      /* Calculate resonance production cross section
      * using the Breit-Wigner distribution as probability amplitude.
      * Integrate over the allowed resonance mass range. */
      IntegrandParameters params = {&type_resonance, second_type.mass(), s};
      printd("Process: %s %s -> %s %s\n", type_particle1.name().c_str(),
      type_particle2.name().c_str(), second_type.name().c_str(),
      type_resonance.name().c_str());
      printd("Limits: %g %g \n", lower_limit, upper_limit);
      double resonance_integral, integral_error;
      quadrature_1d(&spectral_function_integrand, &params,
                    lower_limit, upper_limit,
                    &resonance_integral, &integral_error);
      printd("Integral value: %g Error: %g \n", resonance_integral,
        integral_error);

      /* Cross section for 2->2 process with one resonance in final state.
       * Based on Eq. (46) in PhD thesis of J. Weil
       * (https://gibuu.hepforge.org/trac/chrome/site/files/phd/weil.pdf) */
      float xsection
           = isospin_factor * isospin_factor
             * resonance_integral / (s * std::sqrt(cm_momentum_squared()))
             * nn_to_resonance_matrix_element(s, type_resonance, second_type);

      if (xsection > really_small) {
        process_list.push_back(ProcessBranch(type_resonance.pdgcode(),
                                            second_type.pdgcode(), xsection));
        printd("Found 2->2 creation process for resonance %s (%s).\n",
              type_resonance.pdgcode().string().c_str(),
              type_resonance.name().c_str());
        printd("2->2 with original particles: %s %s Charges: %i %i \n",
              type_particle1.name().c_str(), type_particle2.name().c_str(),
              type_particle1.charge(), type_particle2.charge());
      }
    }
  }
  return process_list;
}


ProcessBranchList ScatterActionBaryonBaryon::NucRes_to_NucNuc (
                            const ParticleType &type_particle1,
                            const ParticleType &type_particle2) {

  const ParticleType *type_resonance, *type_nucleon;
  ProcessBranchList process_list;

  if (type_particle1.pdgcode().iso_multiplet() == 0x1112) {
    type_nucleon = &type_particle1;
    type_resonance = &type_particle2;
  } else if (type_particle2.pdgcode().iso_multiplet() == 0x1112) {
    type_nucleon = &type_particle2;
    type_resonance = &type_particle1;
  } else {
    throw std::runtime_error("Error: no nucleon found in NucRes_to_NucNuc!");
  }

  const double s = mandelstam_s();
  /* CM momentum in final state */
  double p_cm_final = sqrt(s - 4.*type_nucleon->mass_sqr())/2.;

  /* Loop over all nucleon charge states. */
  for (const ParticleType &nuc1 : ParticleType::list_nucleons()) {
    for (const ParticleType &nuc2 : ParticleType::list_nucleons()) {
      /* Check for charge conservation. */
      if (type_resonance->charge() + type_nucleon->charge()
          != nuc1.charge() + nuc2.charge()) {
        continue;
      }

      int I_z = type_resonance->isospin3() + type_nucleon->isospin3();

      /* Compute total isospin range with given initial and final particles. */
      int I_max = std::min(type_resonance->isospin() + type_nucleon->isospin(),
                          nuc1.isospin() + nuc2.isospin());
      int I_min = std::max(abs(type_resonance->isospin()
                               - type_nucleon->isospin()),
                          abs(nuc1.isospin() - nuc2.isospin()));
      I_min = std::max(I_min, abs(I_z));

      /* Loop over total isospin in allowed range.
      * Use decrement of 2, since isospin is multiplied by 2. */
      double isospin_factor = 0.;
      for (int I_tot = I_max; I_tot >= I_min; I_tot -= 2) {
        isospin_factor = isospin_factor +
            isospin_clebsch_gordan(nuc1, nuc2, I_tot, I_z)
          * isospin_clebsch_gordan(*type_resonance, *type_nucleon, I_tot, I_z);
      }

      /* If Clebsch-Gordan coefficient is zero, don't bother with the rest */
      if (std::abs(isospin_factor) < really_small) {
        continue;
      }

      /* Cross section for 2->2 resonance absorption.
       * See eqs. (B.6), (B.9) and (181) in the GiBUU review paper. */
      float xsection
        = isospin_factor * isospin_factor
          * p_cm_final / ( s * std::sqrt(cm_momentum_squared()))
          * nn_to_resonance_matrix_element(s, *type_resonance, *type_nucleon);

      if (xsection > really_small) {
        process_list.push_back(ProcessBranch(nuc1.pdgcode(), nuc2.pdgcode(),
                                             xsection));
        printd("Found 2->2 absorption process for resonance %s (%s).\n",
              type_resonance->pdgcode().string().c_str(),
              type_resonance->name().c_str());
        printd("2->2 with original particles: %s %s Charges: %i %i \n",
              type_particle1.name().c_str(), type_particle2.name().c_str(),
              type_particle1.charge(), type_particle2.charge());
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
