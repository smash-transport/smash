/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/action.h"

#include "include/angles.h"
#include "include/constants.h"
#include "include/logging.h"
#include "include/parametrizations.h"
#include "include/pdgcode.h"
#include "include/random.h"
#include "include/resonances.h"

namespace Smash {

ScatterAction::ScatterAction(const ParticleData &in_part1,
                             const ParticleData &in_part2,
                             float time_of_execution)
    : Action({in_part1, in_part2}, time_of_execution) {}


void ScatterAction::perform(Particles *particles, size_t &id_process) {
  const auto &log = logger<LogArea::ScatterAction>();

  /* Relevant particle IDs for the collision. */
  int id1 = incoming_particles_[0].id();
  int id2 = incoming_particles_[1].id();

  log.debug("Process ", id_process, " particles:\n", incoming_particles_[0],
            incoming_particles_[1]);

  /* Decide for a particular final state. */
  outgoing_particles_ = choose_channel();

  if (is_elastic()) {
    /* 2->2 elastic scattering */
    log.debug("Process: Elastic collision.");

    momenta_exchange();

    // store the process id in the Particle data
    outgoing_particles_[0].set_id_process(id_process);
    outgoing_particles_[1].set_id_process(id_process);

    particles->data(id1) = outgoing_particles_[0];
    particles->data(id2) = outgoing_particles_[1];
  } else {
    /* resonance formation */
    log.debug("Process: Resonance formation.");

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

      log.debug("Resonance: ", new_particle);

      new_particle.set_id(particles->add_data(new_particle));
    }

    /* Remove the initial particles */
    particles->remove(id1);
    particles->remove(id2);

    log.debug("Particle map now has ", particles->size(), " elements.");
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
  const auto &log = logger<LogArea::ScatterAction>();
  // local copy of particles (since we need to boost them)
  ParticleData p1 = incoming_particles_[0];
  ParticleData p2 = incoming_particles_[1];
  /* Boost particles to center-of-momentum frame. */
  ThreeVector velocity = beta_cm();
  p1.boost(velocity);
  p2.boost(velocity);
  ThreeVector pos_diff = p1.position().threevec() - p2.position().threevec();
  ThreeVector mom_diff = p1.momentum().threevec() - p2.momentum().threevec();
  log.debug("Particle ", incoming_particles_, " position difference [fm]: ",
            pos_diff, ", momentum difference [GeV]: ", mom_diff);
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
  const auto &log = logger<LogArea::ScatterAction>();
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

      log.debug("Found resonance: ", type_resonance);
      log.debug("2->1 with original particles: ", type_particle1,
                type_particle2);
    }
  }
  return resonance_process_list;
}


void ScatterAction::momenta_exchange() {
  const auto &log = logger<LogArea::ScatterAction>();
  outgoing_particles_[0] = incoming_particles_[0];
  outgoing_particles_[1] = incoming_particles_[1];

  ParticleData *p1 = &outgoing_particles_[0];
  ParticleData *p2 = &outgoing_particles_[1];

  /* Boost to CM frame. */
  ThreeVector velocity_CM = beta_cm();
  p1->boost(velocity_CM);
  p2->boost(velocity_CM);

  /* debug output */
  log.debug("center of momenta 1", p1->momentum());
  log.debug("center of momenta 2", p2->momentum());

  /* We are in the center of momentum,
     hence this is equal for both particles. */
  const double momentum_radial = p1->momentum().abs3();

  /* Particle exchange momenta and scatter to random direction.
   * XXX: Angles should be sampled from differential cross section
   * of this process. */
  Angles phitheta;
  phitheta.distribute_isotropically();
  log.debug("Random momentum: ", momentum_radial, " ", phitheta);

  /* Only direction of 3-momentum, not magnitude, changes in CM frame.
   * Thus particle energies remain the same (Lorentz boost will change them for
   * computational frame, however). */
  p1->set_3momentum(phitheta.threevec() * momentum_radial);
  p2->set_3momentum(-phitheta.threevec() * momentum_radial);

  /* debug output */
  log.debug("exchanged momenta 1", p1->momentum());
  log.debug("exchanged momenta 2", p2->momentum());

  /* Boost back. */
  p1->boost(-velocity_CM);
  p2->boost(-velocity_CM);
}


void ScatterAction::resonance_formation() {
  const auto &log = logger<LogArea::ScatterAction>();

  switch (outgoing_particles_.size()) {
  case 1:
    /* 1 particle in final state: Center-of-momentum frame of initial
     * particles is the rest frame of the resonance.
     */
    outgoing_particles_[0].set_4momentum(FourVector(sqrt_s(), 0., 0., 0.));

    log.debug("Momentum of the new particle: ",
              outgoing_particles_[0].momentum());
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
  ProcessBranchList resonance_process_list;
  ParticleType type_particle1 = incoming_particles_[0].type(),
               type_particle2 = incoming_particles_[1].type();

  const double s = mandelstam_s();
  const double p_cm_sqr = cm_momentum_squared();

  /* Find all the possible resonances */
  for (const ParticleType &type_resonance : ParticleType::list_all()) {
    /* Not a resonance, go to next type of particle */
    if (type_resonance.is_stable()) {
      continue;
    }

    size_t two_to_two_processes
        = two_to_two_formation(type_particle1, type_particle2, type_resonance,
                               s, p_cm_sqr, &resonance_process_list);
    if (two_to_two_processes > 0) {
      const auto &log = logger<LogArea::ScatterAction>();
      log.debug("Found ", two_to_two_processes,
                " 2->2 processes for resonance ", type_resonance);
      log.debug("2->2 with original particles: ", type_particle1,
                type_particle2);
    }
  }

  return resonance_process_list;
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
