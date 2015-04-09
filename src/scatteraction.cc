/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteraction.h"

#include "include/angles.h"
#include "include/constants.h"
#include "include/cxx14compat.h"
#include "include/distributions.h"
#include "include/kinematics.h"
#include "include/logging.h"
#include "include/pdgcode.h"
#include "include/random.h"

namespace Smash {

ScatterAction::ScatterAction(const ParticleData &in_part_a,
                             const ParticleData &in_part_b,
                             float time_of_execution)
    : Action({in_part_a, in_part_b}, time_of_execution) {}


void ScatterAction::generate_final_state() {
  const auto &log = logger<LogArea::ScatterAction>();

  log.debug("Incoming particles: ", incoming_particles_);

  /* Decide for a particular final state. */
  const ProcessBranch* proc = choose_channel();
  process_type_ = proc->get_type();
  outgoing_particles_ = proc->particle_list();

  log.debug("Chosen channel: ", process_type_, outgoing_particles_);

  /* The production point of the new particles.  */
  FourVector middle_point = get_interaction_point();

  switch (process_type_) {
    case ProcessBranch::Elastic:
      /* 2->2 elastic scattering */
      log.debug("Process: Elastic collision.", process_type_);
      elastic_scattering();
      break;
    case ProcessBranch::TwoToOne:
      /* resonance formation */
      log.debug("Process: Resonance formation.", process_type_);
      /* processes computed in the center of momenta */
      resonance_formation();
      break;
    case ProcessBranch::TwoToTwo:
      /* 2->2 inelastic scattering */
      log.debug("Process: Inelastic scattering.", process_type_);
      /* Sample the particle momenta in CM system. */
      sample_cms_momenta();
      break;
    case ProcessBranch::String:
      /* string excitation */
      log.debug("Process: String Excitation.");
      /// string_excitation(incoming_particles_, outgoing_particles_);
      break;
    case ProcessBranch::None:
      log.debug("ProcessType None should not have been selected");
      break;
    case ProcessBranch::Decay:
      log.debug("ProcessType Decay should have been handled as DecayAction");
      break;
    default:
      throw InvalidScatterAction(
        "ScatterAction::perform: Unknown Process Type. "
        "ProcessType " + std::to_string(process_type_) +
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


CollisionBranch* ScatterAction::elastic_cross_section(float elast_par) {
  return new CollisionBranch(incoming_particles_[0].type(),
                             incoming_particles_[1].type(),
                             elast_par, ProcessBranch::Elastic);
}

CollisionBranch* ScatterAction::string_excitation_cross_section() {
  /* Calculate string-excitation cross section:
   * Parametrized total minus all other present channels. */
  /* TODO(weil): This is currently set to zero,
   * since Pythia is not yet implemented. */
  float sig_string = 0.f;
  // = std::max(0.f, total_cross_section() - total_weight_);

  return new CollisionBranch(sig_string, ProcessBranch::String);
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


ProcessBranchList ScatterAction::resonance_cross_sections() {
  const auto &log = logger<LogArea::ScatterAction>();
  ProcessBranchList resonance_process_list;
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
                                        ProcessBranch::TwoToOne));
      log.debug("Found resonance: ", type_resonance);
      log.debug("2->1 with original particles: ", type_particle_a,
                type_particle_b);
    }
  }
  return resonance_process_list;
}


void ScatterAction::elastic_scattering() {
  const auto &log = logger<LogArea::ScatterAction>();
  outgoing_particles_[0] = incoming_particles_[0];
  outgoing_particles_[1] = incoming_particles_[1];

  /* Determine absolute momentum in center-of-mass frame. */
  const double momentum_radial = cm_momentum();

  /* Particles exchange momenta and scatter to random direction
   * (isotropically). */
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


void ScatterAction::format_debug_output(std::ostream &out) const {
  out << "Scatter of " << incoming_particles_;
  if (outgoing_particles_.empty()) {
    out << " (not performed)";
  } else {
    out << " to " << outgoing_particles_;
  }
}


}  // namespace Smash
