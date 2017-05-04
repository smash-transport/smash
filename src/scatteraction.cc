/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteraction.h"

#include "Pythia8/Pythia.h"

#include "include/constants.h"
#include "include/cxx14compat.h"
#include "include/fpenvironment.h"
#include "include/kinematics.h"
#include "include/logging.h"
#include "include/pdgcode.h"
#include "include/random.h"

namespace Smash {

ScatterAction::ScatterAction(const ParticleData &in_part_a,
                             const ParticleData &in_part_b,
                             double time, bool isotropic,
                             float formation_time)
    : Action({in_part_a, in_part_b}, time),
      total_cross_section_(0.), isotropic_(isotropic),
      formation_time_(formation_time) {}

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
      elastic_scattering();
      break;
    case ProcessType::TwoToOne:
      /* resonance formation */
      /* processes computed in the center of momenta */
      resonance_formation();
      break;
    case ProcessType::TwoToTwo:
      /* 2->2 inelastic scattering */
      /* Sample the particle momenta in CM system. */
      sample_2body_phasespace();
      break;
    case ProcessType::String:
      /* string excitation */
      string_excitation();
      break;
    default:
      throw InvalidScatterAction(
        "ScatterAction::generate_final_state: Invalid process type "
        + std::to_string(static_cast<int>(process_type_)) + " was requested. "
        + "(PDGcode1=" + incoming_particles_[0].pdgcode().string()
        + ", PDGcode2=" + incoming_particles_[1].pdgcode().string()
        + ")");
  }

  /* Set positions & boost to computational frame. */
  for (ParticleData &new_particle : outgoing_particles_) {
    if (proc->get_type() != ProcessType::Elastic) {
      new_particle.set_4position(middle_point);
    }
    new_particle.boost_momentum(-beta_cm());
  }
}


void ScatterAction::add_all_processes(float elastic_parameter,
                                      bool two_to_one, bool two_to_two,
                                      double low_snn_cut,
                                      bool strings_switch) {
  if (two_to_one) {
    /* resonance formation (2->1) */
    add_collisions(resonance_cross_sections());
  }
  if (two_to_two) {
    /** Elastic collisions between two nucleons with sqrt_s() below
     * low_snn_cut can not happen*/
    if (!incoming_particles_[0].type().is_nucleon() ||
        !incoming_particles_[1].type().is_nucleon() ||
        !(incoming_particles_[0].type().antiparticle_sign() ==
          incoming_particles_[1].type().antiparticle_sign()) ||
        sqrt_s() >= low_snn_cut) {
        add_collision(elastic_cross_section(elastic_parameter));
    }
    /* 2->2 (inelastic) */
    add_collisions(two_to_two_cross_sections());
  }
  /* string excitation: the sqrt(s) cut-off is the sum of the masses of the
   * incoming particles + 2 GeV, which is given by PYTHIA as the
   * minimum energy that needs to be available for particle production */
  if (strings_switch &&
     sqrt_s() >= incoming_particles_[0].type().mass() +
                 incoming_particles_[1].type().mass() + 2.) {
  /* Only allow string excitation for the particles that PYTHIA can cope
   * with, i.e. p/n, p/nbar, pi+, pi- and pi0 */
    bool a_in_pythia = false;
    bool b_in_pythia = false;
    if (incoming_particles_[0].type().is_nucleon() ||
        incoming_particles_[0].type().pdgcode().is_pion() ) {
        a_in_pythia = true;
    }
    if (incoming_particles_[1].type().is_nucleon() ||
        incoming_particles_[1].type().pdgcode().is_pion() ) {
        b_in_pythia = true;
    }
    if (a_in_pythia && b_in_pythia) {
      add_collision(string_excitation_cross_section());
    }
  }
}


float ScatterAction::raw_weight_value() const {
  return total_cross_section_;
}


ThreeVector ScatterAction::beta_cm() const {
  return (incoming_particles_[0].momentum() +
          incoming_particles_[1].momentum()).velocity();
}

double ScatterAction::gamma_cm() const {
  return (1./std::sqrt(1-beta_cm().sqr()));
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

double ScatterAction::transverse_distance_sqr() const {
  // local copy of particles (since we need to boost them)
  ParticleData p_a = incoming_particles_[0];
  ParticleData p_b = incoming_particles_[1];
  /* Boost particles to center-of-momentum frame. */
  const ThreeVector velocity = beta_cm();
  p_a.boost(velocity);
  p_b.boost(velocity);
  const ThreeVector pos_diff = p_a.position().threevec() -
                               p_b.position().threevec();
  const ThreeVector mom_diff = p_a.momentum().threevec() -
                               p_b.momentum().threevec();

  const auto &log = logger<LogArea::ScatterAction>();
  log.debug("Particle ", incoming_particles_, " position difference [fm]: ",
            pos_diff, ", momentum difference [GeV]: ", mom_diff);

  const double dp2 = mom_diff.sqr();
  const double dr2 = pos_diff.sqr();
  /* Zero momentum leads to infite distance. */
  if (dp2 < really_small) {
    return dr2;
  }
  const double dpdr = pos_diff * mom_diff;

  /** UrQMD squared distance criterion:
   * \iref{Bass:1998ca} (3.27): in center of momentum frame
   * position of particle a: x_a
   * position of particle b: x_b
   * momentum of particle a: p_a
   * momentum of particle b: p_b
   * d^2_{coll} = (x_a - x_b)^2 - ((x_a - x_b) . (p_a - p_b))^2 / (p_a - p_b)^2
   */
  return dr2 - dpdr*dpdr/dp2;
}

CollisionBranchPtr ScatterAction::elastic_cross_section(float elast_par) {
  float elastic_xs;
  if (elast_par >= 0.) {
    // use constant elastic cross section from config file
    elastic_xs = elast_par;
  } else {
    // use parametrization
    elastic_xs = elastic_parametrization();
  }
  return make_unique<CollisionBranch>(incoming_particles_[0].type(),
                                      incoming_particles_[1].type(),
                                      elastic_xs, ProcessType::Elastic);
}

CollisionBranchPtr ScatterAction::string_excitation_cross_section() {
  const auto &log = logger<LogArea::ScatterAction>();
  /* Calculate string-excitation cross section:
   * Parametrized total minus all other present channels. */
  float sig_string = std::max(0.f, total_cross_section() - cross_section());
  log.debug("String cross section is: ", sig_string);
  return make_unique<CollisionBranch>(sig_string, ProcessType::String);
}


double ScatterAction::two_to_one_formation(const ParticleType &type_resonance,
                                          double srts, double cm_momentum_sqr) {
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
  const float partial_width = type_resonance.get_partial_in_width(srts,
                                incoming_particles_[0], incoming_particles_[1]);
  if (partial_width <= 0.f) {
    return 0.;
  }

  /* Calculate spin factor */
  const double spinfactor = static_cast<double>(type_resonance.spin() + 1)
    / ((type_particle_a.spin() + 1) * (type_particle_b.spin() + 1));
  const int sym_factor = (type_particle_a.pdgcode() ==
                          type_particle_b.pdgcode()) ? 2 : 1;
  /** Calculate resonance production cross section
   * using the Breit-Wigner distribution as probability amplitude.
   * See Eq. (176) in \iref{Buss:2011mx}. */
  return spinfactor * sym_factor * 2. * M_PI * M_PI / cm_momentum_sqr
         * type_resonance.spectral_function(srts)
         * partial_width * hbarc * hbarc / fm2_mb;
}


CollisionBranchList ScatterAction::resonance_cross_sections() {
  const auto &log = logger<LogArea::ScatterAction>();
  CollisionBranchList resonance_process_list;
  const ParticleType &type_particle_a = incoming_particles_[0].type();
  const ParticleType &type_particle_b = incoming_particles_[1].type();

  const double srts = sqrt_s();
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
                                                    srts, p_cm_sqr);

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


void ScatterAction::elastic_scattering() {
  // copy initial particles into final state
  outgoing_particles_[0] = incoming_particles_[0];
  outgoing_particles_[1] = incoming_particles_[1];
  // resample momenta
  sample_angles({outgoing_particles_[0].effective_mass(),
                 outgoing_particles_[1].effective_mass()});
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

  /* Set the formation time of the resonance to the larger formation time of the
   * incoming particles */
  const float t0 = incoming_particles_[0].formation_time();
  const float t1 = incoming_particles_[1].formation_time();
  const size_t index_tmax = (t0 > t1) ? 0 : 1;
  const float sc =
    incoming_particles_[index_tmax].cross_section_scaling_factor();
  if (t0 > time_of_execution_ || t1 > time_of_execution_) {
    outgoing_particles_[0].set_formation_time(std::max(t0, t1));
    outgoing_particles_[0].set_cross_section_scaling_factor(sc);
  }
  log.debug("Momentum of the new particle: ",
            outgoing_particles_[0].momentum());
}

/* This function will generate outgoing particles in CM frame
 * from a hard process. */
void ScatterAction::string_excitation() {
  //Check that there are 2 particles in incoming_particles_
  assert(incoming_particles_.size() == 2);
  const auto &log = logger<LogArea::Pythia>();
  // Disable floating point exception trap for Pythia
  {
  DisableFloatTraps guard;
  /* set all necessary parameters for Pythia
   * Create Pythia object */
  log.debug("Creating Pythia object.");
  static thread_local Pythia8::Pythia pythia(PYTHIA_XML_DIR, false);
  /* select only inelastic events: */
  pythia.readString("SoftQCD:inelastic = on");
  /* suppress unnecessary output */
  pythia.readString("Print:quiet = on");
  /* Create output of the Pythia particle list */
  // pythia.readString("Init:showAllParticleData = on");
  /* No resonance decays, since the resonances will be handled by SMASH */
  pythia.readString("HadronLevel:Decay = off");
  /* Set the random seed of the Pythia Random Number Generator.
   * Please note: Here we use a random number generated by the
   * SMASH, since every call of pythia.init should produce
   * different events. */
  pythia.readString("Random:setSeed = on");
  std::stringstream buffer1;
  buffer1 << "Random:seed = " << Random::canonical();
  pythia.readString(buffer1.str());
  /* set the incoming particles */
  std::stringstream buffer2;
  buffer2 << "Beams:idA = " << incoming_particles_[0].type().pdgcode();
  pythia.readString(buffer2.str());
  log.debug("First particle in string excitation: ",
            incoming_particles_[0].type().pdgcode());
  std::stringstream buffer3;
  buffer3 << "Beams:idB = " << incoming_particles_[1].type().pdgcode();
  log.debug("Second particle in string excitation: ",
            incoming_particles_[1].type().pdgcode());
  pythia.readString(buffer3.str());
  /* Calculate the center-of-mass energy of this collision */
  double sqrts = (incoming_particles_[0].momentum() +
                  incoming_particles_[1].momentum()).abs();
  std::stringstream buffer4;
  buffer4 << "Beams:eCM = " << sqrts;
  pythia.readString(buffer4.str());
  log.debug("Pythia call with eCM = ", buffer4.str());
  /* Initialize. */
  pythia.init();
  /* Short notation for Pythia event */
  Pythia8::Event &event = pythia.event;
  pythia.next();
  ParticleList new_intermediate_particles;
  for (int i = 0; i < event.size(); i++) {
    if (event[i].isFinal()) {
      if (event[i].isHadron()) {
        int pythia_id = event[i].id();
        log.debug("PDG ID from Pythia:", pythia_id);
        /* K_short and K_long need to be converted to K0
         * since SMASH only knows K0 */
        if (pythia_id == 310 || pythia_id == 130) {
          const float prob = Random::uniform(0.f, 1.f);
          if (prob <= 0.5f) {
            pythia_id = 311;
          } else {
            pythia_id = -311;
          }
        }
        const std::string s = std::to_string(pythia_id);
        PdgCode pythia_code(s);
        ParticleData new_particle(ParticleType::find(pythia_code));
        FourVector momentum;
        momentum.set_x0(event[i].e());
        momentum.set_x1(event[i].px());
        momentum.set_x2(event[i].py());
        momentum.set_x3(event[i].pz());
        new_particle.set_4momentum(momentum);
        log.debug("4-momentum from Pythia: ", momentum);
        new_intermediate_particles.push_back(new_particle);
        }
      }
    }
    /*
     * sort new_intermediate_particles according to z-Momentum
     */
    std::sort(new_intermediate_particles.begin(),
              new_intermediate_particles.end(),
              [&](ParticleData i, ParticleData j) {
              return fabs(i.momentum().x3()) > fabs(j.momentum().x3());
              });
    for (ParticleData data : new_intermediate_particles) {
      log.debug("Particle momenta after sorting: ", data.momentum());
      /* The hadrons are not immediately formed, currently a formation time of
       * 1 fm is universally applied and cross section is reduced to zero and
       * to a fraction corresponding to the valence quark content. Hadrons
       * containing a valence quark are determined by highest z-momentum. */
      log.debug("The formation time is: ", formation_time_, "fm/c.");
      /* Additional suppression factor to mimic coherence taken as 0.7
       * from UrQMD (CTParam(59) */
      const float suppression_factor = 0.7;
      if (incoming_particles_[0].is_baryon() ||
          incoming_particles_[1].is_baryon()) {
        if (data == 0) {
          data.set_cross_section_scaling_factor(suppression_factor * 0.66);
        } else if (data == 1) {
          data.set_cross_section_scaling_factor(suppression_factor * 0.34);
        } else {
          data.set_cross_section_scaling_factor(suppression_factor * 0.0);
        }
      } else {
        if (data == 0 || data == 1) {
          data.set_cross_section_scaling_factor(suppression_factor * 0.50);
        } else {
          data.set_cross_section_scaling_factor(suppression_factor * 0.0);
        }
      }
      data.set_formation_time(formation_time_ * gamma_cm());
      outgoing_particles_.push_back(data);
    }
    /* If the incoming particles already were unformed, the formation
     * times and cross section scaling factors need to be adjusted */
    const float tform_in = std::max(incoming_particles_[0].formation_time(),
                                    incoming_particles_[1].formation_time());
    if (tform_in > time_of_execution_) {
      const float fin = (incoming_particles_[0].formation_time() >
                         incoming_particles_[1].formation_time()) ?
                         incoming_particles_[0].cross_section_scaling_factor() :
                         incoming_particles_[1].cross_section_scaling_factor();
      for (size_t i = 0; i < outgoing_particles_.size(); i++) {
        const float tform_out = outgoing_particles_[i].formation_time();
        const float fout = outgoing_particles_[i].cross_section_scaling_factor();
        outgoing_particles_[i].set_cross_section_scaling_factor(fin * fout);
        /* If the unformed incoming particles' formation time is larger than
         * the current outgoing particle's formation time, then the latter
         * is overwritten by the former*/
        if (tform_in > tform_out) {
          outgoing_particles_[i].set_formation_time(tform_in);
        }
      }
    }
    /* Check momentum difference for debugging */
    FourVector in_mom = incoming_particles_[0].momentum() +
      incoming_particles_[1].momentum();
    FourVector out_mom;
    for (ParticleData data : outgoing_particles_) {
      out_mom = out_mom + data.momentum();
    }
    log.debug("Incoming momenta string:", in_mom);
    log.debug("Outgoing momenta string:", out_mom);
  }
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
