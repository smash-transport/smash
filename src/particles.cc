/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include <cmath>
#include <cstdio>
#include <iostream>
#include <sstream>

#include "include/angles.h"
#include "include/constants.h"
#include "include/distributions.h"
#include "include/fourvector.h"
#include "include/macros.h"
#include "include/outputroutines.h"
#include "include/particledata.h"
#include "include/particles.h"
#include "include/pdgcode.h"
#include <assert.h>

namespace Smash {

/* boost_CM - boost to center of momentum */
void boost_CM(ParticleData *particle1, ParticleData *particle2,
  ThreeVector &velocity) {
  FourVector momentum1(particle1->momentum()), momentum2(particle2->momentum());
  FourVector position1(particle1->position()), position2(particle2->position());
  double cms_energy = momentum1.x0() + momentum2.x0();

  // determine CMS velocity
  velocity = (momentum1.threevec() + momentum2.threevec()) / cms_energy;

  // Boost the momenta into CMS frame
  momentum1 = momentum1.LorentzBoost(velocity);
  momentum2 = momentum2.LorentzBoost(velocity);

  // Boost the positions into CMS frame
  position1 = position1.LorentzBoost(velocity);
  position2 = position2.LorentzBoost(velocity);

  particle1->set_momentum(momentum1);
  particle1->set_position(position1);
  particle2->set_momentum(momentum2);
  particle2->set_position(position2);
}

/* boost_back_CM - boost back from center of momentum */
void boost_back_CM(ParticleData *particle1, ParticleData *particle2,
                   const ThreeVector &velocity_orig) {
  FourVector momentum1(particle1->momentum()), momentum2(particle2->momentum());
  FourVector position1(particle1->position()), position2(particle2->position());

  ThreeVector velocity = - velocity_orig;

  /* Boost the momenta back to lab frame */
  momentum1 = momentum1.LorentzBoost(velocity);
  momentum2 = momentum2.LorentzBoost(velocity);

  /* Boost the positions back to lab frame */
  position1 = position1.LorentzBoost(velocity);
  position2 = position2.LorentzBoost(velocity);

  particle1->set_momentum(momentum1);
  particle1->set_position(position1);
  particle2->set_momentum(momentum2);
  particle2->set_position(position2);
}

/* particle_distance - measure distance between two particles
 *                     in center of momentum
 */
double particle_distance(ParticleData *particle_orig1,
  ParticleData *particle_orig2) {
  /* Copy the particles in order to boost them and to forget the copy */
  ParticleData particle1 = *particle_orig1, particle2 = *particle_orig2;
  ThreeVector velocity_CM;

  /* boost particles in center of momenta frame */
  boost_CM(&particle1, &particle2, velocity_CM);
  FourVector position_difference = particle1.position() - particle2.position();
  printd("Particle %d<->%d position difference: %g %g %g %g [fm]\n",
    particle1.id(), particle2.id(), position_difference.x0(),
    position_difference.x1(), position_difference.x2(),
    position_difference.x3());

  FourVector momentum_difference = particle1.momentum() - particle2.momentum();
  printd("Particle %d<->%d momentum difference: %g %g %g %g [fm]\n",
    particle1.id(), particle2.id(), momentum_difference.x0(),
    momentum_difference.x1(), momentum_difference.x2(),
    momentum_difference.x3());
  /* zero momentum leads to infite distance */
  if (fabs(momentum_difference.x1()) < really_small
      && fabs(momentum_difference.x2()) < really_small
      && fabs(momentum_difference.x3()) < really_small)
    return  position_difference.sqr3();

  /* UrQMD squared distance criteria:
   * arXiv:nucl-th/9803035 (3.27): in center of momemtum frame
   * position of particle a: x_a
   * position of particle b: x_b
   * velocity of particle a: v_a
   * velocity of particle b: v_b
   * d^2_{coll} = (x_a - x_b)^2 - ((x_a - x_a) . (v_a - v_b))^2 / (v_a - v_b)^2
   */
  return position_difference.sqr3()
    - (position_difference.threevec()*momentum_difference.threevec())
      * (position_difference.threevec()*momentum_difference.threevec())
      / momentum_difference.sqr3();
}

/* time_collision - measure collision time of two particles */
double collision_time(const ParticleData &particle1,
  const ParticleData &particle2) {
  /* UrQMD collision time
   * arXiv:1203.4418 (5.15): in computational frame
   * position of particle a: x_a
   * position of particle b: x_b
   * momentum of particle a: p_a
   * momentum of particle b: p_b
   * t_{coll} = - (x_a - x_b) . (p_a - p_b) / (p_a - p_b)^2
   */
  FourVector position_difference = particle1.position()
    - particle2.position();
  printd("Particle %d<->%d position difference: %g %g %g %g [fm]\n",
    particle1.id(), particle2.id(), position_difference.x0(),
    position_difference.x1(), position_difference.x2(),
    position_difference.x3());
  FourVector velocity_difference = particle1.momentum()
    / particle1.momentum().x0()
    - particle2.momentum() / particle2.momentum().x0();
  printd("Particle %d<->%d velocity difference: %g %g %g %g [fm]\n",
    particle1.id(), particle2.id(), velocity_difference.x0(),
    velocity_difference.x1(), velocity_difference.x2(),
    velocity_difference.x3());
  /* zero momentum leads to infite distance, particles are not approaching */
  if (fabs(velocity_difference.x1()) < really_small
      && fabs(velocity_difference.x2()) < really_small
      && fabs(velocity_difference.x3()) < really_small)
    return -1.0;
  return - position_difference.threevec()*velocity_difference.threevec()
           / velocity_difference.sqr3();
}

/* momenta_exchange - soft scattering */
void momenta_exchange(ParticleData *particle1, ParticleData *particle2) {
  /* debug output */
  printd_momenta("center of momenta 1", *particle1);
  printd_momenta("center of momenta 2", *particle2);

  /* center of momentum hence this is equal for both particles */
  const double momentum_radial = sqrt(particle1->momentum().x1()
    * particle1->momentum().x1() + particle1->momentum().x2() *
    particle1->momentum().x2() + particle1->momentum().x3() *
    particle1->momentum().x3());

  /* particle exchange momenta and scatter to random direction */
  /* XXX: Angles should be sampled from differential cross section
   * of this process
   */
  Angles phitheta;
  phitheta.distribute_isotropically();
  printd("Random momentum: %g %g %g %g \n", momentum_radial, phitheta.phi(),
        phitheta.costheta(), phitheta.sintheta());

  /* Only direction of 3-momentum, not magnitude, changes in CM frame,
   * thus particle energies remain the same (Lorentz boost will change them for
   * computational frame, however)
   */
  const FourVector momentum1(particle1->momentum().x0(),
     momentum_radial * phitheta.x(),
     momentum_radial * phitheta.y(),
     momentum_radial * phitheta.z());
  particle1->set_momentum(momentum1);
  const FourVector momentum2(particle2->momentum().x0(),
    - momentum_radial * phitheta.x(),
    - momentum_radial * phitheta.y(),
    - momentum_radial * phitheta.z());
  particle2->set_momentum(momentum2);

  /* debug output */
  printd_momenta("exchanged momenta 1", *particle1);
  printd_momenta("exchanged momenta 2", *particle2);
}

void sample_cms_momenta(ParticleData *particle1, ParticleData *particle2,
                        const double cms_energy, const double mass1,
                        const double mass2) {
  double energy1 = (cms_energy * cms_energy + mass1 * mass1 - mass2 * mass2) /
                   (2.0 * cms_energy);
  double momentum_radial = sqrt(energy1 * energy1 - mass1 * mass1);
  if (!(momentum_radial > 0.0))
    printf("Warning: radial momenta %g \n", momentum_radial);
  /* XXX: Angles should be sampled from differential cross section
   * of this process
   */
  Angles phitheta;
  phitheta.distribute_isotropically();
  if (!(energy1 > mass1)) {
    printf("Particle %s radial momenta %g phi %g cos_theta %g\n",
           particle1->pdgcode().string().c_str(), momentum_radial,
           phitheta.phi(), phitheta.costheta());
    printf("Etot: %g m_a: %g m_b %g E_a: %g\n", cms_energy, mass1, mass2,
           energy1);
  }
  /* We use fourvector to set 4-momentum, as setting it
   * with doubles requires that particle uses its
   * pole mass, which is not generally true for resonances
   */
  FourVector momentum1(energy1, momentum_radial * phitheta.x(),
                       momentum_radial * phitheta.y(),
                       momentum_radial * phitheta.z());
  particle1->set_momentum(momentum1);

  FourVector momentum2(cms_energy - energy1, -momentum_radial * phitheta.x(),
                       -momentum_radial * phitheta.y(),
                       -momentum_radial * phitheta.z());
  particle2->set_momentum(momentum2);

  printd("p0: %g %g \n", momentum1.x0(), momentum2.x0());
  printd("p1: %g %g \n", momentum1.x1(), momentum2.x1());
  printd("p2: %g %g \n", momentum1.x2(), momentum2.x2());
  printd("p3: %g %g \n", momentum1.x3(), momentum2.x3());
}

Particles::Particles(const std::string &particles,
                     const std::string &decaymodes)
    : types_(load_particle_types(particles)),
      all_decay_modes_(load_decaymodes(decaymodes)) {}

namespace {/*{{{*/
std::string trim(const std::string &s) {
  const auto begin = s.find_first_not_of(" \t\n\r");
  if (begin == std::string::npos) {
    return {};
  }
  const auto end = s.find_last_not_of(" \t\n\r");
  return s.substr(begin, end - begin + 1);
}
struct Line {/*{{{*/
  Line() = default;
  Line(int n, std::string &&t) : number(n), text(std::move(t)) {
  }
  int number;
  std::string text;
};/*}}}*/

std::string build_error_string(std::string message, const Line &line) {/*{{{*/
  return message + " (on line " + std::to_string(line.number) + ": \"" +
         line.text + "\")";
}/*}}}*/

/**
 * Helper function for parsing particles.txt and decaymodes.txt.
 *
 * This function goes through an input stream line by line and removes
 * comments and empty lines. The remaining lines will be returned as a vector
 * of strings and linenumber pairs (Line).
 *
 * \param input an lvalue reference to an input stream
 */
std::vector<Line> line_parser(const std::string &input) {/*{{{*/
  std::istringstream input_stream(input);
  std::vector<Line> lines;
  lines.reserve(50);

  std::string line;
  int line_number = 0;
  while (std::getline(input_stream, line)) {
    ++line_number;
    const auto hash_pos = line.find('#');
    if (hash_pos != std::string::npos) {
      // found a comment, remove it from the line and look further
      line = line.substr(0, hash_pos);
    }
    if (line.find_first_not_of(" \t") == std::string::npos) {
      // only whitespace (or nothing) on this line. Next, please.
      continue;
    }
    lines.emplace_back(line_number, std::move(line));
    line = std::string();
  }
  return std::move(lines);
}/*}}}*/

void ensure_all_read(std::istream &input, const Line &line) {/*{{{*/
  std::string tmp;
  input >> tmp;
  if (!input.eof()) {
    throw Particles::LoadFailure(
        build_error_string("While loading the Particle data:\nGarbage (" + tmp +
                               ") at the remainder of the line.",
                           line));
  }
}/*}}}*/
}  // unnamed namespace/*}}}*/

Particles::ParticleTypeMap Particles::load_particle_types(  //{{{
    const std::string &input) {
  Particles::ParticleTypeMap types;
  for (const Line &line : line_parser(input)) {
    std::istringstream lineinput(line.text);
    std::string name;
    float mass, width;
    PdgCode pdgcode;
    lineinput >> name >> mass >> width >> pdgcode;
    if (lineinput.fail()) {
      throw LoadFailure(build_error_string(
          "While loading the Particle data:\nFailed to convert the input "
          "string to the expected data types.",
          line));
    }
    ensure_all_read(lineinput, line);

    printd("Setting particle type %s mass %g width %g pdgcode %s\n",
           name.c_str(), mass, width, pdgcode.string().c_str());
    printd("Setting particle type %s isospin %i/2 charge %i spin %i/2\n",
           name.c_str(), pdgcode.isospin_total(), pdgcode.charge(),
                                                  pdgcode.spin());

    types.insert(std::make_pair(
        pdgcode,
        ParticleType{name, mass, width, pdgcode}));
  }
  return std::move(types);
}/*}}}*/

Particles::DecayModesMap Particles::load_decaymodes(const std::string &input) {
  Particles::DecayModesMap decaymodes;
  PdgCode pdgcode = PdgCode::invalid();
  DecayModes decay_modes_to_add;
  float ratio_sum = 0.0;

  const auto end_of_decaymodes = [&]() {
    if (pdgcode == PdgCode::invalid()) {  // at the start of the file
      return;
    }
    if (decay_modes_to_add.empty()) {
      throw MissingDecays("No decay modes found for particle " +
                          pdgcode.string());
    }
    // XXX: why not just unconditionally call renormalize? (mkretz)
    /* Check if ratios add to 1 */
    if (fabs(ratio_sum - 1.0) > really_small) {
      /* They didn't; renormalize */
      printf("Particle %s:\n", pdgcode.string().c_str());
      decay_modes_to_add.renormalize(ratio_sum);
    }
    /* Add the list of decay modes for this particle type */
    decaymodes.insert(std::make_pair(pdgcode, decay_modes_to_add));
    /* Clean up the list for the next particle type */
    decay_modes_to_add.clear();
    ratio_sum = 0.0;
  };

  for (const Line &line : line_parser(input)) {
    const auto trimmed = trim(line.text);
    assert(!trimmed.empty());  // trim(line.text) is never empty - else
                               // line_parser is broken
    // if (trimmed.find_first_not_of("-0123456789") ==
    if (trimmed.find_first_of(" \t") ==
        std::string::npos) {  // a single record on one line signifies a new
                              // decay mode section
      end_of_decaymodes();
      pdgcode = PdgCode(trim(line.text));
      if (!is_particle_type_registered(pdgcode)) {
        throw ReferencedParticleNotFound(build_error_string(
            "Inconsistency: The particle with PDG id " +
                pdgcode.string() +
                " was not registered through particles.txt, but "
                "decaymodes.txt referenced it.",
            line));
      }
      assert(pdgcode != PdgCode::invalid());  // special value for start of file
    } else {
      std::istringstream lineinput(line.text);
      std::vector<PdgCode> decay_particles;
      decay_particles.reserve(4);
      float ratio;
      lineinput >> ratio;

      PdgCode pdg;
      lineinput >> pdg;
      while (lineinput) {
        if (!is_particle_type_registered(pdg)) {
          throw ReferencedParticleNotFound(build_error_string(
              "Inconsistency: The particle with PDG id " +
                  pdg.string() +
                  " was not registered through particles.txt, but "
                  "decaymodes.txt referenced it.",
              line));
        }
        decay_particles.push_back(pdg);
        lineinput >> pdg;
      }
      if (pdg != PdgCode::invalid()) {
        decay_particles.push_back(pdg);
      }
      if (lineinput.fail() && !lineinput.eof()) {
        throw LoadFailure(
            build_error_string("Parse error: expected a PdgCode ", line));
      }
      decay_particles.shrink_to_fit();
      decay_modes_to_add.add_mode(std::move(decay_particles), ratio);
      ratio_sum += ratio;
    }
  }
  end_of_decaymodes();
  return std::move(decaymodes);
}

void Particles::reset() {
  id_max_ = -1;
  data_.clear();
}

}  // namespace Smash
