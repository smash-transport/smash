/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/particletype.h"

#include <assert.h>
#include <algorithm>
#include <map>
#include <vector>

#include "include/cxx14compat.h"
#include "include/decaymodes.h"
#include "include/inputfunctions.h"
#include "include/iomanipulators.h"
#include "include/isoparticletype.h"
#include "include/logging.h"
#include "include/particledata.h"
#include "include/pdgcode.h"
#include "include/processbranch.h"
#include "include/stringfunctions.h"

namespace Smash {

#ifdef SMASH_INLINE_LIST_ALL
const ParticleTypeList *all_particle_types = nullptr;
#else
namespace {
/// Global pointer to the Particle Type list.
const ParticleTypeList *all_particle_types = nullptr;
}  // unnamed namespace

const ParticleTypeList &ParticleType::list_all() {
  assert(all_particle_types);
  return *all_particle_types;
}

ParticleTypePtr ParticleType::operator&() const {
  // Calculate the offset via pointer subtraction:
  const auto offset = this - std::addressof(list_all()[0]);
  // Since we're using uint16_t for storing the index better be safe than sorry:
  // The offset must fit into the data type. If this ever fails we got a lot
  // more particle types than initially expected and you have to increase the
  // ParticleTypePtr storage to uint32_t.
  assert(offset >= 0 && offset < 0xffff);
  // After the assertion above the down-cast to uint16_t is safe:
  return ParticleTypePtr(static_cast<uint16_t>(offset));
}
#endif

std::vector<ParticleTypePtr> ParticleType::list_nucleons() {
  return {&find(0x2212), &find(0x2112)};
}

std::vector<ParticleTypePtr> ParticleType::list_baryon_resonances() {
  std::vector<ParticleTypePtr> list;
  list.reserve(10);
  for (const ParticleType &type_resonance : ParticleType::list_all()) {
    /* Only loop over baryon resonances. */
    if (type_resonance.is_stable()
        || type_resonance.pdgcode().baryon_number() != 1) {
      continue;
    }
    list.emplace_back(&type_resonance);
  }
  return list;
}

SMASH_CONST const ParticleType &ParticleType::find(PdgCode pdgcode) {
  const auto found = std::lower_bound(
      all_particle_types->begin(), all_particle_types->end(), pdgcode,
      [](const ParticleType &l, const PdgCode &r) { return l.pdgcode() < r; });
  if (found == all_particle_types->end() || found->pdgcode() != pdgcode) {
    throw PdgNotFoundFailure("PDG code " + pdgcode.string() + " not found!");
  }
  return *found;
}

SMASH_CONST bool ParticleType::exists(PdgCode pdgcode) {
  const auto found = std::lower_bound(
      all_particle_types->begin(), all_particle_types->end(), pdgcode,
      [](const ParticleType &l, const PdgCode &r) { return l.pdgcode() < r; });
  if (found != all_particle_types->end()) {
    return found->pdgcode() == pdgcode;
  }
  return false;
}

ParticleType::ParticleType(std::string n, float m, float w, PdgCode id)
    : name_(n),
      mass_(m),
      width_(w),
      pdgcode_(id),
      isospin_(pdgcode_.isospin_total()),
      charge_(pdgcode_.charge()) {}

/* Construct an antiparticle name-string from the given name-string for the
 * particle and its PDG code. */
static std::string antiname(const std::string &name, PdgCode code) {
  std::string basename, charge;

  if (name.find("⁺⁺") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁺⁺") + 1);
    charge = "⁻⁻";
  } else if (name.find("⁺") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁺") + 1);
    charge = "⁻";
  } else if (name.find("⁻⁻") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁻⁻") + 1);
    charge = "⁺⁺";
  } else if (name.find("⁻") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁻") + 1);
    charge = "⁺";
  } else if (name.find("⁰") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁰") + 1);
    charge = "⁰";
  } else {
    basename = name;
    charge = "";
  }

  // baryons & strange mesons: insert a bar
  constexpr char bar[] = "\u0305";
  if (code.baryon_number() != 0 || code.strangeness() != 0) {
    basename.insert(utf8::sequence_length(basename.begin()), bar);
  }

  return basename+charge;
}

void ParticleType::create_type_list(const std::string &input) {  // {{{
  const auto &log = logger<LogArea::ParticleType>();
  static ParticleTypeList type_list;
  type_list.clear();  // in case LoadFailure was thrown and caught and we should
                      // try again
  for (const Line &line : line_parser(input)) {
    std::istringstream lineinput(line.text);
    std::string name;
    float mass, width;
    PdgCode pdgcode;
    lineinput >> name >> mass >> width >> pdgcode;
    if (lineinput.fail()) {
      throw ParticleType::LoadFailure(build_error_string(
          "While loading the ParticleType data:\nFailed to convert the input "
          "string to the expected data types.",
          line));
    }
    ensure_all_read(lineinput, line);

    type_list.emplace_back(name, mass, width, pdgcode);
    log.debug() << "Setting     particle type: " << type_list.back();
    if (pdgcode.has_antiparticle()) {
      /* add corresponding antiparticle */
      PdgCode anti = pdgcode.get_antiparticle();
      name = antiname(name, pdgcode);
      type_list.emplace_back(name, mass, width, anti);
      log.debug() << "Setting antiparticle type: " << type_list.back();
    }
  }
  type_list.shrink_to_fit();

  /* Sort the type list by PDG code. */
  std::sort(type_list.begin(), type_list.end(),
            [](const ParticleType &l, const ParticleType &r) {
              return l.pdgcode() < r.pdgcode();
            });

  /* Look for duplicates. */
  PdgCode prev_pdg = 0;
  for (const auto& t : type_list) {
    if (t.pdgcode() == prev_pdg) {
      throw ParticleType::LoadFailure("Duplicate PdgCode in particles.txt: " +
                                      t.pdgcode().string());
    }
    prev_pdg = t.pdgcode();
  }

  if (all_particle_types != nullptr) {
    throw std::runtime_error("Error: Type list was already built!");
  }
  all_particle_types = &type_list;  // note that type_list is a function-local
                                    // static and thus will live on until after
                                    // main().

  // create all isospin multiplets
  for (const auto& t : type_list) {
    IsoParticleType::create_multiplet(t);
  }
}/*}}}*/


float ParticleType::minimum_mass() const {
  float minmass = mass();

  /* If the particle happens to be stable, just return the mass. */
  if (is_stable()) {
    return minmass;
  }
  /* Otherwise, find the lowest mass value needed in any decay mode */
  for (const auto &mode : decay_modes().decay_mode_list()) {
    minmass = std::min(minmass, mode->threshold());
  }
  return minmass;
}


float ParticleType::partial_width(const float m,
                                  const DecayBranch *mode) const {
  if (m < mode->threshold()) {
    return 0.;
  }
  float partial_width_at_pole = width_at_pole()*mode->weight();
  return mode->type().width(mass(), partial_width_at_pole, m);
}

const DecayModes &ParticleType::decay_modes() const {
  const auto offset = this - std::addressof(list_all()[0]);
  return (*DecayModes::all_decay_modes)[offset];
}

float ParticleType::total_width(const float m) const {
  float w = 0.;
  if (is_stable()) {
    return w;
  }
  /* Loop over decay modes and sum up all partial widths. */
  const auto &modes = decay_modes().decay_mode_list();
  for (unsigned int i = 0; i < modes.size(); i++) {
    w = w + partial_width(m, modes[i].get());
  }
  return w;
}

void ParticleType::check_consistency() {
  for (const ParticleType &ptype : ParticleType::list_all()) {
    if (!ptype.is_stable() && ptype.decay_modes().is_empty()) {
      throw std::runtime_error("Unstable particle " + ptype.name() +
                               " has no decay chanels!");
    }
  }
}

DecayBranchList ParticleType::get_partial_widths(const float m) const {
  const auto &decay_mode_list = decay_modes().decay_mode_list();
  if (decay_mode_list.size() == 0) {
    return {};
  }

  /* Loop over decay modes and calculate all partial widths. */
  DecayBranchList partial;
  partial.reserve(decay_mode_list.size());
  for (unsigned int i = 0; i < decay_mode_list.size(); i++) {
    const float w = partial_width(m, decay_mode_list[i].get());
    if (w > 0.) {
      partial.push_back(
          make_unique<DecayBranch>(decay_mode_list[i]->type(), w));
    }
  }
  return std::move(partial);
}

DecayBranchList ParticleType::get_partial_widths_hadronic(const float m) const {
  if (is_stable()) {
    return {};
  }
  /* Loop over decay modes and calculate all partial widths. */
  const auto &decay_mode_list = decay_modes().decay_mode_list();
  DecayBranchList partial;
  partial.reserve(decay_mode_list.size());
  for (unsigned int i = 0; i < decay_mode_list.size(); i++) {
    if ((!(is_dilepton(
             decay_mode_list[i]->type().particle_types()[0]->pdgcode(),
             decay_mode_list[i]->type().particle_types()[1]->pdgcode()))) &&
        (!(has_lepton_pair(
                decay_mode_list[i]->type().particle_types()[0]->pdgcode(),
                decay_mode_list[i]->type().particle_types()[1]->pdgcode(),
                decay_mode_list[i]->type().particle_types()[2]->pdgcode())))) {
      const float w = partial_width(m, decay_mode_list[i].get());
      if (w > 0.) {
        partial.push_back(
            make_unique<DecayBranch>(decay_mode_list[i]->type(), w));
      }
    }
  }
  return std::move(partial);
}

DecayBranchList ParticleType::get_partial_widths_dilepton(const float m) const {
  const auto &decay_mode_list = decay_modes().decay_mode_list();
  if (decay_mode_list.size() == 0) {
    return {};
  }
  /* Loop over decay modes and calculate all partial widths. */
  DecayBranchList partial;
  partial.reserve(decay_mode_list.size());
  for (unsigned int i = 0; i < decay_mode_list.size(); i++) {
    if ((is_dilepton(
                decay_mode_list[i]->type().particle_types()[0]->pdgcode(),
                decay_mode_list[i]->type().particle_types()[1]->pdgcode())) ||
        (has_lepton_pair(
                 decay_mode_list[i]->type().particle_types()[0]->pdgcode(),
                 decay_mode_list[i]->type().particle_types()[1]->pdgcode(),
                 decay_mode_list[i]->type().particle_types()[2]->pdgcode()))) {
      const float w = partial_width(m, decay_mode_list[i].get());
      if (w > 0.) {
        partial.push_back(
            make_unique<DecayBranch>(decay_mode_list[i]->type(), w));
      }
    }
  }
  return std::move(partial);
}

float ParticleType::get_partial_in_width(const float m,
                                         const ParticleData &p_a,
                                         const ParticleData &p_b) const {
  /* Get all decay modes. */
  const auto &decaymodes = decay_modes().decay_mode_list();

  /* Find the right one(s) and add up corresponding widths. */
  float w = 0.;
  for (const auto &mode : decaymodes) {
    float partial_width_at_pole = width_at_pole()*mode->weight();
    if (mode->type().has_particles(p_a.type(), p_b.type())) {
      w += mode->type().in_width(mass(), partial_width_at_pole, m,
                                 p_a.effective_mass(), p_b.effective_mass());
    }
  }
  return w;
}

std::ostream &operator<<(std::ostream &out, const ParticleType &type) {
  const PdgCode &pdg = type.pdgcode();
  return out << type.name() << std::setfill(' ') << std::right
             << "[mass:" << field<6> << type.mass()
             << ", width:" << field<6> << type.width_at_pole()
             << ", PDG:" << field<6> << pdg
             << ", Isospin:" << field<2> << pdg.isospin_total()
             << "/2, Charge:" << field<3> << pdg.charge()
             << ", Spin:" << field<2> << pdg.spin() << "/2]";
}

}  // namespace Smash
