/*
 *
 *    Copyright (c) 2014
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

#include "include/decaymodes.h"
#include "include/inputfunctions.h"
#include "include/iomanipulators.h"
#include "include/logging.h"
#include "include/pdgcode.h"
#include "include/stringfunctions.h"
#include "include/width.h"

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
  return {static_cast<uint16_t>(offset)};
}
#endif

std::vector<ParticleTypePtr> ParticleType::list_nucleons() {
  return {&find(0x2212), &find(0x2112)};
}

std::vector<ParticleTypePtr> ParticleType::list_baryon_resonances() {
  std::vector<ParticleTypePtr> list;
  list.reserve(4);  // currently we have only the Delta (with four charge states)

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
    throw std::runtime_error("PDG code "+pdgcode.string()+ " not found!");
  }
  //assert(found->pdgcode() == pdgcode);
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

#ifdef NDEBUG
ParticleType::ParticleType(std::string, float m, float w, PdgCode id)
    :
#else
ParticleType::ParticleType(std::string n, float m, float w, PdgCode id)
    : name_(fill_right(n, 3)),
#endif
      mass_(m),
      width_(w),
      pdgcode_(id),
      isospin_(pdgcode_.isospin_total()),
      charge_(pdgcode_.charge()) {}

/* Construct an antiparticle name-string from the given name-string for the
 * particle and its PDG code. */
static std::string antiname(const std::string &name, PdgCode code) {
  std::string basename, charge;

  if (name.find("++") != std::string::npos) {
    basename = name.substr(0, name.length()-2);
    charge = "--";
  } else if (name.find("⁺⁺") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁺⁺") + 1);
    charge = "⁻⁻";
  } else if (name.find("+") != std::string::npos) {
    basename = name.substr(0, name.length()-1);
    charge = "-";
  } else if (name.find("⁺") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁺") + 1);
    charge = "⁻";
  } else if (name.find("0") != std::string::npos) {
    basename = name.substr(0, name.length()-1);
    charge = "0";
  } else if (name.find("⁰") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁰") + 1);
    charge = "⁰";
  } else if (name.find("-") != std::string::npos) {
    basename = name.substr(0, name.length()-1);
    charge = "+";
  } else if (name.find("⁻") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁻") + 1);
    charge = "⁺";
  } else if (name.find("--") != std::string::npos) {
    basename = name.substr(0, name.length()-2);
    charge = "++";
  } else if (name.find("⁻⁻") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁻⁻") + 1);
    charge = "⁺⁺";
  } else {
    basename = name;
    charge = "";
  }

  constexpr char bar[] = "\u0305";
  if (code.baryon_number() != 0) {
    return basename+bar+charge;  // baryon
  } else if (code.charge() != 0) {
    return basename+charge;        // charged meson
  } else {
    return basename+bar+charge;  // neutral meson
  }
};

void ParticleType::create_type_list(const std::string &input) {  //{{{
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
      throw Particles::LoadFailure(build_error_string(
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

  std::sort(type_list.begin(), type_list.end(),
            [](const ParticleType &l,
               const ParticleType &r) { return l.pdgcode() < r.pdgcode(); });

  assert(nullptr == all_particle_types);
  all_particle_types = &type_list;  // note that type_list is a function-local
                                    // static and thus will live on until after
                                    // main().
}/*}}}*/


float ParticleType::minimum_mass() const {
  float minmass = mass();

  /* If the particle happens to be stable, just return the mass. */
  if (is_stable()) {
    return minmass;
  }
  /* Otherwise, find the lowest mass value needed in any decay mode */
  for (const auto &mode : decay_modes().decay_mode_list()) {
    minmass = std::min(minmass, mode.threshold());
  }
  return minmass;
}


float ParticleType::partial_width(const float m,
                                  const DecayBranch &mode) const {
  if (m < mode.threshold()) {
    return 0.;
  }
  float partial_width_at_pole = width_at_pole()*mode.weight();
  const ParticleType &t_a = *mode.particle_types()[0];
  const ParticleType &t_b = *mode.particle_types()[1];
  if (mode.particle_types().size() == 2) {
    /* two-body decays */
    if (t_a.is_stable() && t_b.is_stable()) {
      /* mass-dependent width for stable decay products */
      return width_Manley_stable(m, mass(), t_a.mass(), t_b.mass(),
                                 mode.angular_momentum(),
                                 partial_width_at_pole);
    }
    else if (t_a.is_stable()) {
      /* mass-dependent width for one unstable daughter */
      return width_Manley_semistable(m, mass(), t_a.mass(), t_b,
                                     mode.angular_momentum(),
                                     partial_width_at_pole);
    }
    else if (t_b.is_stable()) {
      /* mass-dependent width for one unstable daughter */
      return width_Manley_semistable(m, mass(), t_b.mass(), t_a,
                                     mode.angular_momentum(),
                                     partial_width_at_pole);
    }
    else {
      /* two unstable decay products: assume constant width */
      return partial_width_at_pole;
    }
  } else {
    /* three-body decays: asssume constant width */
    return partial_width_at_pole;
  }
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
  for (const auto &mode : decay_modes().decay_mode_list()) {
    w = w + partial_width(m, mode);
  }
  return w;
}

void ParticleType::check_consistency() {
  for (const ParticleType &ptype : ParticleType::list_all()) {
    if (!ptype.is_stable() && ptype.decay_modes().is_empty()) {
      throw std::runtime_error("Unstable particle " +
                                ptype.pdgcode().string() +
                               " has no decay chanels!");
    }
  }
}

ProcessBranchList ParticleType::get_partial_widths(const float m) const {
  float w = 0.;
  if (is_stable()) {
    return {};
  }
  /* Loop over decay modes and calculate all partial widths. */
  const auto &decay_mode_list = decay_modes().decay_mode_list();
  ProcessBranchList partial;
  partial.reserve(decay_mode_list.size());
  for (const auto &mode : decay_mode_list) {
    w = partial_width(m, mode);
    if (w > 0.) {
      partial.emplace_back(mode.particle_types(), w);
    }
  }
  return std::move(partial);
}

float ParticleType::get_partial_in_width(const float m,
                                         const ParticleData &p_a,
                                         const ParticleData &p_b) const {
  PdgCode pdg_a = p_a.type().pdgcode();
  PdgCode pdg_b = p_b.type().pdgcode();
  /* Get all decay modes. */
  const auto &decaymodes = decay_modes().decay_mode_list();

  /* Find the right one. */
  for (const auto &mode : decaymodes) {
    size_t decay_particles = mode.particle_types().size();
    if ( decay_particles > 3 ) {
      logger<LogArea::ParticleType>().warn("Not a 1->2 or 1->3 process!\n",
                                           "Number of decay particles: ",
                                           decay_particles);
    } else {
      if (decay_particles == 2 &&
          ((mode.particle_types()[0]->pdgcode() == pdg_a &&
            mode.particle_types()[1]->pdgcode() == pdg_b) ||
           (mode.particle_types()[0]->pdgcode() == pdg_b &&
            mode.particle_types()[1]->pdgcode() == pdg_a))) {
        /* Found: calculate width. */
        if (m < mode.threshold()) {
          return 0.;
        }
        float partial_width_at_pole = width_at_pole()*mode.weight();
        const ParticleType &t_a = ParticleType::find(pdg_a);
        const ParticleType &t_b = ParticleType::find(pdg_b);
        if (mode.particle_types().size() == 2) {
          /* two-body decays */
          if (t_a.is_stable() && t_b.is_stable()) {
            /* mass-dependent width for stable decay products */
            return width_Manley_stable(m, mass(), t_a.mass(), t_b.mass(),
                                      mode.angular_momentum(),
                                      partial_width_at_pole);
          }
          else if (t_a.is_stable()) {
            /* mass-dependent in-width for one unstable daughter */
            return in_width_Manley_semistable(m, mass(), t_a.mass(),
                                              p_b.effective_mass(), t_b,
                                              mode.angular_momentum(),
                                              partial_width_at_pole);
          }
          else if (t_b.is_stable()) {
            /* mass-dependent in-width for one unstable daughter */
            return in_width_Manley_semistable(m, mass(), t_b.mass(),
                                              p_a.effective_mass(), t_a,
                                              mode.angular_momentum(),
                                              partial_width_at_pole);
          }
          else {
            /* two unstable decay products: assume constant width */
            return partial_width_at_pole;
          }
        } else {
          /* three-body decays: asssume constant width */
          return partial_width_at_pole;
        }
      }
    }
  }
  /* Decay mode not found: width is zero. */
  return 0.;
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
