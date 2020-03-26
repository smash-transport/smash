/*
 *    Copyright (c) 2015-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "smash/isoparticletype.h"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "smash/filelock.h"
#include "smash/integrate.h"
#include "smash/kinematics.h"
#include "smash/logging.h"

namespace smash {
static constexpr int LParticleType = LogArea::ParticleType::id;

static IsoParticleTypeList iso_type_list;
static std::vector<const IsoParticleType *> iso_baryon_resonances;

const std::vector<const IsoParticleType *>
IsoParticleType::list_baryon_resonances() {
  if (iso_baryon_resonances.empty()) {
    // Initialize.
    for (const auto &res : IsoParticleType::list_all()) {
      const auto baryon_number = res.states_[0]->pdgcode().baryon_number();
      if (res.states_[0]->is_stable() || (baryon_number <= 0)) {
        continue;
      }
      iso_baryon_resonances.push_back(&res);
    }
  }
  return iso_baryon_resonances;
}

IsoParticleType::IsoParticleType(const std::string &n, double m, double w,
                                 unsigned int s, Parity p)
    : name_(n), mass_(m), width_(w), spin_(s), parity_(p) {}

const IsoParticleTypeList &IsoParticleType::list_all() { return iso_type_list; }

/// Helper function for IsoParticleType::try_find and friends.
static IsoParticleType *try_find_private(const std::string &name) {
  auto found =
      std::lower_bound(iso_type_list.begin(), iso_type_list.end(), name,
                       [](const IsoParticleType &l, const std::string &r) {
                         return l.name() < r;
                       });
  if (found == iso_type_list.end() || found->name() != name) {
    return {};  // The default constructor creates an invalid pointer.
  }
  return &*found;
}

const IsoParticleType *IsoParticleType::try_find(const std::string &name) {
  return try_find_private(name);
}

const IsoParticleType &IsoParticleType::find(const std::string &name) {
  const auto found = try_find_private(name);
  if (!found) {
    throw ParticleNotFoundFailure("Isospin multiplet " + name + " not found!");
  }
  return *found;
}

IsoParticleType &IsoParticleType::find_private(const std::string &name) {
  auto found = try_find_private(name);
  if (!found) {
    throw ParticleNotFoundFailure("Isospin multiplet " + name +
                                  " not found (privately)!");
  }
  return *found;
}

bool IsoParticleType::exists(const std::string &name) {
  const auto found = try_find_private(name);
  return found;
}

/**
 * Construct the name-string for an isospin multiplet from the given
 * name-string for the particle.
 *
 * \param[in] name name-string of the particle
 * \return the name-string for an isospin multiplet
 */
static std::string multiplet_name(std::string name) {
  if (name.find("⁺⁺") != std::string::npos) {
    return name.substr(0, name.length() - sizeof("⁺⁺") + 1);
  } else if (name.find("⁺") != std::string::npos) {
    return name.substr(0, name.length() - sizeof("⁺") + 1);
  } else if (name.find("⁻⁻") != std::string::npos) {
    return name.substr(0, name.length() - sizeof("⁻⁻") + 1);
  } else if (name.find("⁻") != std::string::npos) {
    return name.substr(0, name.length() - sizeof("⁻") + 1);
  } else if (name.find("⁰") != std::string::npos) {
    return name.substr(0, name.length() - sizeof("⁰") + 1);
  } else {
    return name;
  }
}

const IsoParticleType *IsoParticleType::anti_multiplet() const {
  if (states_[0]->has_antiparticle()) {
    ParticleTypePtr anti = states_[0]->get_antiparticle();
    if (states_[0]->name() != multiplet_name(anti->name())) {
      return anti->iso_multiplet();
    } else {
      return nullptr;
    }
  } else {
    return nullptr;
  }
}

bool IsoParticleType::has_anti_multiplet() const {
  return anti_multiplet() != nullptr;
}

const ParticleTypePtr IsoParticleType::find_state(const std::string &n) {
  const IsoParticleType &multiplet = IsoParticleType::find(multiplet_name(n));
  auto found = std::find_if(multiplet.states_.begin(), multiplet.states_.end(),
                            [&n](ParticleTypePtr p) { return p->name() == n; });
  if (found == multiplet.states_.end()) {
    throw std::runtime_error("Isospin state " + n + " not found!");
  }
  return *found;
}

IsoParticleType *IsoParticleType::find(const ParticleType &type) {
  std::string multiname = multiplet_name(type.name());
  IsoParticleType &multiplet = find_private(multiname);
  return &multiplet;
}

void IsoParticleType::add_state(const ParticleType &type) {
  states_.push_back(&type);

  // check if isospin symmetry is fulfilled
  if (std::abs(mass() - type.mass()) > really_small) {
    logg[LParticleType].warn()
        << "Isospin symmetry is broken by mass of " << type.name() << ": "
        << type.mass() << " vs. " << mass();
  }
  if (std::abs(width() - type.width_at_pole()) > really_small) {
    logg[LParticleType].warn()
        << "Isospin symmetry is broken by width of " << type.name() << ": "
        << type.width_at_pole() << " vs. " << width();
  }
  if (spin() != type.spin()) {
    logg[LParticleType].error()
        << "Isospin symmetry is broken by spin of " << type.name() << ": "
        << type.spin() << " vs. " << spin();
  }
}

void IsoParticleType::create_multiplet(const ParticleType &type) {
  // create multiplet if it does not exist yet
  std::string multiname = multiplet_name(type.name());
  if (!exists(multiname)) {
    iso_type_list.emplace_back(multiname, type.mass(), type.width_at_pole(),
                               type.spin(), type.parity());
    logg[LParticleType].debug()
        << "Creating isospin multiplet " << multiname
        << " [ m = " << type.mass() << ", Γ = " << type.width_at_pole() << " ]";
  }

  // sort the iso-type list by name
  std::sort(iso_type_list.begin(), iso_type_list.end(),
            [](const IsoParticleType &l, const IsoParticleType &r) {
              return l.name() < r.name();
            });

  // add the specific type to the multiplet
  IsoParticleType &multiplet = find_private(multiname);
  multiplet.add_state(type);
}

static Integrator integrate;
static Integrator2dCuhre integrate2d;

/**
 * Tabulation of all N R integrals.
 *
 * Keys are the multiplet names (which are unique).
 */
static std::unordered_map<std::string, Tabulation> NR_tabulations;

/**
 * Tabulation of all pi R integrals.
 *
 * Keys are the multiplet names (which are unique).
 */
static std::unordered_map<std::string, Tabulation> piR_tabulations;

/**
 * Tabulation of all K R integrals.
 *
 * Keys are the multiplet names (which are unique).
 */
static std::unordered_map<std::string, Tabulation> RK_tabulations;

/**
 * Tabulation of all Delta R integrals.
 *
 * Keys are the pairs of multiplet names (which are unique).
 */
static std::unordered_map<std::string, Tabulation> DeltaR_tabulations;

/**
 * Tabulation of all rho rho integrals.
 *
 * Keys are the pairs of multiplet names (which are unique).
 */
static std::unordered_map<std::string, Tabulation> rhoR_tabulations;

static bf::path generate_tabulation_path(const bf::path &dir,
                                         const std::string &prefix,
                                         const std::string &res_name) {
  return dir / (prefix + res_name + ".bin");
}

inline void cache_integral(
    std::unordered_map<std::string, Tabulation> &tabulations,
    const bf::path &dir, sha256::Hash hash, const IsoParticleType &part,
    const IsoParticleType &res, const IsoParticleType *antires, bool unstable) {
  constexpr double spacing = 2.0;
  constexpr double spacing2d = 3.0;
  const auto path = generate_tabulation_path(dir, part.name(), res.name());
  Tabulation integral;
  if (!dir.empty() && bf::exists(path)) {
    std::ifstream file(path.string());
    integral = Tabulation::from_file(file, hash);
    if (!integral.is_empty()) {
      // Only print message if the found tabulation was valid.
      std::cout << "Tabulation found at " << path.filename() << '\r'
                << std::flush;
    }
  }
  if (integral.is_empty()) {
    if (!dir.empty()) {
      std::cout << "Caching tabulation to " << path.filename() << '\r'
                << std::flush;
    }
    if (!unstable) {
      integral = spectral_integral_semistable(integrate, *res.get_states()[0],
                                              *part.get_states()[0], spacing);
    } else {
      integral = spectral_integral_unstable(integrate2d, *res.get_states()[0],
                                            *part.get_states()[0], spacing2d);
    }
    if (!dir.empty()) {
      std::ofstream file(path.string());
      integral.write(file, hash);
    }
  }
  tabulations.emplace(std::make_pair(res.name(), integral));
  if (antires != nullptr) {
    tabulations.emplace(std::make_pair(antires->name(), integral));
  }
}

void IsoParticleType::tabulate_integrals(sha256::Hash hash,
                                         const bf::path &tabulations_path) {
  // To avoid race conditions, make sure we are the only ones currently storing
  // tabulations. Otherwise, we ignore any stored tabulations and don't store
  // our results.
  FileLock lock(tabulations_path / "tabulations.lock");
  const bf::path &dir = lock.acquire() ? tabulations_path : "";

  const auto nuc = IsoParticleType::try_find("N");
  const auto pion = IsoParticleType::try_find("π");
  const auto kaon = IsoParticleType::try_find("K");
  const auto delta = IsoParticleType::try_find("Δ");
  const auto rho = IsoParticleType::try_find("ρ");
  const auto h1 = IsoParticleType::try_find("h₁(1170)");
  for (const auto &res : IsoParticleType::list_baryon_resonances()) {
    const auto antires = res->anti_multiplet();
    if (nuc) {
      cache_integral(NR_tabulations, dir, hash, *nuc, *res, antires, false);
    }
    if (pion) {
      cache_integral(piR_tabulations, dir, hash, *pion, *res, antires, false);
    }
    if (kaon) {
      cache_integral(RK_tabulations, dir, hash, *kaon, *res, antires, false);
    }
    if (delta) {
      cache_integral(DeltaR_tabulations, dir, hash, *delta, *res, antires,
                     true);
    }
  }
  if (rho) {
    cache_integral(rhoR_tabulations, dir, hash, *rho, *rho, nullptr, true);
  }
  if (rho && h1) {
    cache_integral(rhoR_tabulations, dir, hash, *rho, *h1, nullptr, true);
  }
}

double IsoParticleType::get_integral_NR(double sqrts) {
  if (XS_NR_tabulation_ == nullptr) {
    const auto res = states_[0]->iso_multiplet();
    XS_NR_tabulation_ = &NR_tabulations.at(res->name());
  }
  return XS_NR_tabulation_->get_value_linear(sqrts);
}

double IsoParticleType::get_integral_piR(double sqrts) {
  if (XS_piR_tabulation_ == nullptr) {
    const auto res = states_[0]->iso_multiplet();
    XS_piR_tabulation_ = &piR_tabulations.at(res->name());
  }
  return XS_piR_tabulation_->get_value_linear(sqrts);
}

double IsoParticleType::get_integral_RK(double sqrts) {
  if (XS_RK_tabulation_ == nullptr) {
    const auto res = states_[0]->iso_multiplet();
    XS_RK_tabulation_ = &RK_tabulations.at(res->name());
  }
  return XS_RK_tabulation_->get_value_linear(sqrts);
}

double IsoParticleType::get_integral_rhoR(double sqrts) {
  if (XS_rhoR_tabulation_ == nullptr) {
    const auto res = states_[0]->iso_multiplet();
    XS_rhoR_tabulation_ = &rhoR_tabulations.at(res->name());
  }
  return XS_rhoR_tabulation_->get_value_linear(sqrts);
}

double IsoParticleType::get_integral_RR(IsoParticleType *type_res_2,
                                        double sqrts) {
  const auto res = states_[0]->iso_multiplet();
  if (type_res_2->states_[0]->is_Delta()) {
    if (XS_DeltaR_tabulation_ == nullptr) {
      XS_DeltaR_tabulation_ = &DeltaR_tabulations.at(res->name());
    }
    return XS_DeltaR_tabulation_->get_value_linear(sqrts);
  }
  if (type_res_2->name() == "ρ") {
    if (XS_rhoR_tabulation_ == nullptr) {
      XS_rhoR_tabulation_ = &rhoR_tabulations.at(res->name());
    }
    return XS_rhoR_tabulation_->get_value_linear(sqrts);
  }
  if (type_res_2->name() == "h₁(1170)") {
    if (XS_rhoR_tabulation_ == nullptr) {
      XS_rhoR_tabulation_ = &rhoR_tabulations.at(res->name());
    }
    return XS_rhoR_tabulation_->get_value_linear(sqrts);
  }
  std::stringstream err;
  err << "RR=" << name() << type_res_2->name() << " is not implemented";
  throw std::runtime_error(err.str());
}

}  // namespace smash
