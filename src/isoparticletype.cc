/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "include/isoparticletype.h"

#include "include/integrate.h"
#include "include/kinematics.h"
#include "include/logging.h"

namespace Smash {

static IsoParticleTypeList iso_type_list;

IsoParticleType::IsoParticleType(const std::string &n, float m, float w,
                                 unsigned int s)
                                : name_(n), mass_(m), width_(w), spin_(s) {}

const IsoParticleTypeList& IsoParticleType::list_all() { return iso_type_list; }

const IsoParticleType& IsoParticleType::find(const std::string &name) {
  const auto found = std::lower_bound(
      iso_type_list.begin(), iso_type_list.end(), name,
      [](const IsoParticleType &l, const std::string &r) {
        return l.name() < r;
      });
  if (found == iso_type_list.end() || found->name() != name) {
    throw ParticleNotFoundFailure("Isospin multiplet " + name + " not found!");
  }
  return *found;
}

IsoParticleType& IsoParticleType::find_private(const std::string &name) {
  auto found = std::lower_bound(
      iso_type_list.begin(), iso_type_list.end(), name,
      [](const IsoParticleType &l, const std::string &r) {
        return l.name() < r;
      });
  if (found == iso_type_list.end() || found->name() != name) {
    throw std::runtime_error("Isospin multiplet " + name +
                             " not found (privately)!");
  }
  return *found;
}

bool IsoParticleType::exists(const std::string &name) {
  const auto found = std::lower_bound(
      iso_type_list.begin(), iso_type_list.end(), name,
      [](const IsoParticleType &l, const std::string &r) {
        return l.name() < r;
      });
  if (found != iso_type_list.end()) {
    return found->name() == name;
  }
  return false;
}

/* Construct the name-string for an isospin multiplet from the given
 * name-string for the particle. */
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

bool IsoParticleType::has_anti_multiplet() const {
  if (states_[0]->has_antiparticle()) {
    ParticleTypePtr anti = states_[0]->get_antiparticle();
    return multiplet_name(states_[0]->name()) != multiplet_name(anti->name());
  } else {
    return false;
  }
}

const ParticleTypePtr IsoParticleType::find_state(
                                                        const std::string &n) {
  const IsoParticleType &multiplet = IsoParticleType::find(multiplet_name(n));
  auto found = std::find_if(
    multiplet.states_.begin(), multiplet.states_.end(),
    [&n](ParticleTypePtr p) {
      return p->name() == n;
    });
  if (found == multiplet.states_.end() || (*found)->name() != n) {
    throw std::runtime_error("Isospin state " + n + " not found!");
  }
  return *found;
}

IsoParticleType* IsoParticleType::find(const ParticleType &type) {
  std::string multiname = multiplet_name(type.name());
  IsoParticleType &multiplet = find_private(multiname);
  return &multiplet;
}

void IsoParticleType::add_state(const ParticleType &type) {
  states_.push_back(&type);

  // check if isospin symmetry is fulfilled
  const auto &log = logger<LogArea::ParticleType>();
  if (std::abs(mass() - type.mass()) > really_small) {
    log.warn() << "Isospin symmetry is broken by mass of " << type.name()
               << ": " << type.mass() << " vs. " << mass();
  }
  if (std::abs(width() - type.width_at_pole()) > really_small) {
    log.warn() << "Isospin symmetry is broken by width of " << type.name()
               << ": " << type.width_at_pole() << " vs. " << width();
  }
  if (spin() != type.spin()) {
    log.error() << "Isospin symmetry is broken by spin of " << type.name()
                << ": " << type.spin() << " vs. " << spin();
  }
}

void IsoParticleType::create_multiplet(const ParticleType &type) {
  const auto &log = logger<LogArea::ParticleType>();

  // create multiplet if it does not exist yet
  std::string multiname = multiplet_name(type.name());
  if (!exists(multiname)) {
    iso_type_list.emplace_back(multiname, type.mass(), type.width_at_pole(),
                               type.spin());
    log.debug() << "Creating isospin multiplet " << multiname
                << " [ m = " << type.mass()
                << ", Γ = " << type.width_at_pole() << " ]";
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


/**
 * Spectral function integrand for GSL integration, with one resonance in the
 * final state (the second particle is stable).
 *
 * The integrand is \f$ A(m) p_{cm}^f \f$, where \f$ m \f$ is the
 * resonance mass, \f$ A(m) \f$ is the spectral function
 *  and \f$ p_{cm}^f \f$ is the center-of-mass momentum of the final state.
 *
 * \param[in] resonance_mass Actual mass of the resonance.
 * \param[in] sqrts Center-of-mass energy, i.e. sqrt of Mandelstam s.
 * \param[in] stable_mass mass of the stable particle in the final state
 * \param[in] type type of the resonance
 */
static float spec_func_integrand_1res(float resonance_mass, float sqrts,
                               float stable_mass, const ParticleType &type) {
  if (sqrts <= stable_mass + resonance_mass) {
    return 0.;
  }

  /* Integrand is the spectral function weighted by the CM momentum of the
   * final state. */
  return type.spectral_function(resonance_mass)
       * pCM(sqrts, stable_mass, resonance_mass);
}


/**
 * Spectral function integrand for GSL integration, with two resonances in the
 * final state.
 *
 * The integrand is \f$ A_1(m_1) A_2(m_2) p_{cm}^f \f$, where \f$ m_1 \f$ and
 * \f$ m_2 \f$ are the resonance masses, \f$ A_1 \f$ and \f$ A_2 \f$ are the
 * spectral functions and \f$ p_{cm}^f \f$ is the center-of-mass momentum of
 * the final state.
 *
 * \param[in] sqrts Center-of-mass energy, i.e. sqrt of Mandelstam s.
 * \param[in] res_mass_1 Actual mass of the first resonance.
 * \param[in] res_mass_2 Actual mass of the second resonance.
 * \param[in] t1 Type of the first resonance.
 * \param[in] t2 Type of the second resonance.
 */
static float spec_func_integrand_2res(float sqrts,
                              float res_mass_1, float res_mass_2,
                              const ParticleType &t1, const ParticleType &t2) {
  if (sqrts <= res_mass_1 + res_mass_2) {
    return 0.;
  }

  /* Integrand is the product of the spectral function weighted by the
   * CM momentum of the final state. */
  return t1.spectral_function(res_mass_1)
       * t2.spectral_function(res_mass_2)
       * pCM(sqrts, res_mass_1, res_mass_2);
}

static thread_local Integrator integrate;

double IsoParticleType::get_integral_NR(double sqrts) {
  if (XS_NR_tabulation == nullptr) {
    // initialize tabulation
    /* TODO(weil): Move this lazy init to a global initialization function,
      * in order to avoid race conditions in multi-threading. */
    ParticleTypePtr type_res = states_[0];
    ParticleTypePtr nuc = IsoParticleType::find("N").get_states()[0];
    XS_NR_tabulation = make_unique<Tabulation>(
          type_res->minimum_mass() + nuc->mass(), 2.f, 100,
          [&](float srts) {
            return integrate(type_res->minimum_mass(), srts - nuc->mass(),
                             [&](float m) {
                               return spec_func_integrand_1res(m, srts,
                                                        nuc->mass(), *type_res);
                             });
          });
  }
  return XS_NR_tabulation->get_value_linear(sqrts);
}

static thread_local Integrator2d integrate2d(1E4);

double IsoParticleType::get_integral_DR(double sqrts) {
  if (XS_DR_tabulation == nullptr) {
    // initialize tabulation
    /* TODO(weil): Move this lazy init to a global initialization function,
      * in order to avoid race conditions in multi-threading. */
//    ParticleTypePtr type_res = states_[0];
    ParticleTypePtr Delta = IsoParticleType::find("Δ").get_states()[0];
    XS_DR_tabulation = integrate_RR(Delta);
//make_unique<Tabulation>(
//          type_res->minimum_mass() + Delta->minimum_mass(), 2.5f, 100,
//          [&](float srts) {
//            return integrate2d(type_res->minimum_mass(),
//                               srts - Delta->minimum_mass(),
//                               Delta->minimum_mass(),
//                               srts - type_res->minimum_mass(),
//                               [&](float m1, float m2) {
//                                 return spec_func_integrand_2res(srts, m1, m2,
//                                                            *type_res, *Delta);
//                               });
//          });
  }
  return XS_DR_tabulation->get_value_linear(sqrts);
}

double IsoParticleType::get_integral_RR(const ParticleType &type_res_2, double sqrts) {
  for(int i=0;i<XS_RR_tabulations.size(); i++) {
    if (resonances_[i] == find(type_res_2)) {
      std::cout << "WAZZA\n";
      return XS_RR_tabulations[i]->get_value_linear(sqrts);
    }
  }
  std::cout << "Found new resonance pair : " << states_[0]->pdgcode() << " " << type_res_2.pdgcode() << std::endl;
  XS_RR_tabulations.push_back(integrate_RR(find(type_res_2)->get_states()[0]));
  resonances_.push_back(find(type_res_2));
  return XS_RR_tabulations.back()->get_value_linear(sqrts);
}

TabulationPtr IsoParticleType::integrate_RR(ParticleTypePtr &type_res_2) {
  ParticleTypePtr type_res_1 = states_[0];
  return make_unique<Tabulation>(
         type_res_1->minimum_mass() + type_res_2->minimum_mass(), 3.f, 125,
         [&](float srts) {
            return integrate2d(type_res_1->minimum_mass(),
                               srts - type_res_2->minimum_mass(),
                               type_res_2->minimum_mass(),
                               srts - type_res_1->minimum_mass(),
                               [&](float m1, float m2) {
                                  return spec_func_integrand_2res(srts, m1, m2,
                                                      *type_res_1, *type_res_2);
                               });
         });
}

}  // namespace Smash
