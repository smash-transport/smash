/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteractionnucleonnucleon.h"

#include "include/angles.h"
#include "include/cxx14compat.h"
#include "include/integrate.h"
#include "include/kinematics.h"
#include "include/parametrizations.h"
#include "include/random.h"
#include "include/resonances.h"

namespace Smash {


CollisionBranchPtr ScatterActionNucleonNucleon::elastic_cross_section(float) {
  const PdgCode &pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode &pdg_b = incoming_particles_[1].type().pdgcode();

  const double s = mandelstam_s();

  /* Use parametrized cross sections. */
  float sig_el;
  if (pdg_a == pdg_b) {                          /* pp */
    sig_el = pp_elastic(s);
  } else if (pdg_a.is_antiparticle_of(pdg_b)) {  /* ppbar */
    sig_el = ppbar_elastic(s);
  } else {                                       /* np */
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
}


/**
 * Computes the B coefficients from the Cugnon parametrization of the angular
 * distribution in elastic pp scattering, see equation (8) in \iref{Cugnon:1996kh}.
 * Note: The original Cugnon parametrization is only applicable for
 * plab < 6 GeV and keeps rising above that. We add an upper limit of b <= 9,
 * in order to be compatible with high-energy data (up to plab ~ 25 GeV).
 * \param[in] plab Lab momentum in GeV.
 */
static float Cugnon_bpp(float plab) {
  if (plab < 2.) {
    float p8 = std::pow(plab, 8);
    return 5.5*p8 / (7.7+p8);
  } else {
    return std::min(9.0, 5.334 + 0.67*(plab-2.));
  }
}


/**
 * Computes the B coefficients from the Cugnon parametrization of the angular
 * distribution in elastic np scattering, see equation (10) in \iref{Cugnon:1996kh}.
 * \param[in] plab Lab momentum in GeV.
 */
static float Cugnon_bnp(float plab) {
  if (plab < 0.225) {
    return 0.;
  } else if (plab < 0.6) {
    return 16.53*(plab-0.225);
  } else if (plab < 1.6) {
    return -1.63*plab + 7.16;
  } else {
    return Cugnon_bpp(plab);
  }
}


void ScatterActionNucleonNucleon::elastic_scattering() {
  const auto &log = logger<LogArea::ScatterAction>();
  outgoing_particles_[0] = incoming_particles_[0];
  outgoing_particles_[1] = incoming_particles_[1];

  /* Determine absolute momentum in center-of-mass frame. */
  const double momentum_radial = cm_momentum();

  Angles phitheta;
  if (isotropic_) {
    phitheta.distribute_isotropically();
  } else {
    /** Choose angular distribution according to Cugnon parametrization,
     * see \iref{Cugnon:1996kh}. */
    double bb, a, plab = plab_from_s_NN(mandelstam_s());
    if (incoming_particles_[0].type().charge() +
        incoming_particles_[1].type().charge() == 1) {
      // pn
      bb = std::max(Cugnon_bnp(plab), really_small);
      a = (plab < 0.8) ? 1. : 0.64/(plab*plab);
    } else {
      // pp or nn
      bb = std::max(Cugnon_bpp(plab), really_small);
      a = 1.;
    }
    double t, t_max = -4.*momentum_radial*momentum_radial;
    if (Random::canonical() <= 1./(1.+a)) {
      t = Random::expo(bb, 0., t_max);
    } else {
      t = t_max - Random::expo(bb, 0., t_max);
    }

    // determine scattering angles in center-of-mass frame
    phitheta = Angles(2.*M_PI*Random::canonical(), 1. - 2.*t/t_max);
  }

  ThreeVector pscatt = phitheta.threevec();
  // 3-momentum of first incoming particle in center-of-mass frame
  ThreeVector pcm = incoming_particles_[0].momentum().
                    LorentzBoost(beta_cm()).threevec();

  pscatt.rotate_to(pcm);

  /* Set 4-momentum: Masses stay the same, 3-momentum changes. */
  outgoing_particles_[0].set_4momentum(outgoing_particles_[0].effective_mass(),
                                       pscatt * momentum_radial);
  outgoing_particles_[1].set_4momentum(outgoing_particles_[1].effective_mass(),
                                       - pscatt * momentum_radial);

  /* debug output */
  log.debug("exchanged momenta a", outgoing_particles_[0].momentum());
  log.debug("exchanged momenta b", outgoing_particles_[1].momentum());
}


CollisionBranchList ScatterActionNucleonNucleon::two_to_two_cross_sections() {
  const ParticleType &type_particle_a = incoming_particles_[0].type();
  const ParticleType &type_particle_b = incoming_particles_[1].type();

  /* Find all resonance-production channels. */
  CollisionBranchList process_list = nuc_nuc_to_nuc_res(type_particle_a,
                                                      type_particle_b);

  return process_list;
}


CollisionBranchList ScatterActionNucleonNucleon::nuc_nuc_to_nuc_res(
                            const ParticleType &type_particle_a,
                            const ParticleType &type_particle_b) {
  const auto &log = logger<LogArea::ScatterAction>();
  CollisionBranchList process_list;
  const double s = mandelstam_s();
  const double srts = sqrt_s();

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
      double upper_limit = srts - second_type->mass();
      /* Check the available energy (requiring it to be a little above the
      * threshold, because the integration will not work if it's too close). */
      if (upper_limit - lower_limit < 1E-3) {
        continue;
      }

      /* Calculate matrix element. */
      const float matrix_element =
            nn_to_resonance_matrix_element(srts, *type_resonance, *second_type);
      if (matrix_element <= 0.) {
        continue;
      }

      /* Calculate resonance production cross section
       * using the Breit-Wigner distribution as probability amplitude.
       * Integrate over the allowed resonance mass range. */

      const int res_id = type_resonance->pdgcode().iso_multiplet();
      if (XS_tabulation[res_id] == NULL) {
        // initialize tabulation, we need one per resonance multiplet
        /* TODO(weil): Move this lazy init to a global initialization function,
         * in order to avoid race conditions in multi-threading. */
        Integrator integrate;
        XS_tabulation[res_id] = make_unique<Tabulation>(
              type_resonance->minimum_mass() + second_type->mass(), 2.f, 100,
              [&](float sqrts) {
                return integrate(type_resonance->minimum_mass(),
                                  sqrts - second_type->mass(),
                                  [&](float m) {
                                    return spectral_function_integrand(m, sqrts,
                                                            second_type->mass(),
                                                            *type_resonance);
                                  });
              });
      }
      const double resonance_integral =
                    XS_tabulation[res_id]->get_value_linear(srts);

      /** Cross section for 2->2 process with one resonance in final state.
       * Based on Eq. (46) in \iref{Weil:2013mya}. */
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


void ScatterActionNucleonNucleon::sample_cms_momenta() {
  const auto &log = logger<LogArea::ScatterAction>();
  /* This function only operates on 2-particle final states. */
  assert(outgoing_particles_.size() == 2);

  ParticleData *p_a = &outgoing_particles_[0];
  ParticleData *p_b = &outgoing_particles_[1];

  const ParticleType &t_a = p_a->type();
  const ParticleType &t_b = p_b->type();

  double mass_a = t_a.mass();
  double mass_b = t_b.mass();

  const double cms_energy = sqrt_s();

  if (cms_energy < t_a.minimum_mass() + t_b.minimum_mass()) {
    throw InvalidResonanceFormation("resonance_formation: not enough energy! " +
      std::to_string(cms_energy) + " " + std::to_string(t_a.minimum_mass()) +
      " " + std::to_string(t_b.minimum_mass()) + " " +
      p_a->pdgcode().string() + " " + p_b->pdgcode().string());
  }

  /* If one of the particles is a resonance, sample its mass. */
  /* TODO: Other particle assumed stable! */
  if (!t_a.is_stable()) {
    mass_a = sample_resonance_mass(t_a, t_b.mass(), cms_energy);
  } else if (!t_b.is_stable()) {
    mass_b = sample_resonance_mass(t_b, t_a.mass(), cms_energy);
  }

  const std::array<double, 2> t_range = get_t_range(cms_energy, mN, mN,
                                                    mass_a, mass_b);
  Angles phitheta;
  if (t_a.pdgcode().iso_multiplet() == 0x1114 &&
      t_b.pdgcode().iso_multiplet() == 0x1112 && !isotropic_) {
    /** NN->NDelta: Sample scattering angles in center-of-mass frame from an
     * anisotropic angular distribution, using the same distribution as for
     * elastic pp scattering, as suggested in \iref{Cugnon:1996kh}. */
    double plab = plab_from_s_NN(mandelstam_s());
    double bb = std::max(Cugnon_bpp(plab), really_small);
    double t = Random::expo(bb, t_range[0], t_range[1]);
    if (Random::canonical() > 0.5) {
      t = t_range[0] + t_range[1] - t;  // symmetrize
    }
    phitheta = Angles(2.*M_PI*Random::canonical(),
                      1. - 2.*(t-t_range[0])/(t_range[1]-t_range[0]));
  } else if (t_b.pdgcode().iso_multiplet() == 0x1112 && !isotropic_ &&
             (t_a.pdgcode().is_Nstar() || t_a.pdgcode().is_Deltastar())) {
    /** NN->NR: Fit to HADES data, see \iref{Agakishiev:2014wqa}. */
    const std::array<float, 4> p { 1.46434, 5.80311, -6.89358, 1.94302 };
    const double a = p[0] + mass_a * (p[1] + mass_a * (p[2] + mass_a * p[3]));
    double t = Random::power(-a, t_range[0], t_range[1]);
    if (Random::canonical() > 0.5) {
      t = t_range[0] + t_range[1] - t;  // symmetrize
    }
    phitheta = Angles(2.*M_PI*Random::canonical(),
                      1. - 2.*(t-t_range[0])/(t_range[1]-t_range[0]));
  } else {
    /* isotropic angular distribution */
    phitheta.distribute_isotropically();
  }

  ThreeVector pscatt = phitheta.threevec();
  // 3-momentum of first incoming particle in center-of-mass frame
  ThreeVector pcm = incoming_particles_[0].momentum().
                    LorentzBoost(beta_cm()).threevec();
  pscatt.rotate_to(pcm);

  // final-state CM momentum
  const double p_f = pCM(cms_energy, mass_a, mass_b);
  if (!(p_f > 0.0)) {
    log.warn("Particle: ", t_a.pdgcode(),
             " radial momentum: ", p_f);
    log.warn("Etot: ", cms_energy, " m_a: ", mass_a, " m_b: ", mass_b);
  }
  p_a->set_4momentum(mass_a,  pscatt * p_f);
  p_b->set_4momentum(mass_b, -pscatt * p_f);

  log.debug("p_a: ", *p_a, "\np_b: ", *p_b);
}


}  // namespace Smash
