/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteractionsfinder.h"

#include <algorithm>

#include "include/configuration.h"
#include "include/constants.h"
#include "include/cxx14compat.h"
#include "include/decaymodes.h"
#include "include/experimentparameters.h"
#include "include/isoparticletype.h"
#include "include/logging.h"
#include "include/macros.h"
#include "include/particles.h"
#include "include/scatteraction.h"
#include "include/scatteractionbaryonbaryon.h"
#include "include/scatteractionbaryonmeson.h"
#include "include/scatteractiondeltakaon.h"
#include "include/scatteractionhyperonpion.h"
#include "include/scatteractionmesonmeson.h"
#include "include/scatteractionnucleonkaon.h"
#include "include/scatteractionnucleonnucleon.h"
#include "include/scatteractionphoton.h"
#include "include/stringfunctions.h"

namespace Smash {
/*!\Userguide
* \page input_collision_term_ Collision_Term
* \key Elastic_Cross_Section (float, optional, default = -1.0 [mb]) \n
* If a non-negative value is given, it will override the parametrized
* elastic cross sections (which are energy-dependent) with a constant value.
* This constant elastic cross section is used for all collisions.
*
* \key Isotropic (bool, optional, default = false) \n
* Do all collisions isotropically.
* \key Strings (bool, optional, default = false): \n
* true - string excitation is enabled\n
* false - string excitation is disabled
* \key Formation_Time (float, optional, default = 1.0) \n
* Parameter for formation time in string fragmentation in fm/c
* \key low_snn_cut (double) in GeV \n
* The elastic collisions betwen two nucleons with sqrt_s below
* low_snn_cut cannot happen.
* <1.88 - below the threshold energy of the elastic collsion, no effect
* >2.02 - beyond the threshold energy of the inelastic collision NN->NNpi, not suggested
*/

ScatterActionsFinder::ScatterActionsFinder(
    Configuration config, const ExperimentParameters &parameters,
    const std::vector<bool> &nucleon_has_interacted, int N_tot, int N_proj,
     int n_fractional_photons = 1)
    : elastic_parameter_(config.take({"Collision_Term",
                                      "Elastic_Cross_Section"}, -1.0f)),
      testparticles_(parameters.testparticles),
      isotropic_(config.take({"Collision_Term", "Isotropic"}, false)),
      two_to_one_(parameters.two_to_one),
      two_to_two_(parameters.two_to_two),
      low_snn_cut_(parameters.low_snn_cut),
      strings_switch_(parameters.strings_switch),
      nucleon_has_interacted_(nucleon_has_interacted),
      N_tot_(N_tot),
      N_proj_(N_proj),
      formation_time_(config.take({"Collision_Term", "Formation_Time"}, 1.0f)),
      photons_(parameters.photons_switch),
      n_fractional_photons_(n_fractional_photons) {
        if (is_constant_elastic_isotropic()) {
          const auto &log = logger<LogArea::FindScatter>();
          log.info("Constant elastic isotropic cross-section mode:",
          " using ", elastic_parameter_, " mb as maximal cross-section.");
        }
      }

ScatterActionsFinder::ScatterActionsFinder(
    float elastic_parameter, int testparticles,
    const std::vector<bool> &nucleon_has_interacted, bool two_to_one)
    : elastic_parameter_(elastic_parameter),
      testparticles_(testparticles),
      isotropic_(false),
      two_to_one_(two_to_one),
      two_to_two_(true),
      low_snn_cut_(0.0),
      strings_switch_(true),
      nucleon_has_interacted_(nucleon_has_interacted),
      N_tot_(0),
      N_proj_(0),
      formation_time_(1.0f),
      photons_(false),
      n_fractional_photons_(1) {}

ScatterActionPtr ScatterActionsFinder::construct_scatter_action(
                                            const ParticleData &data_a,
                                            const ParticleData &data_b,
                                            double time_until_collision)
                                            const {
  const auto &pdg_a = data_a.pdgcode();
  const auto &pdg_b = data_b.pdgcode();
  ScatterActionPtr act;
  if (data_a.is_baryon() && data_b.is_baryon()) {
    if ((pdg_a.is_nucleon() && pdg_b.is_nucleon()) &&
        (pdg_a.antiparticle_sign() == pdg_b.antiparticle_sign())) {
      act = make_unique<ScatterActionNucleonNucleon>(data_a, data_b,
                                              time_until_collision, isotropic_,
                                              formation_time_);
    } else {
      act = make_unique<ScatterActionBaryonBaryon>(data_a, data_b,
                                              time_until_collision, isotropic_,
                                              formation_time_);
    }
  } else if (data_a.is_baryon() || data_b.is_baryon()) {
    if ((pdg_a.is_nucleon() && pdg_b.is_kaon()) ||
        (pdg_b.is_nucleon() && pdg_a.is_kaon())) {
      act = make_unique<ScatterActionNucleonKaon>(data_a, data_b,
                                              time_until_collision, isotropic_,
                                              formation_time_);
    } else if ((pdg_a.is_hyperon() && pdg_b.is_pion()) ||
               (pdg_b.is_hyperon() && pdg_a.is_pion())) {
      act = make_unique<ScatterActionHyperonPion>(data_a, data_b,
                                              time_until_collision, isotropic_,
                                              formation_time_);
    } else if ((pdg_a.is_Delta() && pdg_b.is_kaon()) ||
               (pdg_b.is_Delta() && pdg_a.is_kaon())) {
      act = make_unique<ScatterActionDeltaKaon>(data_a, data_b,
                                              time_until_collision, isotropic_,
                                              formation_time_);
    } else {
      act = make_unique<ScatterActionBaryonMeson>(data_a, data_b,
                                              time_until_collision, isotropic_,
                                              formation_time_);
    }
  } else {
    act = make_unique<ScatterActionMesonMeson>(data_a, data_b,
                                              time_until_collision, isotropic_,
                                              formation_time_);
  }
  return act;
}

ActionPtr ScatterActionsFinder::check_collision(
    const ParticleData &data_a, const ParticleData &data_b, double dt) const {
#ifndef NDEBUG
  const auto &log = logger<LogArea::FindScatter>();
#endif

  /* just collided with this particle */
  if (data_a.id_process() > 0 && data_a.id_process() == data_b.id_process()) {
#ifndef NDEBUG
    log.debug("Skipping collided particles at time ", data_a.position().x0(),
              " due to process ", data_a.id_process(),
              "\n    ", data_a,
              "\n<-> ", data_b);
#endif
    return nullptr;
  }
  /** If the two particles
    * 1) belong to the two colliding nuclei
    * 2) are within the same nucleus
    * 3) both of them have never experienced any collisons,
    * then the collision between them are banned. */
  assert(data_a.id() >= 0);
  assert(data_b.id() >= 0);
  if (data_a.id() < N_tot_ && data_b.id() < N_tot_ &&
      ((data_a.id() < N_proj_ && data_b.id() < N_proj_) ||
       (data_a.id() > N_proj_ && data_b.id() > N_proj_)) &&
       !(nucleon_has_interacted_[data_a.id()] ||
         nucleon_has_interacted_[data_b.id()])) {
    return nullptr;
  }

  /* Determine time of collision. */
  const double time_until_collision = collision_time(data_a, data_b);

  /* Check that collision happens in this timestep. */
  if (time_until_collision < 0.f || time_until_collision >= dt) {
    return nullptr;
  }

  /* Create ScatterAction object. */
  ScatterActionPtr act = construct_scatter_action(data_a, data_b,
                                                  time_until_collision);
  const double distance_squared = act->transverse_distance_sqr();

  /* Don't calculate cross section if the particles are very far apart. */
  if (distance_squared >= max_transverse_distance_sqr(testparticles_)) {
    return nullptr;
  }

  /* Add various subprocesses.  */
  act->add_all_processes(elastic_parameter_, two_to_one_,
                         two_to_two_, low_snn_cut_, strings_switch_);

  /* Add photons to collision finding if necessary */
  double photon_cross_section = 0.0;
  if (photons_ &&
    (ScatterActionPhoton::is_photon_reaction(act->incoming_particles())
      != ScatterActionPhoton::ReactionType::no_reaction) ) {
        ScatterActionPhoton photon_act(act->incoming_particles(), 0.0,
                                   n_fractional_photons_);
        photon_act.add_single_channel();
        photon_cross_section = photon_act.cross_section();
  }
  /* Cross section for collision criterion */
  float cross_section_criterion = (act->cross_section() + photon_cross_section)
                                  * fm2_mb * M_1_PI
                                  / static_cast<float>(testparticles_);
  /* Consider cross section scaling factors only if the particles
   * are not formed yet at the prospective time of the interaction */
  if (data_a.formation_time() > data_a.position().x0() + time_until_collision) {
    cross_section_criterion *= data_a.cross_section_scaling_factor();
  }
  if (data_b.formation_time() > data_b.position().x0() + time_until_collision) {
    cross_section_criterion *= data_b.cross_section_scaling_factor();
  }

  /* distance criterion according to cross_section */
  if (distance_squared >= cross_section_criterion) {
    return nullptr;
  }

#ifndef NDEBUG
  log.debug("particle distance squared: ", distance_squared,
            "\n    ", data_a,
            "\n<-> ", data_b);
#endif

  return std::move(act);
}

ActionList ScatterActionsFinder::find_actions_in_cell(
    const ParticleList &search_list, float dt) const {
  std::vector<ActionPtr> actions;
  for (const ParticleData &p1 : search_list) {
    for (const ParticleData &p2 : search_list) {
      if (p1.id() < p2.id()) {
        // Check if a collision is possible.
        ActionPtr act = check_collision(p1, p2, dt);
        if (act) {
          actions.push_back(std::move(act));
        }
      }
    }
  }
  return actions;
}

ActionList ScatterActionsFinder::find_actions_with_neighbors(
    const ParticleList &search_list, const ParticleList &neighbors_list,
    float dt) const {
  std::vector<ActionPtr> actions;
  for (const ParticleData &p1 : search_list) {
    for (const ParticleData &p2 : neighbors_list) {
      assert(p1.id() != p2.id());
      // Check if a collision is possible.
      ActionPtr act = check_collision(p1, p2, dt);
      if (act) {
        actions.push_back(std::move(act));
      }
    }
  }
  return actions;
}

ActionList ScatterActionsFinder::find_actions_with_surrounding_particles(
    const ParticleList &search_list, const Particles &surrounding_list,
    float dt) const {
  std::vector<ActionPtr> actions;
  for (const ParticleData &p2 : surrounding_list) {
    // don't look for collisions if the particle from the surrounding list is
    // also in the search list
    auto result = std::find_if(
        search_list.begin(), search_list.end(),
        [&p2](const ParticleData &p) { return p.id() == p2.id(); });
    if (result != search_list.end()) {
      continue;
    }
    for (const ParticleData &p1 : search_list) {
      // Check if a collision is possible.
      ActionPtr act = check_collision(p1, p2, dt);
      if (act) {
        actions.push_back(std::move(act));
      }
    }
  }
  return actions;
}

void ScatterActionsFinder::dump_reactions() const {
  constexpr double time = 0.0;

  const size_t N_isotypes = IsoParticleType::list_all().size();
  const size_t N_pairs = N_isotypes * (N_isotypes - 1) / 2;

  std::cout << N_isotypes << " iso-particle types." << std::endl;
  std::cout << "They can make " << N_pairs << " pairs." << std::endl;
  std::vector<double> momentum_scan_list = {0.1, 0.3, 0.5, 1.0, 2.0,
                                            3.0, 5.0, 10.0};
  for (const IsoParticleType &A_isotype : IsoParticleType::list_all()) {
    for (const IsoParticleType &B_isotype : IsoParticleType::list_all()) {
      if (&A_isotype > &B_isotype) {
        continue;
      }
      bool any_nonzero_cs = false;
      std::vector<std::string> r_list;
      for (const ParticleTypePtr A_type : A_isotype.get_states()) {
        for (const ParticleTypePtr B_type : B_isotype.get_states()) {
          if (A_type > B_type) {
            continue;
          }
          ParticleData A(*A_type), B(*B_type);
          for (auto mom : momentum_scan_list) {
            A.set_4momentum(A.pole_mass(), mom, 0.0, 0.0);
            B.set_4momentum(B.pole_mass(), -mom, 0.0, 0.0);
            ScatterActionPtr act = construct_scatter_action(A, B, time);
            act->add_all_processes(elastic_parameter_, two_to_one_,
                                   two_to_two_, low_snn_cut_, strings_switch_);
            const float total_cs = act->cross_section();
            if (total_cs <= 0.0) {
              continue;
            }
            any_nonzero_cs = true;
            for (const auto& channel : act->collision_channels()) {
              std::string r;
              if (channel->get_type() == ProcessType::String) {
                r =  A_type->name() + B_type->name()
                     + std::string(" → strings");
              } else {
                std::string r_type =
                  (channel->get_type() == ProcessType::Elastic) ?
                  std::string(" (el)") :
                        (channel->get_type() == ProcessType::TwoToTwo) ?
                        std::string(" (inel)") :
                             std::string(" (?)");
                r = A_type->name() + B_type->name()
                      + std::string(" → ")
                      + channel->particle_types()[0]->name()
                      + channel->particle_types()[1]->name()
                      + r_type;
              }
              isoclean(r);
              r_list.push_back(r);
            }
          }
        }
      }
      std::sort(r_list.begin(), r_list.end());
      r_list.erase(std::unique(r_list.begin(), r_list.end()), r_list.end() );
      if (any_nonzero_cs) {
        for (auto r : r_list) {
          std::cout << r;
          if (r_list.back() != r) {
            std::cout <<  ", ";
          }
        }
        std::cout << std::endl;
      }
    }
  }
}

void ScatterActionsFinder::dump_cs(const ParticleType &a,
                                   const ParticleType &b,
                                   float m_a,
                                   float m_b) const {
  const ParticleTypePtrList l = {&a, &b};
  std::vector<ParticleTypePtr> ab_products;

  for (const ParticleType &resonance : ParticleType::list_all()) {
    const auto &decaymodes = resonance.decay_modes().decay_mode_list();
    for (const auto &mode : decaymodes) {
      if (mode->type().has_particles(l)) {
        ab_products.push_back(&resonance);
      }
    }
  }

  std::string description = "# sqrt(s) [GeV]";
  std::string ab_string = " " + a.name() + b.name() + "->";
  for (const ParticleTypePtr resonance : ab_products) {
    description += ab_string + resonance->name();
  }
  std::cout << description << std::endl;

  if (ab_products.size() > 0) {
    ParticleData a_data(a), b_data(b);
    ScatterAction act(a_data, b_data, 0.0, false, 0.0);
    constexpr int n_points = 1000;

    std::cout << std::fixed;
    std::cout << std::setprecision(8);
    for (int i = 1; i < n_points; i++) {
      const float sqrts = m_a + m_b + 0.01f * i;
      std::cout << sqrts << " ";
      for (const ParticleTypePtr resonance : ab_products) {
        const double p_cm_sqr = pCM_sqr(sqrts, m_a, m_b);
        const double cs = (sqrts < resonance->minimum_mass()) ? 0.0 :
            act.two_to_one_formation(*resonance, sqrts, p_cm_sqr);
        std::cout << cs << " ";
      }
      std::cout << std::endl;
    }
  }
}

}  // namespace Smash
