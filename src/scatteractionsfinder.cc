/*
 *
 *    Copyright (c) 2014-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/scatteractionsfinder.h"

#include <algorithm>
#include <map>
#include <vector>

#include "smash/constants.h"
#include "smash/decaymodes.h"
#include "smash/logging.h"
#include "smash/parametrizations.h"
#include "smash/scatteraction.h"
#include "smash/scatteractionmulti.h"
#include "smash/scatteractionphoton.h"
#include "smash/stringfunctions.h"

namespace smash {
static constexpr int LFindScatter = LogArea::FindScatter::id;
/*!\Userguide
 * \page doxypage_input_conf_ct_string_parameters
 *
 *
 *
 */

ScatterActionsFinder::ScatterActionsFinder(
    Configuration& config, const ExperimentParameters& parameters)
    : finder_parameters_(config, parameters),
      isotropic_(config.take(InputKeys::collTerm_isotropic)),
      box_length_(parameters.box_length),
      string_formation_time_(
          config.take(InputKeys::collTerm_stringParam_formationTime)) {
  if (is_constant_elastic_isotropic()) {
    logg[LFindScatter].info(
        "Constant elastic isotropic cross-section mode:", " using ",
        finder_parameters_.elastic_parameter, " mb as maximal cross-section.");
  }
  if (finder_parameters_.included_multi.any() &&
      finder_parameters_.coll_crit != CollisionCriterion::Stochastic) {
    throw std::invalid_argument(
        "Multi-body reactions (like e.g. 3->1 or 3->2) are only possible with "
        "the stochastic "
        "collision "
        "criterion. Change your config accordingly.");
  }

  if (finder_parameters_
              .included_multi[IncludedMultiParticleReactions::Deuteron_3to2] ==
          1 &&
      (finder_parameters_
               .included_2to2[IncludedReactions::PiDeuteron_to_pidprime] == 1 ||
       finder_parameters_
               .included_2to2[IncludedReactions::NDeuteron_to_Ndprime] == 1)) {
    throw std::invalid_argument(
        "To prevent double counting it is not possible to enable deuteron 3->2 "
        "reactions\nand reactions involving the d' at the same time\ni.e. to "
        "include \"Deuteron_3to2\" in `Multi_Particle_Reactions` and\n "
        "\"PiDeuteron_to_pidprime\" "
        "or \"NDeuteron_to_Ndprime\" in `Included_2to2` at the same time.\n"
        "Change your config accordingly.");
  }

  if (finder_parameters_
              .included_multi[IncludedMultiParticleReactions::Deuteron_3to2] ==
          1 &&
      ParticleType::try_find(pdg::dprime)) {
    throw std::invalid_argument(
        "Do not use the d' resonance and enable \"Deuteron_3to2\" "
        "`Multi_Particle_Reactions` at the same time. Either use the direct "
        "3-to-2 reactions or the d' together with \"PiDeuteron_to_pidprime\" "
        "and \"NDeuteron_to_Ndprime\" in `Included_2to2`. Otherwise the "
        "deuteron 3-to-2 reactions would be double counted.");
  }

  if ((finder_parameters_.nnbar_treatment == NNbarTreatment::TwoToFive &&
       finder_parameters_
               .included_multi[IncludedMultiParticleReactions::NNbar_5to2] !=
           1) ||
      (finder_parameters_
               .included_multi[IncludedMultiParticleReactions::NNbar_5to2] ==
           1 &&
       finder_parameters_.nnbar_treatment != NNbarTreatment::TwoToFive)) {
    throw std::invalid_argument(
        "In order to conserve detailed balance, when \"NNbar_5to2\" is "
        "included in\n`Multi_Particle_Reactions`, the `NNbarTreatment` has to "
        "be set to \"two to five\" and vice versa.");
  }

  if (finder_parameters_.nnbar_treatment == NNbarTreatment::Resonances &&
      finder_parameters_.included_2to2[IncludedReactions::NNbar] != 1) {
    throw std::invalid_argument(
        "'NNbar' has to be in the list of allowed 2 to 2 processes "
        "to enable annihilation to go through resonances");
  }

  if (finder_parameters_.strings_switch) {
    string_process_interface_ = std::make_unique<StringProcess>(
        config.take(InputKeys::collTerm_stringParam_stringTension),
        string_formation_time_,
        config.take(InputKeys::collTerm_stringParam_gluonBeta),
        config.take(InputKeys::collTerm_stringParam_gluonPMin),
        config.take(InputKeys::collTerm_stringParam_quarkAlpha),
        config.take(InputKeys::collTerm_stringParam_quarkBeta),
        config.take(InputKeys::collTerm_stringParam_strangeSuppression),
        config.take(InputKeys::collTerm_stringParam_diquarkSuppression),
        config.take(InputKeys::collTerm_stringParam_sigmaPerp),
        config.take(InputKeys::collTerm_stringParam_stringZALeading),
        config.take(InputKeys::collTerm_stringParam_stringZBLeading),
        config.take(InputKeys::collTerm_stringParam_stringZA),
        config.take(InputKeys::collTerm_stringParam_stringZB),
        config.take(InputKeys::collTerm_stringParam_stringSigmaT),
        config.take(InputKeys::collTerm_stringParam_formTimeFactor),
        config.take(InputKeys::collTerm_stringParam_mDependentFormationTimes),
        config.take(InputKeys::collTerm_stringParam_probabilityPToDUU),
        config.take(InputKeys::collTerm_stringParam_separateFragmentBaryon),
        config.take(InputKeys::collTerm_stringParam_popcornRate),
        config.take(InputKeys::collTerm_stringParam_useMonashTune,
                    parameters.use_monash_tune_default.value()));
  }
}

static StringTransitionParameters create_string_transition_parameters(
    Configuration& config) {
  auto sqrts_range_Npi = config.take(InputKeys::collTerm_stringTrans_rangeNpi);
  auto sqrts_range_NN = config.take(InputKeys::collTerm_stringTrans_rangeNN);

  if (sqrts_range_Npi.first < nucleon_mass + pion_mass) {
    sqrts_range_Npi.first = nucleon_mass + pion_mass;
    if (sqrts_range_Npi.second < sqrts_range_Npi.first)
      sqrts_range_Npi.second = sqrts_range_Npi.first;
    logg[LFindScatter].warn(
        "Lower bound of Sqrts_Range_Npi too small, setting it to mass "
        "threshold. New range is [",
        sqrts_range_Npi.first, ',', sqrts_range_Npi.second, "] GeV");
  }
  if (sqrts_range_NN.first < 2 * nucleon_mass) {
    sqrts_range_NN.first = 2 * nucleon_mass;
    if (sqrts_range_NN.second < sqrts_range_NN.first)
      sqrts_range_NN.second = sqrts_range_NN.first;
    logg[LFindScatter].warn(
        "Lower bound of Sqrts_Range_NN too small, setting it to mass "
        "threshold. New range is [",
        sqrts_range_NN.first, ',', sqrts_range_NN.second, "] GeV.");
  }

  return {sqrts_range_Npi,
          sqrts_range_NN,
          config.take(InputKeys::collTerm_stringTrans_lower),
          config.take(InputKeys::collTerm_stringTrans_range_width),
          config.take(InputKeys::collTerm_stringTrans_pipiOffset),
          config.take(InputKeys::collTerm_stringTrans_KNOffset)};
}

ScatterActionsFinderParameters::ScatterActionsFinderParameters(
    Configuration& config, const ExperimentParameters& parameters)
    : elastic_parameter(config.take(InputKeys::collTerm_elasticCrossSection)),
      low_snn_cut(parameters.low_snn_cut),
      scale_xs(parameters.scale_xs),
      additional_el_xs(
          config.take(InputKeys::collTerm_additionalElasticCrossSection)),
      maximum_cross_section(parameters.maximum_cross_section),
      coll_crit(parameters.coll_crit),
      nnbar_treatment(parameters.nnbar_treatment),
      included_2to2(parameters.included_2to2),
      included_multi(parameters.included_multi),
      testparticles(parameters.testparticles),
      two_to_one(parameters.two_to_one),
      allow_collisions_within_nucleus(
          config.take(InputKeys::modi_collider_collisionWithinNucleus)),
      spin_interaction_type(parameters.spin_interaction_type),
      strings_switch(parameters.strings_switch),
      use_AQM(config.take(InputKeys::collTerm_useAQM)),
      strings_with_probability(
          config.take(InputKeys::collTerm_stringsWithProbability)),
      only_warn_for_high_prob(
          config.take(InputKeys::collTerm_onlyWarnForHighProbability)),
      transition_high_energy{create_string_transition_parameters(config)},
      total_xs_strategy(config.take(InputKeys::collTerm_totXsStrategy)),
      pseudoresonance_method(config.take(InputKeys::collTerm_pseudoresonance)),
      AQM_charm_suppression(
          config.take(InputKeys::collTerm_HF_AQMcSuppression)),
      AQM_bottom_suppression(
          config.take(InputKeys::collTerm_HF_AQMbSuppression)) {
  if (total_xs_strategy == TotalCrossSectionStrategy::BottomUp) {
    logg[LFindScatter].info(
        "Evaluating total cross sections from partial processes.");
  } else if (parameters.included_2to2[IncludedReactions::Elastic] == 1 &&
             parameters.included_2to2.count() == 1) {
    throw std::invalid_argument(
        "The BottomUp strategy for total cross section evaluation is needed to "
        "have only elastic interactions, please change the configuration "
        "accordingly.");
  } else if (total_xs_strategy == TotalCrossSectionStrategy::TopDown) {
    logg[LFindScatter].info(
        "Evaluating total cross sections from parametrizations.");
  } else if (total_xs_strategy == TotalCrossSectionStrategy::TopDownMeasured) {
    logg[LFindScatter].info(
        "Evaluating total cross sections from parametrizations only for "
        "measured processes.");
  }

  if (AQM_charm_suppression < 0 || AQM_bottom_suppression < 0 ||
      AQM_charm_suppression > 1 || AQM_bottom_suppression > 1) {
    throw std::invalid_argument(
        "Suppression factors for AQM should be between 0 and 1.");
  }
}

ActionPtr ScatterActionsFinder::check_collision_two_part(
    const ParticleData& data_a, const ParticleData& data_b, double dt,
    const std::vector<FourVector>& beam_momentum,
    const double gcell_vol) const {
  /* If the two particles
   * 1) belong to one of the two colliding nuclei, and
   * 2) both of them have never experienced any collisions,
   * then the collisions between them are banned. */
  if (!finder_parameters_.allow_collisions_within_nucleus) {
    assert(data_a.id() >= 0);
    assert(data_b.id() >= 0);
    bool in_same_nucleus = (data_a.belongs_to() == BelongsTo::Projectile &&
                            data_b.belongs_to() == BelongsTo::Projectile) ||
                           (data_a.belongs_to() == BelongsTo::Target &&
                            data_b.belongs_to() == BelongsTo::Target);
    bool never_interacted_before =
        data_a.get_history().collisions_per_particle == 0 &&
        data_b.get_history().collisions_per_particle == 0;
    if (in_same_nucleus && never_interacted_before) {
      return nullptr;
    }
  }

  // No grid or search in cell means no collision for stochastic criterion
  if (finder_parameters_.coll_crit == CollisionCriterion::Stochastic &&
      gcell_vol < really_small) {
    return nullptr;
  }

  // Determine time of collision.
  const double time_until_collision =
      collision_time(data_a, data_b, dt, beam_momentum);

  // Check that collision happens in this timestep.
  if (time_until_collision < 0. || time_until_collision >= dt) {
    return nullptr;
  }

  // Determine which total cross section to use
  bool incoming_parametrized = (finder_parameters_.total_xs_strategy ==
                                TotalCrossSectionStrategy::TopDown);
  if (finder_parameters_.total_xs_strategy ==
      TotalCrossSectionStrategy::TopDownMeasured) {
    const PdgCode& pdg_a = data_a.type().pdgcode();
    const PdgCode& pdg_b = data_b.type().pdgcode();
    incoming_parametrized = parametrization_exists(pdg_a, pdg_b);
  }

  // Create ScatterAction object.
  ScatterActionPtr act = std::make_unique<ScatterAction>(
      data_a, data_b, time_until_collision, isotropic_, string_formation_time_,
      box_length_, incoming_parametrized,
      finder_parameters_.spin_interaction_type);

  if (finder_parameters_.coll_crit == CollisionCriterion::Stochastic) {
    act->set_stochastic_pos_idx();
  }

  if (finder_parameters_.strings_switch) {
    act->set_string_interface(string_process_interface_.get());
  }

  // Distance squared calculation not needed for stochastic criterion
  const double distance_squared =
      (finder_parameters_.coll_crit == CollisionCriterion::Geometric)
          ? act->transverse_distance_sqr()
      : (finder_parameters_.coll_crit == CollisionCriterion::Covariant)
          ? act->cov_transverse_distance_sqr()
          : 0.0;

  // Don't calculate cross section if the particles are very far apart.
  // Not needed for stochastic criterion because of cell structure.
  if (finder_parameters_.coll_crit != CollisionCriterion::Stochastic &&
      distance_squared >=
          max_transverse_distance_sqr(finder_parameters_.testparticles)) {
    return nullptr;
  }

  if (incoming_parametrized) {
    act->set_parametrized_total_cross_section(finder_parameters_);
  } else {
    // Add various subprocesses.
    act->add_all_scatterings(finder_parameters_);
  }

  double xs = act->cross_section() * fm2_mb /
              static_cast<double>(finder_parameters_.testparticles);

  // Take cross section scaling factors into account
  xs *= data_a.xsec_scaling_factor(time_until_collision);
  xs *= data_b.xsec_scaling_factor(time_until_collision);

  if (finder_parameters_.coll_crit == CollisionCriterion::Stochastic) {
    const double v_rel = act->relative_velocity();
    /* Collision probability for 2-particle scattering, see
     * \iref{Staudenmaier:2021lrg}. */
    const double prob = xs * v_rel * dt / gcell_vol;

    logg[LFindScatter].debug(
        "Stochastic collison criterion parameters (2-particles):\nprob = ",
        prob, ", xs = ", xs, ", v_rel = ", v_rel, ", dt = ", dt,
        ", gcell_vol = ", gcell_vol,
        ", testparticles = ", finder_parameters_.testparticles);

    if (prob > 1.) {
      std::stringstream err;
      err << "Probability larger than 1 for stochastic rates. ( P_22 = " << prob
          << " )\n"
          << data_a.type().name() << data_b.type().name() << " with masses "
          << data_a.momentum().abs() << " and " << data_b.momentum().abs()
          << " at sqrts[GeV] = " << act->sqrt_s()
          << " with xs[fm^2]/Ntest = " << xs
          << "\nConsider using smaller timesteps.";
      if (finder_parameters_.only_warn_for_high_prob) {
        logg[LFindScatter].warn(err.str());
      } else {
        throw std::runtime_error(err.str());
      }
    }

    // probability criterion
    double random_no = random::uniform(0., 1.);
    if (random_no > prob) {
      return nullptr;
    }

  } else if (finder_parameters_.coll_crit == CollisionCriterion::Geometric ||
             finder_parameters_.coll_crit == CollisionCriterion::Covariant) {
    // just collided with this particle
    if (data_a.id_process() > 0 && data_a.id_process() == data_b.id_process()) {
      logg[LFindScatter].debug("Skipping collided particles at time ",
                               data_a.position().x0(), " due to process ",
                               data_a.id_process(), "\n    ", data_a, "\n<-> ",
                               data_b);

      return nullptr;
    }

    // Cross section for collision criterion
    const double cross_section_criterion = xs * M_1_PI;

    // distance criterion according to cross_section
    if (distance_squared >= cross_section_criterion) {
      return nullptr;
    }

    logg[LFindScatter].debug("particle distance squared: ", distance_squared,
                             "\n    ", data_a, "\n<-> ", data_b);
  }

  // Include possible outgoing branches
  if (incoming_parametrized) {
    act->add_all_scatterings(finder_parameters_);
  }

  return act;
}

ActionPtr ScatterActionsFinder::check_collision_multi_part(
    const ParticleList& plist, double dt, const double gcell_vol) const {
  /* If all particles
   * 1) belong to the two colliding nuclei
   * 2) are within the same nucleus
   * 3) have never experienced any collisons,
   * then the collision between them are banned also for multi-particle
   * interactions. */
  if (!finder_parameters_.allow_collisions_within_nucleus) {
    bool all_projectile =
        std::all_of(plist.begin(), plist.end(), [&](const ParticleData& data) {
          return data.belongs_to() == BelongsTo::Projectile;
        });
    bool all_target =
        std::all_of(plist.begin(), plist.end(), [&](const ParticleData& data) {
          return data.belongs_to() == BelongsTo::Target;
        });
    bool none_collided =
        std::all_of(plist.begin(), plist.end(), [&](const ParticleData& data) {
          return data.get_history().collisions_per_particle == 0;
        });
    if ((all_projectile || all_target) && none_collided) {
      return nullptr;
    }
  }
  // No grid or search in cell
  if (gcell_vol < really_small) {
    return nullptr;
  }

  /* Optimisation for later: Already check here at the beginning
   * if collision with plist is possible before constructing actions. */

  // 1. Determine time of collision.
  const double time_until_collision = dt * random::uniform(0., 1.);

  // 2. Create ScatterAction object.
  ScatterActionMultiPtr act =
      std::make_unique<ScatterActionMulti>(plist, time_until_collision);

  act->set_stochastic_pos_idx();

  // 3. Add possible final states (dt and gcell_vol for probability calculation)
  act->add_possible_reactions(dt, gcell_vol, finder_parameters_.included_multi);

  /* 4. Return total collision probability
   *    Scales with 1 over the number of testpartciles to the power of the
   *    number of incoming particles - 1 */
  const double prob =
      act->get_total_weight() /
      std::pow(finder_parameters_.testparticles, plist.size() - 1);

  // 5. Check that probability is smaller than one
  if (prob > 1.) {
    std::stringstream err;
    err << "Probability " << prob << " larger than 1 for stochastic rates for ";
    for (const ParticleData& data : plist) {
      err << data.type().name();
    }
    err << " at sqrts[GeV] = " << act->sqrt_s()
        << "\nConsider using smaller timesteps.";
    if (finder_parameters_.only_warn_for_high_prob) {
      logg[LFindScatter].warn(err.str());
    } else {
      throw std::runtime_error(err.str());
    }
  }

  // 6. Perform probability decisions
  double random_no = random::uniform(0., 1.);
  if (random_no > prob) {
    return nullptr;
  }

  return act;
}

ActionList ScatterActionsFinder::find_actions_in_cell(
    const ParticleList& search_list, double dt, const double gcell_vol,
    const std::vector<FourVector>& beam_momentum) const {
  std::vector<ActionPtr> actions;
  for (const ParticleData& p1 : search_list) {
    for (const ParticleData& p2 : search_list) {
      // Check for 2 particle scattering
      if (p1.id() < p2.id()) {
        ActionPtr act =
            check_collision_two_part(p1, p2, dt, beam_momentum, gcell_vol);
        if (act) {
          actions.push_back(std::move(act));
        }
      }
      if (finder_parameters_.included_multi.any()) {
        // Also, check for 3 particle scatterings with stochastic criterion
        for (const ParticleData& p3 : search_list) {
          if (finder_parameters_.included_multi
                      [IncludedMultiParticleReactions::Deuteron_3to2] == 1 ||
              finder_parameters_.included_multi
                      [IncludedMultiParticleReactions::Meson_3to1] == 1) {
            if (p1.id() < p2.id() && p2.id() < p3.id()) {
              ActionPtr act =
                  check_collision_multi_part({p1, p2, p3}, dt, gcell_vol);
              if (act) {
                actions.push_back(std::move(act));
              }
            }
          }
          for (const ParticleData& p4 : search_list) {
            if (finder_parameters_.included_multi
                    [IncludedMultiParticleReactions::A3_Nuclei_4to2]) {
              if (p1.id() < p2.id() && p2.id() < p3.id() && p3.id() < p4.id()) {
                ActionPtr act =
                    check_collision_multi_part({p1, p2, p3, p4}, dt, gcell_vol);
                if (act) {
                  actions.push_back(std::move(act));
                }
              }
            }
            if (finder_parameters_.included_multi
                        [IncludedMultiParticleReactions::NNbar_5to2] == 1 &&
                search_list.size() >= 5) {
              for (const ParticleData& p5 : search_list) {
                if ((p1.id() < p2.id() && p2.id() < p3.id() &&
                     p3.id() < p4.id() && p4.id() < p5.id()) &&
                    (p1.is_pion() && p2.is_pion() && p3.is_pion() &&
                     p4.is_pion() && p5.is_pion())) {
                  // at the moment only pure pion 5-body reactions
                  ActionPtr act = check_collision_multi_part(
                      {p1, p2, p3, p4, p5}, dt, gcell_vol);
                  if (act) {
                    actions.push_back(std::move(act));
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return actions;
}

ActionList ScatterActionsFinder::find_actions_with_neighbors(
    const ParticleList& search_list, const ParticleList& neighbors_list,
    double dt, const std::vector<FourVector>& beam_momentum) const {
  std::vector<ActionPtr> actions;
  if (finder_parameters_.coll_crit == CollisionCriterion::Stochastic) {
    // Only search in cells
    return actions;
  }
  for (const ParticleData& p1 : search_list) {
    for (const ParticleData& p2 : neighbors_list) {
      assert(p1.id() != p2.id());
      // Check if a collision is possible.
      ActionPtr act = check_collision_two_part(p1, p2, dt, beam_momentum);
      if (act) {
        actions.push_back(std::move(act));
      }
    }
  }
  return actions;
}

ActionList ScatterActionsFinder::find_actions_with_surrounding_particles(
    const ParticleList& search_list, const Particles& surrounding_list,
    double dt, const std::vector<FourVector>& beam_momentum) const {
  std::vector<ActionPtr> actions;
  if (finder_parameters_.coll_crit == CollisionCriterion::Stochastic) {
    // Only search in cells
    return actions;
  }
  for (const ParticleData& p2 : surrounding_list) {
    /* don't look for collisions if the particle from the surrounding list is
     * also in the search list */
    auto result = std::find_if(
        search_list.begin(), search_list.end(),
        [&p2](const ParticleData& p) { return p.id() == p2.id(); });
    if (result != search_list.end()) {
      continue;
    }
    for (const ParticleData& p1 : search_list) {
      // Check if a collision is possible.
      ActionPtr act = check_collision_two_part(p1, p2, dt, beam_momentum);
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
  std::vector<double> momentum_scan_list = {0.1, 0.3, 0.5, 1.0,
                                            2.0, 3.0, 5.0, 10.0};
  for (const IsoParticleType& A_isotype : IsoParticleType::list_all()) {
    for (const IsoParticleType& B_isotype : IsoParticleType::list_all()) {
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
            ScatterActionPtr act = std::make_unique<ScatterAction>(
                A, B, time, isotropic_, string_formation_time_, -1, false,
                finder_parameters_.spin_interaction_type);
            if (finder_parameters_.strings_switch) {
              act->set_string_interface(string_process_interface_.get());
            }
            act->add_all_scatterings(finder_parameters_);
            const double total_cs = act->cross_section();
            if (total_cs <= 0.0) {
              continue;
            }
            any_nonzero_cs = true;
            for (const auto& channel : act->collision_channels()) {
              const auto type = channel->get_type();
              std::string r;
              if (is_string_soft_process(type) ||
                  type == ProcessType::StringHard) {
                r = A_type->name() + B_type->name() + std::string(" → strings");
              } else {
                std::string r_type =
                    (type == ProcessType::Elastic) ? std::string(" (el)")
                    : (channel->get_type() == ProcessType::TwoToTwo)
                        ? std::string(" (inel)")
                        : std::string(" (?)");
                r = A_type->name() + B_type->name() + std::string(" → ") +
                    channel->particle_types()[0]->name() +
                    channel->particle_types()[1]->name() + r_type;
              }
              isoclean(r);
              r_list.push_back(r);
            }
          }
        }
      }
      std::sort(r_list.begin(), r_list.end());
      r_list.erase(std::unique(r_list.begin(), r_list.end()), r_list.end());
      if (any_nonzero_cs) {
        for (auto r : r_list) {
          std::cout << r;
          if (r_list.back() != r) {
            std::cout << ", ";
          }
        }
        std::cout << std::endl;
      }
    }
  }
}

/// Represent a final-state cross section.
struct FinalStateCrossSection {
  /// Name of the final state.
  std::string name_;

  /// Corresponding cross section in mb.
  double cross_section_;

  /// Total mass of final state particles.
  double mass_;

  /**
   * Construct a final-state cross section.
   *
   * \param name Name of the final state.
   * \param cross_section Corresponding cross section in mb.
   * \param mass Total mass of final state particles.
   * \return Constructed object.
   */
  FinalStateCrossSection(const std::string& name, double cross_section,
                         double mass)
      : name_(name), cross_section_(cross_section), mass_(mass) {}
};

namespace decaytree {

/**
 * Node of a decay tree, representing a possible action (2-to-2 or 1-to-2).
 *
 * This data structure can be used to build a tree going from the initial
 * state (a collision of two particles) to all possible final states by
 * recursively performing all possible decays. The tree can be used to
 * calculate the final state cross sections.
 *
 * The initial actions are 2-to-2 or 2-to-1 scatterings, all other actions are
 * 1-to-2 decays.
 */
struct Node {
 public:
  /// Name for printing.
  std::string name_;

  /// Weight (cross section or branching ratio).
  double weight_;

  /// Initial-state particle types in this action.
  ParticleTypePtrList initial_particles_;

  /// Final-state particle types in this action.
  ParticleTypePtrList final_particles_;

  /// Particle types corresponding to the global state after this action.
  ParticleTypePtrList state_;

  /// Possible actions after this action.
  std::vector<Node> children_;

  /// Cannot be copied
  Node(const Node&) = delete;
  /// Move constructor
  Node(Node&&) = default;

  /**
   * \return A new decay tree node.
   *
   * \param name Name for printing.
   * \param weight Cross section or branching ratio.
   * \param initial_particles Initial-state particle types in this node.
   * \param final_particles Final-state particle types in this node.
   * \param state Curent particle types of the system.
   * \param children Possible actions after this action.
   */
  Node(const std::string& name, double weight,
       ParticleTypePtrList&& initial_particles,
       ParticleTypePtrList&& final_particles, ParticleTypePtrList&& state,
       std::vector<Node>&& children)
      : name_(name),
        weight_(weight),
        initial_particles_(std::move(initial_particles)),
        final_particles_(std::move(final_particles)),
        state_(std::move(state)),
        children_(std::move(children)) {}

  /**
   * Add an action to the children of this node.
   *
   * The current particle state of the new action is automatically calculated.
   *
   * \param name Name of the action used for output.
   * \param weight Cross section/branching ratio of the action.
   * \param initial_particles Initial-state particle types of the action.
   * \param final_particles Final-state particle types of the action.
   * \return Newly added node by reference.
   */
  Node& add_action(const std::string& name, double weight,
                   ParticleTypePtrList&& initial_particles,
                   ParticleTypePtrList&& final_particles) {
    // Copy parent state and update it.
    ParticleTypePtrList state(state_);
    for (const auto& p : initial_particles) {
      state.erase(std::find(state.begin(), state.end(), p));
    }
    for (const auto& p : final_particles) {
      state.push_back(p);
    }
    // Sort the state to normalize the output.
    std::sort(state.begin(), state.end(),
              [](ParticleTypePtr a, ParticleTypePtr b) {
                return a->name() < b->name();
              });
    // Push new node to children.
    Node new_node(name, weight, std::move(initial_particles),
                  std::move(final_particles), std::move(state), {});
    children_.emplace_back(std::move(new_node));
    return children_.back();
  }

  /// Print the decay tree starting with this node.
  void print() const { print_helper(0); }

  /**
   * \return Final-state cross sections.
   */
  std::vector<FinalStateCrossSection> final_state_cross_sections() const {
    std::vector<FinalStateCrossSection> result;
    final_state_cross_sections_helper(0, result, "", 1.);
    return result;
  }

 private:
  /**
   * Internal helper function for `print`, to be called recursively to print
   * all nodes.
   *
   * \param depth Recursive call depth.
   */
  void print_helper(uint64_t depth) const {
    for (uint64_t i = 0; i < depth; i++) {
      std::cout << " ";
    }
    std::cout << name_ << " " << weight_ << std::endl;
    for (const auto& child : children_) {
      child.print_helper(depth + 1);
    }
  }

  /**
   * Internal helper function for `final_state_cross_sections`, to be called
   * recursively to calculate all final-state cross sections.
   *
   * \param depth Recursive call depth.
   * \param result Pairs of process names and exclusive cross sections.
   * \param name Current name.
   * \param weight current Weight/cross section.
   * \param show_intermediate_states Whether intermediate states should be
   * shown.
   */
  void final_state_cross_sections_helper(
      uint64_t depth, std::vector<FinalStateCrossSection>& result,
      const std::string& name, double weight,
      bool show_intermediate_states = false) const {
    // The first node corresponds to the total cross section and has to be
    // ignored. The second node corresponds to the partial cross section. All
    // further nodes correspond to branching ratios.
    if (depth > 0) {
      weight *= weight_;
    }

    std::string new_name;
    double mass = 0.;

    if (show_intermediate_states) {
      new_name = name;
      if (!new_name.empty()) {
        new_name += "->";
      }
      new_name += name_;
      new_name += "{";
    } else {
      new_name = "";
    }
    for (const auto& s : state_) {
      new_name += s->name();
      mass += s->mass();
    }
    if (show_intermediate_states) {
      new_name += "}";
    }

    if (children_.empty()) {
      result.emplace_back(FinalStateCrossSection(new_name, weight, mass));
      return;
    }
    for (const auto& child : children_) {
      child.final_state_cross_sections_helper(depth + 1, result, new_name,
                                              weight, show_intermediate_states);
    }
  }
};

/**
 * Generate name for decay and update final state.
 *
 * \param[in] res_name Name of resonance.
 * \param[in] decay Decay branch.
 * \param[out] final_state Final state of decay.
 * \return Name of decay.
 */
static std::string make_decay_name(const std::string& res_name,
                                   const DecayBranchPtr& decay,
                                   ParticleTypePtrList& final_state) {
  std::stringstream name;
  name << "[" << res_name << "->";
  for (const auto& p : decay->particle_types()) {
    name << p->name();
    final_state.push_back(p);
  }
  name << "]";
  return name.str();
}

/**
 * Add nodes for all decays possible from the given node and all of its
 * children.
 *
 * \param node Starting node.
 * \param[in] sqrts center-of-mass energy.
 */
static void add_decays(Node& node, double sqrts) {
  // If there is more than one unstable particle in the current state, then
  // there will be redundant paths in the decay tree, corresponding to
  // reorderings of the decays. To avoid double counting, we normalize by the
  // number of possible decay orderings. Normalizing by the number of unstable
  // particles recursively corresponds to normalizing by the factorial that
  // gives the number of reorderings.
  //
  // Ideally, the redundant paths should never be added to the decay tree, but
  // we never have more than two redundant paths, so it probably does not
  // matter much.
  uint32_t n_unstable = 0;
  double sqrts_minus_masses = sqrts;
  for (const ParticleTypePtr ptype : node.state_) {
    if (!ptype->is_stable()) {
      n_unstable += 1;
    }
    sqrts_minus_masses -= ptype->mass();
  }
  const double norm =
      n_unstable != 0 ? 1. / static_cast<double>(n_unstable) : 1.;

  for (const ParticleTypePtr ptype : node.state_) {
    if (!ptype->is_stable()) {
      const double sqrts_decay = sqrts_minus_masses + ptype->mass();
      bool can_decay = false;
      for (const auto& decay : ptype->decay_modes().decay_mode_list()) {
        // Make sure to skip kinematically impossible decays.
        // In principle, we would have to integrate over the mass of the
        // resonance, but as an approximation we just assume it at its pole.
        double final_state_mass = 0.;
        for (const auto& p : decay->particle_types()) {
          final_state_mass += p->mass();
        }
        if (final_state_mass > sqrts_decay) {
          continue;
        }
        can_decay = true;

        ParticleTypePtrList parts;
        const auto name = make_decay_name(ptype->name(), decay, parts);
        auto& new_node = node.add_action(name, norm * decay->weight(), {ptype},
                                         std::move(parts));
        add_decays(new_node, sqrts_decay);
      }
      if (!can_decay) {
        // Remove final-state cross sections with resonances that cannot
        // decay due to our "mass = pole mass" approximation.
        node.weight_ = 0;
        return;
      }
    }
  }
}

}  // namespace decaytree

/**
 * Deduplicate the final-state cross sections by summing.
 *
 * \param[inout] final_state_xs Final-state cross sections.
 */
static void deduplicate(std::vector<FinalStateCrossSection>& final_state_xs) {
  std::sort(final_state_xs.begin(), final_state_xs.end(),
            [](const FinalStateCrossSection& a,
               const FinalStateCrossSection& b) { return a.name_ < b.name_; });
  auto current = final_state_xs.begin();
  while (current != final_state_xs.end()) {
    auto adjacent = std::adjacent_find(
        current, final_state_xs.end(),
        [](const FinalStateCrossSection& a, const FinalStateCrossSection& b) {
          return a.name_ == b.name_;
        });
    current = adjacent;
    if (adjacent != final_state_xs.end()) {
      adjacent->cross_section_ += (adjacent + 1)->cross_section_;
      final_state_xs.erase(adjacent + 1);
    }
  }
}

void ScatterActionsFinder::dump_cross_sections(
    const ParticleType& a, const ParticleType& b, double m_a, double m_b,
    bool final_state, std::vector<double>& plab) const {
  typedef std::vector<std::pair<double, double>> xs_saver;
  std::map<std::string, xs_saver> xs_dump;
  std::map<std::string, double> outgoing_total_mass;

  ParticleData a_data(a), b_data(b);
  int n_momentum_points = 200;
  constexpr double momentum_step = 0.02;
  if (plab.size() > 0) {
    n_momentum_points = plab.size();
    // Remove duplicates.
    std::sort(plab.begin(), plab.end());
    plab.erase(std::unique(plab.begin(), plab.end()), plab.end());
  }
  for (int i = 0; i < n_momentum_points; i++) {
    double momentum;
    if (plab.size() > 0) {
      momentum = pCM_from_s(s_from_plab(plab.at(i), m_a, m_b), m_a, m_b);
    } else {
      momentum = momentum_step * (i + 1);
    }
    a_data.set_4momentum(m_a, momentum, 0.0, 0.0);
    b_data.set_4momentum(m_b, -momentum, 0.0, 0.0);
    const double sqrts = (a_data.momentum() + b_data.momentum()).abs();
    const ParticleList incoming = {a_data, b_data};
    ScatterActionPtr act = std::make_unique<ScatterAction>(
        a_data, b_data, 0., isotropic_, string_formation_time_, -1, false,
        finder_parameters_.spin_interaction_type);
    if (finder_parameters_.strings_switch) {
      act->set_string_interface(string_process_interface_.get());
    }
    act->add_all_scatterings(finder_parameters_);
    decaytree::Node tree(a.name() + b.name(), act->cross_section(), {&a, &b},
                         {&a, &b}, {&a, &b}, {});
    const CollisionBranchList& processes = act->collision_channels();
    for (const auto& process : processes) {
      const double xs = process->weight();
      if (xs <= 0.0) {
        continue;
      }
      if (!final_state) {
        std::stringstream process_description_stream;
        process_description_stream << *process;
        const std::string& description = process_description_stream.str();
        double m_tot = 0.0;
        for (const auto& ptype : process->particle_types()) {
          m_tot += ptype->mass();
        }
        outgoing_total_mass[description] = m_tot;
        if (!xs_dump[description].empty() &&
            std::abs(xs_dump[description].back().first - sqrts) <
                really_small) {
          xs_dump[description].back().second += xs;
        } else {
          xs_dump[description].push_back(std::make_pair(sqrts, xs));
        }
      } else {
        std::stringstream process_description_stream;
        process_description_stream << *process;
        const std::string& description = process_description_stream.str();
        ParticleTypePtrList initial_particles = {&a, &b};
        ParticleTypePtrList final_particles = process->particle_types();
        auto& process_node =
            tree.add_action(description, xs, std::move(initial_particles),
                            std::move(final_particles));
        decaytree::add_decays(process_node, sqrts);
      }
    }
    xs_dump["total"].push_back(std::make_pair(sqrts, act->cross_section()));
    // Total cross-section should be the first in the list -> negative mass
    outgoing_total_mass["total"] = -1.0;
    if (final_state) {
      // tree.print();
      auto final_state_xs = tree.final_state_cross_sections();
      deduplicate(final_state_xs);
      for (const auto& p : final_state_xs) {
        // Don't print empty columns.
        //
        // FIXME(steinberg): The better fix would be to not have them in the
        // first place.
        if (p.name_ == "") {
          continue;
        }
        outgoing_total_mass[p.name_] = p.mass_;
        xs_dump[p.name_].push_back(std::make_pair(sqrts, p.cross_section_));
      }
    }
  }
  // Get rid of cross sections that are zero.
  // (This only happens if their is a resonance in the final state that cannot
  // decay with our simplified assumptions.)
  for (auto it = begin(xs_dump); it != end(xs_dump);) {
    // Sum cross section over all energies.
    const xs_saver& xs = (*it).second;
    double sum = 0;
    for (const auto& p : xs) {
      sum += p.second;
    }
    if (sum == 0.) {
      it = xs_dump.erase(it);
    } else {
      ++it;
    }
  }

  // Nice ordering of channels by summed pole mass of products
  std::vector<std::string> all_channels;
  for (const auto& channel : xs_dump) {
    all_channels.push_back(channel.first);
  }
  std::sort(all_channels.begin(), all_channels.end(),
            [&](const std::string& str_a, const std::string& str_b) {
              return outgoing_total_mass[str_a] < outgoing_total_mass[str_b];
            });

  // Print header
  std::cout << "# Dumping partial " << a.name() << b.name()
            << " cross-sections in mb, energies in GeV" << std::endl;
  std::cout << "   sqrt_s";
  // Align everything to 24 unicode characters.
  // This should be enough for the longest channel name: 7 final-state
  // particles, or 2 of the longest named resonances (currently Ds0*(2317)⁺).
  for (const auto& channel : all_channels) {
    std::cout << utf8::fill_left(channel, 24, ' ');
  }
  std::cout << std::endl;

  // Print out all partial cross-sections in mb
  for (int i = 0; i < n_momentum_points; i++) {
    double momentum;
    if (plab.size() > 0) {
      momentum = pCM_from_s(s_from_plab(plab.at(i), m_a, m_b), m_a, m_b);
    } else {
      momentum = momentum_step * (i + 1);
    }
    a_data.set_4momentum(m_a, momentum, 0.0, 0.0);
    b_data.set_4momentum(m_b, -momentum, 0.0, 0.0);
    const double sqrts = (a_data.momentum() + b_data.momentum()).abs();
    std::printf("%9.6f", sqrts);
    for (const auto& channel : all_channels) {
      const xs_saver energy_and_xs = xs_dump[channel];
      size_t j = 0;
      for (; j < energy_and_xs.size() && energy_and_xs[j].first < sqrts; j++) {
      }
      double xs = 0.0;
      if (j < energy_and_xs.size() &&
          std::abs(energy_and_xs[j].first - sqrts) < really_small) {
        xs = energy_and_xs[j].second;
      }
      std::printf("%16.6f", xs);  // Same alignment as in the header.
    }
    std::printf("\n");
  }
}

}  // namespace smash
