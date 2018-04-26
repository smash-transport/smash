/*
 *
 *    Copyright (c) 2014-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteractionsfinder.h"

#include <algorithm>
#include <map>
#include <vector>

#include "include/configuration.h"
#include "include/constants.h"
#include "include/crosssections.h"
#include "include/cxx14compat.h"
#include "include/decaymodes.h"
#include "include/experimentparameters.h"
#include "include/isoparticletype.h"
#include "include/logging.h"
#include "include/macros.h"
#include "include/particles.h"
#include "include/scatteraction.h"
#include "include/scatteractionphoton.h"
#include "include/stringfunctions.h"

namespace smash {
/*!\Userguide
 * \page input_collision_term_ Collision_Term
 * \key Elastic_Cross_Section (double, optional, default = -1.0 [mb]) \n
 * If a non-negative value is given, it will override the parametrized
 * elastic cross sections (which are energy-dependent) with a constant value.
 * This constant elastic cross section is used for all collisions.
 *
 * \key Isotropic (bool, optional, default = \key false) \n
 * Do all collisions isotropically.
 *
 * \key Elastic_NN_Cutoff_Sqrts (double, optional, default = 1.98): \n
 * The elastic collisions betwen two nucleons with sqrt_s below
 * Elastic_NN_Cutoff_Sqrts, in GeV, cannot happen. \n
 * \li \key Elastic_NN_Cutoff_Sqrts < 1.88 - Below the threshold energy of the
 * elastic collsion, no effect \n
 * \li \key Elastic_NN_Cutoff_Sqrts > 2.02 - Beyond the threshold energy of the
 * inelastic collision NN->NNpi, not suggested
 *
 * \key Strings (bool, optional, default = \key true for each setup except box): \n
 * \li \key true - String excitation is enabled\n
 * \li \key false - String excitation is disabled
 *
 * \key String_Formation_Time (double, optional, default = 1.0): \n
 * Parameter for formation time in string fragmentation, in fm/c.
 *
 * \subpage string_parameters
 *
 * \page string_parameters String_Parameters
 * A set of parameters with which the string fragmentation can be modified.
 * (TODO(Mohs): Explain them in one sentence and add formulas.) \n
 *
 * \key String_Tension (double, optional, default = 1.0 GeV/fm) \n
 *
 * \key Gluon_Beta (double, optional, default = 0.5) \n
 *
 * \key Gluon_Pmin (double, optional, default = 0.001 GeV) \n
 *
 * \key Quark_Alpha (double, optional, default = 1.0) \n
 *
 * \key Quark_Beta (double, optional, default = 2.5) \n
 *
 * \key Strange_Supp (double, optional, default = 0.217 as in Pythia) \n
 * Strangeness suppression factor.
 *
 * \key Diquark_Supp (double, optional, default = 0.081 as in Pythia) \n
 * Diquark suppression factor.
 *
 * \key Sigma_Perp (double, optional, default = 0.7 ) \n
 *
 * \key StringZ_A (double, optional, default = 0.68 as in Pythia) \n
 *
 * \key StringZ_B (double, optional, default = 0.98 as in Pythia) \n
 *
 * \key String_Sigma_T (double, optional, default = 0.25)
 *
 */

ScatterActionsFinder::ScatterActionsFinder(
    Configuration config, const ExperimentParameters &parameters,
    const std::vector<bool> &nucleon_has_interacted, int N_tot, int N_proj,
    int n_fractional_photons = 1)
    : elastic_parameter_(
          config.take({"Collision_Term", "Elastic_Cross_Section"}, -1.)),
      testparticles_(parameters.testparticles),
      isotropic_(config.take({"Collision_Term", "Isotropic"}, false)),
      two_to_one_(parameters.two_to_one),
      incl_set_(parameters.included_2to2),
      low_snn_cut_(parameters.low_snn_cut),
      strings_switch_(parameters.strings_switch),
      use_AQM_(parameters.use_AQM),
      use_transition_probability_(parameters.use_transition_probability),
      nnbar_treatment_(parameters.nnbar_treatment),
      nucleon_has_interacted_(nucleon_has_interacted),
      N_tot_(N_tot),
      N_proj_(N_proj),
      string_formation_time_(config.take(
          {"Collision_Term", "String_Parameters", "Formation_Time"}, 1.)),
      photons_(parameters.photons_switch),
      n_fractional_photons_(n_fractional_photons) {
  if (is_constant_elastic_isotropic()) {
    const auto &log = logger<LogArea::FindScatter>();
    log.info("Constant elastic isotropic cross-section mode:", " using ",
             elastic_parameter_, " mb as maximal cross-section.");
  }
  if (strings_switch_) {
    auto subconfig = config["Collision_Term"]["String_Parameters"];
    string_process_interface_ =
        make_unique<StringProcess>(subconfig.take({"String_Tension"}, 1.0),
                                   string_formation_time_,
                                   subconfig.take({"Gluon_Beta"}, 0.5),
                                   subconfig.take({"Gluon_Pmin"}, 0.001),
                                   subconfig.take({"Quark_Alpha"}, 1.0),
                                   subconfig.take({"Quark_Beta"}, 2.5),
                                   subconfig.take({"Strange_Supp"}, 0.217),
                                   subconfig.take({"Diquark_Supp"}, 0.081),
                                   subconfig.take({"Sigma_Perp"}, 0.7),
                                   subconfig.take({"StringZ_A"}, 0.68),
                                   subconfig.take({"StringZ_B"}, 0.98),
                                   subconfig.take({"String_Sigma_T"}, 0.25));
  }
}

ActionPtr ScatterActionsFinder::check_collision(const ParticleData &data_a,
                                                const ParticleData &data_b,
                                                double dt) const {
#ifndef NDEBUG
  const auto &log = logger<LogArea::FindScatter>();
#endif

  // just collided with this particle
  if (data_a.id_process() > 0 && data_a.id_process() == data_b.id_process()) {
#ifndef NDEBUG
    log.debug("Skipping collided particles at time ", data_a.position().x0(),
              " due to process ", data_a.id_process(), "\n    ", data_a,
              "\n<-> ", data_b);
#endif
    return nullptr;
  }
  /* If the two particles
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

  // Determine time of collision.
  const double time_until_collision = collision_time(data_a, data_b);

  // Check that collision happens in this timestep.
  if (time_until_collision < 0. || time_until_collision >= dt) {
    return nullptr;
  }

  // Create ScatterAction object.
  ScatterActionPtr act = make_unique<ScatterAction>(
      data_a, data_b, time_until_collision, isotropic_, string_formation_time_);
  if (strings_switch_) {
    act->set_string_interface(string_process_interface_.get());
  }

  const double distance_squared = act->transverse_distance_sqr();

  // Don't calculate cross section if the particles are very far apart.
  if (distance_squared >= max_transverse_distance_sqr(testparticles_)) {
    return nullptr;
  }

  // Add various subprocesses.
  act->add_all_scatterings(elastic_parameter_, two_to_one_, incl_set_,
                           low_snn_cut_, strings_switch_, use_AQM_,
                           use_transition_probability_,  nnbar_treatment_);

  // Cross section for collision criterion
  double cross_section_criterion = act->cross_section() * fm2_mb * M_1_PI /
                                   static_cast<double>(testparticles_);
  /* Consider cross section scaling factors only if the particles
   * are not formed yet at the prospective time of the interaction */
  if (data_a.formation_time() > data_a.position().x0() + time_until_collision) {
    cross_section_criterion *= data_a.cross_section_scaling_factor();
  }
  if (data_b.formation_time() > data_b.position().x0() + time_until_collision) {
    cross_section_criterion *= data_b.cross_section_scaling_factor();
  }

  // distance criterion according to cross_section
  if (distance_squared >= cross_section_criterion) {
    return nullptr;
  }

#ifndef NDEBUG
  log.debug("particle distance squared: ", distance_squared, "\n    ", data_a,
            "\n<-> ", data_b);
#endif

  return std::move(act);
}

ActionList ScatterActionsFinder::find_actions_in_cell(
    const ParticleList &search_list, double dt) const {
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
    double dt) const {
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
    double dt) const {
  std::vector<ActionPtr> actions;
  for (const ParticleData &p2 : surrounding_list) {
    /* don't look for collisions if the particle from the surrounding list is
     * also in the search list */
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
  std::vector<double> momentum_scan_list = {0.1, 0.3, 0.5, 1.0,
                                            2.0, 3.0, 5.0, 10.0};
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
            ScatterActionPtr act = make_unique<ScatterAction>(
                A, B, time, isotropic_, string_formation_time_);
            act->add_all_scatterings(elastic_parameter_, two_to_one_, incl_set_,
                                     low_snn_cut_, strings_switch_, use_AQM_,
                                     use_transition_probability_,
                                     nnbar_treatment_);
            const double total_cs = act->cross_section();
            if (total_cs <= 0.0) {
              continue;
            }
            any_nonzero_cs = true;
            for (const auto &channel : act->collision_channels()) {
              const auto type = channel->get_type();
              std::string r;
              if (type == ProcessType::StringSoft) {
                r = A_type->name() + B_type->name() +
                    std::string(" → strings (soft)");
              } else if (type == ProcessType::StringHard) {
                r = A_type->name() + B_type->name() +
                    std::string(" → strings (hard)");
              } else {
                std::string r_type =
                    (type == ProcessType::Elastic)
                        ? std::string(" (el)")
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

void ScatterActionsFinder::dump_cross_sections(const ParticleType &a,
                                               const ParticleType &b,
                                               double m_a, double m_b) const {
  typedef std::vector<std::pair<double, double>> xs_saver;
  std::map<std::string, xs_saver> xs_dump;
  std::map<std::string, double> outgoing_total_mass;

  ParticleData a_data(a), b_data(b);
  constexpr int n_momentum_points = 200;
  constexpr double momentum_step = 0.02;
  for (int i = 1; i < n_momentum_points; i++) {
    const double momentum = momentum_step * i;
    a_data.set_4momentum(m_a, momentum, 0.0, 0.0);
    b_data.set_4momentum(m_b, -momentum, 0.0, 0.0);
    const double sqrts = (a_data.momentum() + b_data.momentum()).abs();
    const ParticleList incoming = {a_data, b_data};
    cross_sections xs_class(incoming, sqrts);
    // Create ScatterAction object.
    ScatterActionPtr act = make_unique<ScatterAction>(
        a_data, b_data, 0., isotropic_, string_formation_time_);
    if (strings_switch_) {
      act->set_string_interface(string_process_interface_.get());
    }
    act->add_all_scatterings(elastic_parameter_, two_to_one_, incl_set_,
                             low_snn_cut_, strings_switch_, use_AQM_,
                             use_transition_probability_,
                             nnbar_treatment_);
    const CollisionBranchList& processes = act->collision_channels();
    for (const auto &process : processes) {
      const double xs = process->weight();
      if (xs <= 0.0) {
        continue;
      }
      std::stringstream process_description_stream;
      process_description_stream << *process;
      std::string description = process_description_stream.str();
      double m_tot = 0.0;
      for (const auto ptype : process->particle_types()) {
        m_tot += ptype->mass();
      }
      outgoing_total_mass[description] = m_tot;
      if (!xs_dump[description].empty() &&
          std::abs(xs_dump[description].back().first - sqrts) < really_small) {
        xs_dump[description].back().second += xs;
      } else {
        xs_dump[description].push_back(std::make_pair(sqrts, xs));
      }
    }
  }

  // Nice ordering of channels by summed pole mass of products
  std::vector<std::string> all_channels;
  for (const auto channel : xs_dump) {
    all_channels.push_back(channel.first);
  }
  std::sort(all_channels.begin(), all_channels.end(),
            [&](const std::string &str_a, const std::string &str_b) {
              return outgoing_total_mass[str_a] < outgoing_total_mass[str_b];
            });

  // Print header
  std::cout << "# Dumping partial cross-sections in mb" << std::endl;
  std::cout << "# sqrt(s) [GeV], " << a.name() << b.name() << "→ ";
  for (const auto channel : all_channels) {
    std::cout << utf8::fill_left(channel, 12, ' ');
  }
  std::cout << std::endl;

  // Print out all partial cross-sections in mb
  for (int i = 1; i < n_momentum_points; i++) {
    const double momentum = momentum_step * i;
    a_data.set_4momentum(m_a, momentum, 0.0, 0.0);
    b_data.set_4momentum(m_b, -momentum, 0.0, 0.0);
    const double sqrts = (a_data.momentum() + b_data.momentum()).abs();
    printf("%17.5f      ", sqrts);
    for (const auto channel : all_channels) {
      const xs_saver energy_and_xs = xs_dump[channel];
      size_t j = 0;
      for (; j < energy_and_xs.size() && energy_and_xs[j].first < sqrts; j++) {
      }
      double xs = 0.0;
      if (j < energy_and_xs.size() &&
          std::abs(energy_and_xs[j].first - sqrts) < really_small) {
        xs = energy_and_xs[j].second;
      }
      printf("%12.6f", xs);
    }
    printf("\n");
  }
}

}  // namespace smash
