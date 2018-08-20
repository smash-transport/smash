/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/scatteractionsfinder.h"

#include <algorithm>
#include <map>
#include <vector>

#include "smash/configuration.h"
#include "smash/constants.h"
#include "smash/crosssections.h"
#include "smash/cxx14compat.h"
#include "smash/decaymodes.h"
#include "smash/experimentparameters.h"
#include "smash/isoparticletype.h"
#include "smash/logging.h"
#include "smash/macros.h"
#include "smash/particles.h"
#include "smash/scatteraction.h"
#include "smash/scatteractionphoton.h"
#include "smash/stringfunctions.h"

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
 * \key Strings (bool, optional, default = \key true for each setup except box):
 * \n \li \key true - String excitation is enabled\n \li \key false - String
 * excitation is disabled
 *
 * \key String_Formation_Time (double, optional, default = 1.0): \n
 * Parameter for formation time in string fragmentation, in fm/c.
 *
 * \subpage string_parameters
 *
 * \page string_parameters String_Parameters
 * A set of parameters with which the string fragmentation can be modified.
 *
 * \key String_Tension (double, optional, default = 1.0 GeV/fm) \n
 * String tension \f$\kappa\f$ connecting massless quarks in Hamiltonian:
 * \f[H=|p_1|+|p_2|+\kappa |x_1-x_2|\f]
 * This parameter is only used to determine particles' formation times
 * according to the yo-yo formalism (in the soft string routine for now).
 *
 * \key Gluon_Beta (double, optional, default = 0.5) \n
 * Parameter \f$\beta\f$ in parton distribution function for gluons:
 * \f[\mathrm{PDF}_g(x) \propto \frac{1}{x}(1-x)^{\beta+1}\f]
 *
 * \key Gluon_Pmin (double, optional, default = 0.001 GeV) \n
 * Smallest possible scale for gluon lightcone momentum.
 * This is divided by sqrts to get the minimum fraction to be sampled
 * from PDF shown above.
 *
 * \key Quark_Alpha (double, optional, default = 2.0) \n
 * Parameter \f$\alpha\f$ in parton distribution function for quarks:
 * \f[\mathrm{PDF}_q\propto x^{\alpha-1}(1-x)^{\beta-1}\f]
 *
 * \key Quark_Beta (double, optional, default = 5.0) \n
 * Parameter \f$\beta\f$ in PDF for quarks shown above.
 *
 * \key Strange_Supp (double, optional, default = 0.165) \n
 * Strangeness suppression factor \f$\lambda\f$:
 * \f[\lambda=\frac{P(s\bar{s})}{P(u\bar{u})}=\frac{P(s\bar{s})}{P(d\bar{d})}\f]
 * Defines the probability to produce a \f$s\bar{s}\f$ pair relative to produce
 * a light \f$q\bar{q}\f$ pair
 *
 * \key Diquark_Supp (double, optional, default = 0.042) \n
 * Diquark suppression factor. Defines the probability to produce a diquark
 * antidiquark pair relative to producing a qurk antiquark pair.
 *
 * \key Sigma_Perp (double, optional, default = 0.4) \n
 * Parameter \f$\sigma_\perp\f$ in distribution for transverse momentum
 * transfer between colliding hadrons \f$p_\perp\f$ and string mass \f$M_X\f$:
 * \f[\frac{d^3N}{dM^2_Xd^2\mathbf{p_\perp}}\propto \frac{1}{M_X^2}
 * \exp\left(-\frac{p_\perp^2}{\sigma_\perp^2}\right)\f]
 *
 * \key StringZ_A (double, optional, default = 0.68 as in Pythia) \n
 * Parameter a in pythia fragmentation function \f$f(z)\f$:
 * \f[f(z) = \frac{1}{z} (1-z)^a \exp\left(-b\frac{m_T^2}{z}\right)\f]
 *
 * \key StringZ_B (double, optional, default = 0.5) \n
 * Parameter \f$b\f$ in pythia fragmentation function shown above.
 *
 * \key String_Sigma_T (double, optional, default = 0.5)
 * Standard deviation in Gaussian for transverse momentum distributed to
 * string fragments during fragmentation.
 *
 * \key Form_Time_Factor (double, optional, default = 1.0) \n
 * Factor to be multiplied with the formation time of string fragments from
 * the soft string routine.
 * \key Power_Particle_Formation (double, optional, default = -1.0) \n
 * If positive, the power with which the cross section scaling factor of
 * string fragments grows in time until it reaches 1. If negative, the scaling
 * factor will be constant and jump to 1 once the particle forms.
 *
 * \key Use_Yoyo_Model (bool, optional, default = true)
 * Calculate the string fragments' formation times from the yoyo-model if
 * enabled. Else constant proper time defined by Formation_Time is used.
 *
 * \page input_collision_term_ Collision_Term
 * \n
 * Example: Configuring the Collision Term
 * --------------
 * The following example configures SMASH to include all but
 * strangeness exchange involving 2 <--> 2 scatterings, to treat N + Nbar
 * processes as resonance formations and to not force decays at the end of the
 * simulation. The elastic cross section is globally set to 30 mbarn and the
 * \f$ \sqrt{s} \f$ cutoff for elastic nucleon + nucleon collisions is 1.93 GeV.
 * All collisions are performed isotropically and 2 <--> 1 processes are
 * forbidden.
 *
 *\verbatim
 Collision_Term:
     Included_2to2:    ["Elastic", "NN_to_NR", "NN_to_DR", "KN_to_KN",
 "KN_to_KDelta"] Two_to_One: True Force_Decays_At_End: False NNbar_Treatment:
 "resonances" Elastic_Cross_Section: 30.0 Elastic_NN_Cutoff_Sqrts: 1.93
     Isotropic: True
 \endverbatim
 *
 * If necessary, all collisions can be turned off by inserting
 *\verbatim
     No_Collisions: True
 \endverbatim
 * in the configuration file. \n
 * \n
 * Additionally, string fragmentation can be activated. If desired, the user can
 * also configure the string parameters.
 *
 *\verbatim
     Strings: True
     String_Parameters:
         String_Tension: 1.0
         Gluon_Beta: 0.5
         Gluon_Pmin: 0.001
         Quark_Alpha: 1.0
         Quark_Beta: 2.5
         Strange_Supp: 0.217
         Diquark_Supp: 0.081
         Sigma_Perp: 0.7
         StringZ_A: 0.68
         StringZ_B: 0.98
         String_Sigma_T: 0.25
 \endverbatim
 *
 * Pauli Blocking can further be activated by means of the following subsection
 *\verbatim
     Pauli_Blocking:
         Spatial_Averaging_Radius: 1.86
         Momentum_Averaging_Radius: 0.08
         Gaussian_Cutoff: 2.2
 \endverbatim
 */

ScatterActionsFinder::ScatterActionsFinder(
    Configuration config, const ExperimentParameters& parameters,
    const std::vector<bool>& nucleon_has_interacted, int N_tot, int N_proj)
    : elastic_parameter_(
          config.take({"Collision_Term", "Elastic_Cross_Section"}, -1.)),
      testparticles_(parameters.testparticles),
      isotropic_(config.take({"Collision_Term", "Isotropic"}, false)),
      two_to_one_(parameters.two_to_one),
      incl_set_(parameters.included_2to2),
      low_snn_cut_(parameters.low_snn_cut),
      strings_switch_(parameters.strings_switch),
      use_AQM_(parameters.use_AQM),
      strings_with_probability_(parameters.strings_with_probability),
      nnbar_treatment_(parameters.nnbar_treatment),
      nucleon_has_interacted_(nucleon_has_interacted),
      N_tot_(N_tot),
      N_proj_(N_proj),
      string_formation_time_(config.take(
          {"Collision_Term", "String_Parameters", "Formation_Time"}, 1.)) {
  if (is_constant_elastic_isotropic()) {
    const auto& log = logger<LogArea::FindScatter>();
    log.info("Constant elastic isotropic cross-section mode:", " using ",
             elastic_parameter_, " mb as maximal cross-section.");
  }
  if (strings_switch_) {
    auto subconfig = config["Collision_Term"]["String_Parameters"];
    string_process_interface_ = make_unique<StringProcess>(
        subconfig.take({"String_Tension"}, 1.0), string_formation_time_,
        subconfig.take({"Gluon_Beta"}, 0.5),
        subconfig.take({"Gluon_Pmin"}, 0.001),
        subconfig.take({"Quark_Alpha"}, 2.0),
        subconfig.take({"Quark_Beta"}, 5.0),
        subconfig.take({"Strange_Supp"}, 0.165),
        subconfig.take({"Diquark_Supp"}, 0.042),
        subconfig.take({"Sigma_Perp"}, 0.4),
        subconfig.take({"StringZ_A"}, 0.68), subconfig.take({"StringZ_B"}, 0.5),
        subconfig.take({"String_Sigma_T"}, 0.5),
        subconfig.take({"Form_Time_Factor"}, 1.0),
        subconfig.take({"Use_Yoyo_Model"}, true));
  }
}

ActionPtr ScatterActionsFinder::check_collision(const ParticleData& data_a,
                                                const ParticleData& data_b,
                                                double dt) const {
#ifndef NDEBUG
  const auto& log = logger<LogArea::FindScatter>();
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
                           strings_with_probability_, nnbar_treatment_);

  // Cross section for collision criterion
  double cross_section_criterion = act->cross_section() * fm2_mb * M_1_PI /
                                   static_cast<double>(testparticles_);
  // Take cross section scaling factors into account
  cross_section_criterion *= data_a.xsec_scaling_factor(time_until_collision);
  cross_section_criterion *= data_b.xsec_scaling_factor(time_until_collision);

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
    const ParticleList& search_list, double dt) const {
  std::vector<ActionPtr> actions;
  for (const ParticleData& p1 : search_list) {
    for (const ParticleData& p2 : search_list) {
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
    const ParticleList& search_list, const ParticleList& neighbors_list,
    double dt) const {
  std::vector<ActionPtr> actions;
  for (const ParticleData& p1 : search_list) {
    for (const ParticleData& p2 : neighbors_list) {
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
    const ParticleList& search_list, const Particles& surrounding_list,
    double dt) const {
  std::vector<ActionPtr> actions;
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
            ScatterActionPtr act = make_unique<ScatterAction>(
                A, B, time, isotropic_, string_formation_time_);
            if (strings_switch_) {
              act->set_string_interface(string_process_interface_.get());
            }
            act->add_all_scatterings(elastic_parameter_, two_to_one_, incl_set_,
                                     low_snn_cut_, strings_switch_, use_AQM_,
                                     strings_with_probability_,
                                     nnbar_treatment_);
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
 * This data structure can be used to build a tree going from the initial state
 * (a collision of two particles) to all possible final states by recursively
 * performing all possible decays. The tree can be used to calculate the final
 * state cross sections.
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
   * Internal helper function for `print`, to be called recursively to print all
   * nodes.
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
 */
static void add_decays(Node& node) {
  // If there is more than one unstable particle in the current state, then
  // there will be redundant paths in the decay tree, corresponding to
  // reorderings of the decays. To avoid double counting, we normalize by the
  // number of possible decay orderings. Normalizing by the number of unstable
  // particles recursively corresponds to normalizing by the factorial that
  // gives the number of reorderings.
  //
  // Ideally, the redundant paths should never be added to the decay tree, but
  // we never have more than two redundant paths, so it probably does not matter
  // much.
  uint32_t n_unstable = 0;
  for (const ParticleTypePtr ptype : node.state_) {
    if (!ptype->is_stable()) {
      n_unstable += 1;
    }
  }
  const double norm =
      n_unstable != 0 ? 1. / static_cast<double>(n_unstable) : 1.;

  for (const ParticleTypePtr ptype : node.state_) {
    if (!ptype->is_stable()) {
      for (const auto& decay : ptype->decay_modes().decay_mode_list()) {
        ParticleTypePtrList parts;
        const auto name = make_decay_name(ptype->name(), decay, parts);
        auto& new_node = node.add_action(name, norm * decay->weight(), {ptype},
                                         std::move(parts));
        add_decays(new_node);
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

void ScatterActionsFinder::dump_cross_sections(const ParticleType& a,
                                               const ParticleType& b,
                                               double m_a, double m_b,
                                               bool final_state) const {
  typedef std::vector<std::pair<double, double>> xs_saver;
  std::map<std::string, xs_saver> xs_dump;
  std::map<std::string, double> outgoing_total_mass;

  ParticleData a_data(a), b_data(b);
  constexpr int n_momentum_points = 200;
  constexpr double momentum_step = 0.02;
  for (int i = 1; i < n_momentum_points + 1; i++) {
    const double momentum = momentum_step * i;
    a_data.set_4momentum(m_a, momentum, 0.0, 0.0);
    b_data.set_4momentum(m_b, -momentum, 0.0, 0.0);
    const double sqrts = (a_data.momentum() + b_data.momentum()).abs();
    const ParticleList incoming = {a_data, b_data};
    ScatterActionPtr act = make_unique<ScatterAction>(
        a_data, b_data, 0., isotropic_, string_formation_time_);
    if (strings_switch_) {
      act->set_string_interface(string_process_interface_.get());
    }
    act->add_all_scatterings(elastic_parameter_, two_to_one_, incl_set_,
                             low_snn_cut_, strings_switch_, use_AQM_,
                             strings_with_probability_, nnbar_treatment_);
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
        decaytree::add_decays(process_node);
      }
    }
    xs_dump["total"].push_back(std::make_pair(sqrts, act->cross_section()));
    // Total cross-section should be the first in the list -> negative mass
    outgoing_total_mass["total"] = -1.0;
    // tree.print();
    auto final_state_xs = tree.final_state_cross_sections();
    deduplicate(final_state_xs);
    for (const auto& p : final_state_xs) {
      outgoing_total_mass[p.name_] = p.mass_;
      xs_dump[p.name_].push_back(std::make_pair(sqrts, p.cross_section_));
    }
  }

  // Nice ordering of channels by summed pole mass of products
  std::vector<std::string> all_channels;
  for (const auto channel : xs_dump) {
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
  // Align everything to 16 unicode characters.
  // This should be enough for the longest channel name (7 final-state
  // particles).
  for (const auto channel : all_channels) {
    std::cout << utf8::fill_left(channel, 16, ' ');
  }
  std::cout << std::endl;

  // Print out all partial cross-sections in mb
  for (int i = 1; i < n_momentum_points; i++) {
    const double momentum = momentum_step * i;
    a_data.set_4momentum(m_a, momentum, 0.0, 0.0);
    b_data.set_4momentum(m_b, -momentum, 0.0, 0.0);
    const double sqrts = (a_data.momentum() + b_data.momentum()).abs();
    printf("%9.5f", sqrts);
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
      printf("%16.6f", xs);  // Same alignment as in the header.
    }
    printf("\n");
  }
}

}  // namespace smash
