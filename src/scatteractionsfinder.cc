/*
 *
 *    Copyright (c) 2014-2020
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
#include "smash/cxx14compat.h"
#include "smash/decaymodes.h"
#include "smash/logging.h"
#include "smash/scatteraction.h"
#include "smash/scatteractionmulti.h"
#include "smash/scatteractionphoton.h"
#include "smash/stringfunctions.h"

namespace smash {
static constexpr int LFindScatter = LogArea::FindScatter::id;
/*!\Userguide
 * \page input_collision_term_ Collision_Term
 *
 * \key Elastic_Cross_Section (double, optional, default = -1.0 [mb]) \n
 * If a non-negative value is given, it will override the parametrized
 * elastic cross sections (which are energy-dependent) with a constant value.
 * This constant elastic cross section is used for all collisions.
 *
 * \key Isotropic (bool, optional, default = \key false) \n
 * Do all collisions isotropically.
 *
 * \key Maximum_Cross_Section (double, optional, default = 200 mb or 2000 mb
 *                             in case d' is present) \n
 * The maximal cross section for which it is guaranteed that all
 * collisions with this cross section will be found.
 * This means that all particle pairs, where the transverse distance
 * is smaller or equal to \f$ \sqrt{\sigma_{max}/\pi} \f$,
 * will be checked for collions.
 * The maximum around 200 mb occurs in the Delta peak of the pi+p
 * cross section. Many SMASH cross sections diverge close at the threshold,
 * these divergent parts are effectively cut off. If deuteron production
 * via d' is considered, then the default should be increased to 2000 mb
 * to function correctly (see \iref{Oliinychenko:2018ugs}). The maximal cross
 * section is scaled with \key Cross_Section_Scaling factor.
 *
 * \key Cross_Section_Scaling (double, optional, default = 1.0) \n
 * Scale all cross sections by a global factor. WARNING: Most cross sections are
 * constrained by experimental data. Scaling them will therefore lead to
 * unphysical results and is only meant for explorative studies.
 *
 * \key Additional_Elastic_Cross_Section (double, optional, default = 0.0 mb) \n
 * Add an additional constant contribution to the elastic cross section.
 * WARNING: Most elastic cross sections are constrained by experimental data.
 * Adding an additional contribution to them will therefore lead to unphysical
 * results and is only meant for explorative studies.
 *
 * \key Elastic_NN_Cutoff_Sqrts (double, optional, default = 1.98): \n
 * The elastic collisions betwen two nucleons with sqrt_s below
 * Elastic_NN_Cutoff_Sqrts, in GeV, cannot happen. \n
 * \li \key Elastic_NN_Cutoff_Sqrts < 1.88 - Below the threshold energy of the
 * elastic collsion, no effect \n
 * \li \key Elastic_NN_Cutoff_Sqrts > 2.02 - Beyond the threshold energy of the
 * inelastic collision NN->NNpi, not suggested
 *
 *
 * \key Strings (bool, optional, default = \key true for each setup except box):
 * \n \li \key true - String excitation is enabled\n \li \key false - String
 * excitation is disabled
 *
 * \key Collision_Criterion (string, optional, default = "Covariant"): \n
 * Choose collision criterion. For more information see
 * \subpage collision_criterion
 *
 * \key Only_Warn_For_High_Probability (bool, optional, default = \key false):
 * \n Only warn and not error for reaction probabilities higher than 1.
 * This switch is meant for very long production runs with the stochastic
 * criterion. It has no effect on the other criteria. If enabled the users
 * for themself have to make sure that the warning, that the probability has
 * slipped above 1, is printed very rarely.
 *
 * For information about more configuration options see the
 * following subpages \n
 * \li \subpage pauliblocker
 * \li \subpage string_parameters
 * \li \subpage input_dileptons
 * \li \subpage input_photons
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
 * \key Quark_Beta (double, optional, default = 7.0) \n
 * Parameter \f$\beta\f$ in PDF for quarks shown above.
 *
 * \key Strange_Supp (double, optional, default = 0.16) \n
 * Strangeness suppression factor \f$\lambda\f$:
 * \f[\lambda=\frac{P(s\bar{s})}{P(u\bar{u})}=\frac{P(s\bar{s})}{P(d\bar{d})}\f]
 * Defines the probability to produce a \f$s\bar{s}\f$ pair relative to produce
 * a light \f$q\bar{q}\f$ pair
 *
 * \key Diquark_Supp (double, optional, default = 0.036) \n
 * Diquark suppression factor. Defines the probability to produce a diquark
 * antidiquark pair relative to producing a qurk antiquark pair.
 *
 * \key Sigma_Perp (double, optional, default = 0.42 GeV) \n
 * Parameter \f$\sigma_\perp\f$ in distribution for transverse momentum
 * transfer between colliding hadrons \f$p_\perp\f$ and string mass \f$M_X\f$:
 * \f[\frac{d^3N}{dM^2_Xd^2\mathbf{p_\perp}}\propto \frac{1}{M_X^2}
 * \exp\left(-\frac{p_\perp^2}{\sigma_\perp^2}\right)\f]
 *
 * \key StringZ_A (double, optional, default = 2.0) \n
 * Parameter a in pythia fragmentation function \f$f(z)\f$:
 * \f[f(z) = \frac{1}{z} (1-z)^a \exp\left(-b\frac{m_T^2}{z}\right)\f]
 *
 * \key StringZ_B (double, optional, default = 0.55 1/GeV²) \n
 * Parameter \f$b\f$ in pythia fragmentation function shown above.
 *
 * \key Separate_Fragment_Baryon (bool, optional, default = True) \n
 * Whether to use a separate fragmentation function for leading baryons in
 * non-diffractive string processes.
 *
 * \key StringZ_A_Leading (double, optional, default = 0.2) \n
 * Parameter a in Lund fragmentation function used to sample the light cone
 * momentum fraction of leading baryons in non-diffractive string processes.
 *
 * \key StringZ_B_Leading (double, optional, default = 2.0 1/GeV²) \n
 * Parameter b in Lund fraghmentation function used to sample the light cone
 * momentum fraction of leading baryons in non-diffractive string processes.
 *
 * \key String_Sigma_T (double, optional, default = 0.5 GeV) \n
 * Standard deviation in Gaussian for transverse momentum distributed to
 * string fragments during fragmentation.
 *
 * \key Form_Time_Factor (double, optional, default = 1.0) \n
 * Factor to be multiplied with the formation time of string fragments from
 * the soft string routine.
 *
 * \key Power_Particle_Formation (double, optional, default = 1.0 or -1 for
 collisions at sqrt(s) = 200 GeV or higher) \n
 * If positive, the power with which the cross section scaling factor of
 * string fragments grows in time until it reaches 1. If negative, the scaling
 * factor will be constant and jump to 1 once the particle forms.
 *
 * \key Formation_Time (double, optional, default = 1.0 fm): \n
 * Parameter for formation time in string fragmentation, in fm/c.
 *
 * \key Mass_Dependent_Formation_Times (bool, optional, default = False) \n
 * Whether the formation time of string fragments should depend on their mass.
 * If it is set to true, the formation time is calculated as
 * \f$ \tau = \sqrt{2}\frac{m}{\kappa} \f$.
 *
 * \key Prob_proton_to_d_uu (double, optional, default = 1/3) \n
 * Probability of splitting an (anti)nucleon into the quark it has only once
 * and the diquark it contains twice in terms of flavour in the soft string
 * routine.
 *
 * \key Popcorn_Rate (double, optional, default = 0.15) \n
 * Parameter StringFlav:popcornRate, which determines production rate of
 * popcorn mesons in string fragmentation.
 * It is possible to produce a popcorn meson from the diquark end of a string
 * with certain probability (i.e., diquark to meson + diquark).
 *
 * **Examples: Configuring the String Paramters**\n
 *
 * String fragmentation is activated and if desired, the string parameters can
 * be altered as shown below.
 *
 *\verbatim
 Collision_Term:
     Strings: True
     String_Parameters:
         String_Tension: 1.0
         Gluon_Beta: 0.5
         Gluon_Pmin: 0.001
         Quark_Alpha: 2.0
         Quark_Beta: 7.0
         Strange_Supp: 0.16
         Diquark_Supp: 0.036
         Sigma_Perp: 0.42
         StringZ_A_Leading: 0.2
         StringZ_B_Leading: 2.0
         StringZ_A: 2.0
         StringZ_B: 0.55
         String_Sigma_T: 0.5
         Prob_proton_to_d_uu: 0.33
         Separate_Fragment_Baryon: True
         Popcorn_Rate: 0.15
  \endverbatim
 *
 *
 * \page collision_criterion Collision_Criterion
 * \key "Geometric" - Geometric collision criterion \n
 * The geometric collision criterion calculates the two-particle impact
 * parameter as the closest approach distance in the two-particle
 * center-of-momentum frame by boosting to the respective frame. The collision
 * time used for the ordering is calculated as the time of the closest approach
 * in the computational frame. \n For further details, see \iref{Bass:1998ca}.
 *
 * \key "Stochastic" - Stochastic collision criterion \n
 * The stochastic collision criterion employs a probability to decide whether
 * particles collide inside a given space-time cell. The probability is derived
 * directly from the scattering rate given by the Boltzmann equation. The
 * stochastic criterion is the only criterion that allows to treat
 * multi-particle reactions.\n
 * For more details, see e.g. Lang et al. (1993) \cite Lang1993.
 *
 * \note The stochastic criterion is only applicable within limits. For example,
 * it might not lead to reasonable results for very dilute systems like pp
 * collisions. Futhermore, the fixed time step mode is required. The
 * assumption for the criterion is that only one reaction per particle per
 * timestep occurs. Therefore, small enough timesteps (\key Delta_Time) have to
 * be used. In doubt, test if the results change with smaller timesteps. Since
 * the probability value is not by defintion limited to 1 in case of large
 * timesteps, an error is thrown if it gets larger than 1.
 *
 * \key "Covariant" - Covariant collision criterion \n
 * The covariant collision criterion uses a covariant expression of the
 * two-particle impact parameter in the two-particle center-of-momentum frame,
 * which allows for its calculation in the computational frame without boosting.
 * Furthermore, it calculates the collision times used for the collision
 * ordering in the two-particle center-of-momentum frame. \n
 * Further details are described in \iref{Hirano:2012yy}.
 *
 * \page input_collision_term_ Collision_Term
 * \n
 * **Examples: Configuring the Collision Term**\n
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
     Included_2to2:["Elastic","NN_to_NR","NN_to_DR","KN_to_KN","KN_to_KDelta"]
     Two_to_One: True
     Force_Decays_At_End: False
     NNbar_Treatment: "resonances"
     Elastic_Cross_Section: 30.0
     Elastic_NN_Cutoff_Sqrts: 1.93
     Isotropic: True
 \endverbatim
 *
 * If necessary, all collisions can be turned off by adding
 *\verbatim
     No_Collisions: True
 \endverbatim
 * in the configuration file. \n\n
 *
 * **Examples: Configuring deuteron multi-particle reactions**\n
 * The following example configures SMASH to includes deuteron multi-particle
 * reactions scatterings. Note, that all 2-to-2 reactions, in particular \key
 * "PiDeuteron_to_NN", are included except the d' reactions, since they
 * effectively yield the same reaction. Before using the example, if important,
 * check for completeness of all 2-to-2 reactions. The list of 2-to-2 reactions
 * might have grown and the example is therefore potentially out of date.
 *
 *\verbatim
 Collision_Term:
     Collision_Criterion: Stochastic
     # All (check for completeness, if important) 2-to-2 reactions except d'
     Included_2to2: ["Elastic","NN_to_NR", "NN_to_DR", "KN_to_KN",
                     "Strangeness_exchange", "NNbar", "PiDeuteron_to_NN"]
     Multi_Particle_Reactions: ["Deuteron_3to2"]
 \endverbatim
 *

 */

ScatterActionsFinder::ScatterActionsFinder(
    Configuration config, const ExperimentParameters& parameters)
    : coll_crit_(parameters.coll_crit),
      elastic_parameter_(
          config.take({"Collision_Term", "Elastic_Cross_Section"}, -1.)),
      testparticles_(parameters.testparticles),
      isotropic_(config.take({"Collision_Term", "Isotropic"}, false)),
      two_to_one_(parameters.two_to_one),
      incl_set_(parameters.included_2to2),
      incl_multi_set_(parameters.included_multi),
      scale_xs_(parameters.scale_xs),
      additional_el_xs_(parameters.additional_el_xs),
      low_snn_cut_(parameters.low_snn_cut),
      strings_switch_(parameters.strings_switch),
      use_AQM_(parameters.use_AQM),
      strings_with_probability_(parameters.strings_with_probability),
      nnbar_treatment_(parameters.nnbar_treatment),
      box_length_(parameters.box_length),
      string_formation_time_(config.take(
          {"Collision_Term", "String_Parameters", "Formation_Time"}, 1.)),
      maximum_cross_section_(parameters.maximum_cross_section),
      allow_first_collisions_within_nucleus_(
          parameters.allow_collisions_within_nucleus),
      only_warn_for_high_prob_(config.take(
          {"Collision_Term", "Only_Warn_For_High_Probability"}, false)) {
  if (is_constant_elastic_isotropic()) {
    logg[LFindScatter].info(
        "Constant elastic isotropic cross-section mode:", " using ",
        elastic_parameter_, " mb as maximal cross-section.");
  }
  if (incl_multi_set_.any() && coll_crit_ != CollisionCriterion::Stochastic) {
    throw std::invalid_argument(
        "Multi-body reactions (like e.g. 3->1 or 3->2) are only possible with "
        "the stochastic "
        "collision "
        "criterion. Change your config accordingly.");
  }

  if (incl_multi_set_[IncludedMultiParticleReactions::Deuteron_3to2] == 1 &&
      (incl_set_[IncludedReactions::PiDeuteron_to_pidprime] == 1 ||
       incl_set_[IncludedReactions::NDeuteron_to_Ndprime] == 1)) {
    throw std::invalid_argument(
        "To prevent double counting it is not possible to enable deuteron 3->2 "
        "reactions\nand reactions involving the d' at the same time\ni.e. to "
        "include `Deuteron_3to2` in `Multi_Particle_Reactions` and\n "
        "\"PiDeuteron_to_pidprime\" "
        "or \"NDeuteron_to_Ndprime\" in `Included_2to2` at the same time.\n"
        "Change your config accordingly.");
  }

  // TODO(stdnmr) Check that NNbar_5to2 is only set in sensible ways

  if (nnbar_treatment_ == NNbarTreatment::Resonances &&
      incl_set_[IncludedReactions::NNbar] != 1) {
    throw std::invalid_argument(
        "'NNbar' has to be in the list of allowed 2 to 2 processes "
        "to enable annihilation to go through resonances");
  }

  if (strings_switch_) {
    auto subconfig = config["Collision_Term"]["String_Parameters"];
    string_process_interface_ = make_unique<StringProcess>(
        subconfig.take({"String_Tension"}, 1.0), string_formation_time_,
        subconfig.take({"Gluon_Beta"}, 0.5),
        subconfig.take({"Gluon_Pmin"}, 0.001),
        subconfig.take({"Quark_Alpha"}, 2.0),
        subconfig.take({"Quark_Beta"}, 7.0),
        subconfig.take({"Strange_Supp"}, 0.16),
        subconfig.take({"Diquark_Supp"}, 0.036),
        subconfig.take({"Sigma_Perp"}, 0.42),
        subconfig.take({"StringZ_A_Leading"}, 0.2),
        subconfig.take({"StringZ_B_Leading"}, 2.0),
        subconfig.take({"StringZ_A"}, 2.0), subconfig.take({"StringZ_B"}, 0.55),
        subconfig.take({"String_Sigma_T"}, 0.5),
        subconfig.take({"Form_Time_Factor"}, 1.0),
        subconfig.take({"Mass_Dependent_Formation_Times"}, false),
        subconfig.take({"Prob_proton_to_d_uu"}, 1. / 3.),
        subconfig.take({"Separate_Fragment_Baryon"}, true),
        subconfig.take({"Popcorn_Rate"}, 0.15));
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
  if (!allow_first_collisions_within_nucleus_) {
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
  if (coll_crit_ == CollisionCriterion::Stochastic &&
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

  // Create ScatterAction object.
  ScatterActionPtr act = make_unique<ScatterAction>(
      data_a, data_b, time_until_collision, isotropic_, string_formation_time_,
      box_length_);

  if (coll_crit_ == CollisionCriterion::Stochastic) {
    act->set_stochastic_pos_idx();
  }

  if (strings_switch_) {
    act->set_string_interface(string_process_interface_.get());
  }

  // Distance squared calculation not needed for stochastic criterion
  const double distance_squared =
      (coll_crit_ == CollisionCriterion::Geometric)
          ? act->transverse_distance_sqr()
          : (coll_crit_ == CollisionCriterion::Covariant)
                ? act->cov_transverse_distance_sqr()
                : 0.0;

  // Don't calculate cross section if the particles are very far apart.
  // Not needed for stochastic criterion because of cell structure.
  if (coll_crit_ != CollisionCriterion::Stochastic &&
      distance_squared >= max_transverse_distance_sqr(testparticles_)) {
    return nullptr;
  }

  // Add various subprocesses.
  act->add_all_scatterings(elastic_parameter_, two_to_one_, incl_set_,
                           incl_multi_set_, low_snn_cut_, strings_switch_,
                           use_AQM_, strings_with_probability_,
                           nnbar_treatment_, scale_xs_, additional_el_xs_);

  double xs =
      act->cross_section() * fm2_mb / static_cast<double>(testparticles_);

  // Take cross section scaling factors into account
  xs *= data_a.xsec_scaling_factor(time_until_collision);
  xs *= data_b.xsec_scaling_factor(time_until_collision);

  if (coll_crit_ == CollisionCriterion::Stochastic) {
    const double v_rel = act->relative_velocity();
    /* Collision probability for 2-particle scattering, see e.g.
     * \iref{Xu:2004mz}, eq. (11) */
    const double prob = xs * v_rel * dt / gcell_vol;

    logg[LFindScatter].debug(
        "Stochastic collison criterion parameters (2-particles):\nprob = ",
        prob, ", xs = ", xs, ", v_rel = ", v_rel, ", dt = ", dt,
        ", gcell_vol = ", gcell_vol, ", testparticles = ", testparticles_);

    if (prob > 1.) {
      std::stringstream err;
      err << "Probability larger than 1 for stochastic rates. ( P_22 = " << prob
          << " )\nConsider using smaller timesteps.";
      if (only_warn_for_high_prob_) {
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

    ////////////////////////////////////////////////////////////////////////////

    // BLOCKING ALL 2-BODY BACKREACTIONS //

    // if (data_a.id_process() > 0 && data_a.id_process() == data_b.id_process()) {
    //   logg[LFindScatter].warn("# Accepted 2-body scattering originating from same process ( type: ",static_cast<int>(data_a.get_history().process_type),", id: ", data_a.id_process()," ). This should not happen often.");
    //   std::cout << data_a << '\n';
    //   std::cout << data_b << '\n';
    //   // std::cout << ". Blocking..." << '\n';
    //   // return nullptr;
    // }

    // BLOCKING ONLY 2-TO-5 BACKREACTIONS //

    if (data_a.id_process() > 0 && data_a.id_process() == data_b.id_process()) {
      if (data_a.get_history().process_type == ProcessType::MultiParticleFiveToTwo) {
        logg[LFindScatter].warn("# Blocking accepted 2-body scattering coming from same 5-to-2 process.");
        return nullptr;
      }
    }

    ////////////////////////////////////////////////////////////////////////////

  } else if (coll_crit_ == CollisionCriterion::Geometric ||
             coll_crit_ == CollisionCriterion::Covariant) {
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

  // Using std::move here is redundant with newer compilers, but required for
  // supporting GCC 4.8. Once we drop this support, std::move should be removed.
  return std::move(act);
}

ActionPtr ScatterActionsFinder::check_collision_multi_part(
    const ParticleList& plist, double dt, const double gcell_vol) const {
  /* If the two particles
   * 1) belong to the two colliding nuclei
   * 2) are within the same nucleus
   * 3) both of them have never experienced any collisons,
   * then the collision between them are banned also for multi-particle
   * interactions. */
  if (!allow_first_collisions_within_nucleus_) {
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
      make_unique<ScatterActionMulti>(plist, time_until_collision);

  act->set_stochastic_pos_idx();

  // 3. Add possible final states (dt and gcell_vol for probability calculation)
  act->add_possible_reactions(dt, gcell_vol, incl_multi_set_);

  /* 4. Return total collision probability
   *    Scales with 1 over the number of testpartciles to the power of the
   *    number of incoming particles - 1 */
  const double prob =
      act->get_total_weight() / std::pow(testparticles_, plist.size() - 1);

  // 5. Check that probability is smaller than one
  if (prob > 1.) {
    std::stringstream err;
    err << "Probability larger than 1 for stochastic rates. ( P_nm = " << prob
        << " )\nConsider using smaller timesteps.";
    if (only_warn_for_high_prob_) {
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

  //////////////////////////////////////////////////////////////////////////////


  // BLOCKING 5-BODY BACKREACTION WHEN ALL PARTICLES FROM SAME REACTION (INEFFECTIVE)  //

  bool all_5_particles_from_parent = false;
  if (same_parent_2to5_process && std::all_of(std::next(plist.begin()), plist.end(),
          [&](const ParticleData& data) { return plist[0].id_process() == data.id_process(); })) {
    all_5_particles_from_parent = true;
  }

  // if (all_5_particles_from_parent) {
  //   logg[LFindScatter].warn("# Accepted ", plist.size() ,"-body scattering originating from same process: ", plist[0].id_process(), ".");
  //   std::cout << ". Blocking..." << '\n';
  //   return nullptr;
  // }


  // BLOCKING 5-TO-2 BACKREACTION WHEN ANY PARTICLES FROM SAME REACTION (MORE EFFECTIVE)  //

  bool same_parent_2to5_process = false;
  for (int i = 0; i<plist.size(); i++) {
    for (int j = i+1; j<plist.size(); j++) {
      if (plist[i].id_process() == plist[j].id_process()) {
        // count_same_origin++;
        if (plist[i].get_history().process_type == ProcessType::TwoToFive) {
          same_parent_2to5_process = true;
          break;
        }
      }
    }
  }

  if (same_parent_2to5_process) {
    logg[LFindScatter].warn("# Blocking accepted 5-body process with particles coming from same 2-to-5 process.");
    if (all_5_particles_from_parent) {
      std::cout << "All from same parent." << '\n';
    }
    std::cout << "Parent ids: ";
    for (int i = 0; i<plist.size(); i++) {
      std::cout << plist[i].id_process() << " (type: " << static_cast<int>(plist[i].get_history().process_type) << ") ";
    }
    std::cout << '\n';
    return nullptr;
  }

  
  ////////////////////////////////////////////////////////////////////////////
  return std::move(act);
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
      if (incl_multi_set_.any()) {
        // Also, check for 3 particle scatterings with stochastic criterion
        for (const ParticleData& p3 : search_list) {
          if (incl_multi_set_[IncludedMultiParticleReactions::Deuteron_3to2] == 1 ||
              incl_multi_set_[IncludedMultiParticleReactions::Meson_3to1] == 1) {
            if (p1.id() < p2.id() && p2.id() < p3.id()) {
              ActionPtr act =
                  check_collision_multi_part({p1, p2, p3}, dt, gcell_vol);
              if (act) {
                actions.push_back(std::move(act));
              }
            }
          }
          // TODO(stdnmr) Introduce config option here
          if (incl_multi_set_[IncludedMultiParticleReactions::NNbar_5to2] == 1 &&
              search_list.size() >= 5) {
            for (const ParticleData& p4 : search_list) {
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
  if (coll_crit_ == CollisionCriterion::Stochastic) {
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
  if (coll_crit_ == CollisionCriterion::Stochastic) {
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
            ScatterActionPtr act = make_unique<ScatterAction>(
                A, B, time, isotropic_, string_formation_time_);
            if (strings_switch_) {
              act->set_string_interface(string_process_interface_.get());
            }
            act->add_all_scatterings(
                elastic_parameter_, two_to_one_, incl_set_, incl_multi_set_,
                low_snn_cut_, strings_switch_, use_AQM_,
                strings_with_probability_, nnbar_treatment_, scale_xs_,
                additional_el_xs_);
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
  /*
  // Round to output precision.
  for (auto& p : plab) {
    p = std::floor((p * 100000) + 0.5) / 100000;
  }
  */
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
    ScatterActionPtr act = make_unique<ScatterAction>(
        a_data, b_data, 0., isotropic_, string_formation_time_);
    if (strings_switch_) {
      act->set_string_interface(string_process_interface_.get());
    }
    act->add_all_scatterings(elastic_parameter_, two_to_one_, incl_set_,
                             incl_multi_set_, low_snn_cut_, strings_switch_,
                             use_AQM_, strings_with_probability_,
                             nnbar_treatment_, scale_xs_, additional_el_xs_);
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
  // Align everything to 16 unicode characters.
  // This should be enough for the longest channel name (7 final-state
  // particles).
  for (const auto& channel : all_channels) {
    std::cout << utf8::fill_left(channel, 16, ' ');
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
