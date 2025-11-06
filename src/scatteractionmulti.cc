/*
 *
 *    Copyright (c) 2020-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/scatteractionmulti.h"

#include <map>

#include "gsl/gsl_sf_ellint.h"

#include "smash/crosssections.h"
#include "smash/integrate.h"
#include "smash/logging.h"
#include "smash/parametrizations.h"
#include "smash/pow.h"

namespace smash {
static constexpr int LScatterActionMulti = LogArea::ScatterActionMulti::id;

ScatterActionMulti::ScatterActionMulti(
    const ParticleList& in_plist, double time,
    const SpinInteractionType spin_interaction_type)
    : Action(in_plist, time),
      total_probability_(0.),
      spin_interaction_type_(spin_interaction_type) {}

void ScatterActionMulti::add_reaction(CollisionBranchPtr p) {
  add_process<CollisionBranch>(p, reaction_channels_, total_probability_);
}

void ScatterActionMulti::add_reactions(CollisionBranchList pv) {
  add_processes<CollisionBranch>(std::move(pv), reaction_channels_,
                                 total_probability_);
}

double ScatterActionMulti::get_total_weight() const {
  double xsec_scaling = 1.0;
  for (const ParticleData& in_part : incoming_particles_) {
    xsec_scaling *= in_part.xsec_scaling_factor();
  }
  return total_probability_ * xsec_scaling;
}

double ScatterActionMulti::get_partial_weight() const {
  double xsec_scaling = 1.0;
  for (const ParticleData& in_part : incoming_particles_) {
    xsec_scaling *= in_part.xsec_scaling_factor();
  }
  return partial_probability_ * xsec_scaling;
}

void ScatterActionMulti::add_possible_reactions(
    double dt, const double gcell_vol,
    const MultiParticleReactionsBitSet incl_multi) {
  // 3 -> m
  if (incoming_particles_.size() == 3) {
    // 3 -> 1
    if (incl_multi[IncludedMultiParticleReactions::Meson_3to1] == 1) {
      if (three_different_pions(incoming_particles_[0], incoming_particles_[1],
                                incoming_particles_[2])) {
        // 3pi -> omega
        const ParticleTypePtr type_omega = ParticleType::try_find(0x223);
        if (type_omega) {
          add_reaction(std::make_unique<CollisionBranch>(
              *type_omega,
              probability_three_to_one(*type_omega, dt, gcell_vol,
                                       type_omega->spin_degeneracy()),
              ProcessType::MultiParticleThreeMesonsToOne));
        }
        // 3pi -> phi
        const ParticleTypePtr type_phi = ParticleType::try_find(0x333);
        if (type_phi) {
          add_reaction(std::make_unique<CollisionBranch>(
              *type_phi,
              probability_three_to_one(*type_phi, dt, gcell_vol,
                                       type_phi->spin_degeneracy()),
              ProcessType::MultiParticleThreeMesonsToOne));
        }
      } else if (two_pions_eta(incoming_particles_[0], incoming_particles_[1],
                               incoming_particles_[2])) {
        // eta2pi -> eta-prime
        const ParticleTypePtr type_eta_prime = ParticleType::try_find(0x331);

        int sym_factor_in = 1;
        if (incoming_particles_[0].type() == incoming_particles_[1].type() ||
            incoming_particles_[1].type() == incoming_particles_[2].type() ||
            incoming_particles_[2].type() == incoming_particles_[0].type()) {
          sym_factor_in = 2;  // 2 factorial
        }

        if (type_eta_prime) {
          add_reaction(std::make_unique<CollisionBranch>(
              *type_eta_prime,
              probability_three_to_one(
                  *type_eta_prime, dt, gcell_vol,
                  sym_factor_in * type_eta_prime->spin_degeneracy()),
              ProcessType::MultiParticleThreeMesonsToOne));
        }
      }
    }
    // 3 -> 2
    if (incl_multi[IncludedMultiParticleReactions::Deuteron_3to2] == 1) {
      const PdgCode pdg_a = incoming_particles_[0].pdgcode();
      const PdgCode pdg_b = incoming_particles_[1].pdgcode();
      const PdgCode pdg_c = incoming_particles_[2].pdgcode();
      const ParticleTypePtr type_deuteron =
          ParticleType::try_find(pdg::deuteron);
      const ParticleTypePtr type_anti_deuteron =
          ParticleType::try_find(pdg::antideuteron);

      const int spin_factor_inc = pdg_a.spin_degeneracy() *
                                  pdg_b.spin_degeneracy() *
                                  pdg_c.spin_degeneracy();

      if (type_deuteron && type_anti_deuteron) {
        // πpn → πd
        if ((pdg_a.is_pion() && pdg_b == pdg::p && pdg_c == pdg::n) ||
            (pdg_a.is_pion() && pdg_b == pdg::n && pdg_c == pdg::p) ||
            (pdg_a == pdg::p && pdg_b.is_pion() && pdg_c == pdg::n) ||
            (pdg_a == pdg::n && pdg_b.is_pion() && pdg_c == pdg::p) ||
            (pdg_a == pdg::p && pdg_b == pdg::n && pdg_c.is_pion()) ||
            (pdg_a == pdg::n && pdg_b == pdg::p && pdg_c.is_pion())) {
          // Get type of incoming π
          ParticleList::iterator it = std::find_if(
              incoming_particles_.begin(), incoming_particles_.end(),
              [](ParticleData x) { return x.is_pion(); });
          const ParticleType& type_pi = it->type();

          const double spin_degn =
              react_degen_factor(spin_factor_inc, type_pi.spin_degeneracy(),
                                 type_deuteron->spin_degeneracy());

          add_reaction(std::make_unique<CollisionBranch>(
              type_pi, *type_deuteron,
              probability_three_to_two(type_pi, *type_deuteron, dt, gcell_vol,
                                       spin_degn),
              ProcessType::MultiParticleThreeToTwo));
        }

        // πp̅n̅ → πd̅
        if ((pdg_a.is_pion() && pdg_b == -pdg::p && pdg_c == -pdg::n) ||
            (pdg_a.is_pion() && pdg_b == -pdg::n && pdg_c == -pdg::p) ||
            (pdg_a == -pdg::p && pdg_b.is_pion() && pdg_c == -pdg::n) ||
            (pdg_a == -pdg::n && pdg_b.is_pion() && pdg_c == -pdg::p) ||
            (pdg_a == -pdg::p && pdg_b == -pdg::n && pdg_c.is_pion()) ||
            (pdg_a == -pdg::n && pdg_b == -pdg::p && pdg_c.is_pion())) {
          // Get type of incoming π
          ParticleList::iterator it = std::find_if(
              incoming_particles_.begin(), incoming_particles_.end(),
              [](ParticleData x) { return x.is_pion(); });
          const ParticleType& type_pi = it->type();

          const double spin_degn =
              react_degen_factor(spin_factor_inc, type_pi.spin_degeneracy(),
                                 type_anti_deuteron->spin_degeneracy());

          add_reaction(std::make_unique<CollisionBranch>(
              type_pi, *type_anti_deuteron,
              probability_three_to_two(type_pi, *type_anti_deuteron, dt,
                                       gcell_vol, spin_degn),
              ProcessType::MultiParticleThreeToTwo));
        }

        // Nnp → Nd, N̅np → N̅d
        if ((pdg_a.is_nucleon() && pdg_b == pdg::p && pdg_c == pdg::n) ||
            (pdg_a.is_nucleon() && pdg_b == pdg::n && pdg_c == pdg::p) ||
            (pdg_a == pdg::p && pdg_b.is_nucleon() && pdg_c == pdg::n) ||
            (pdg_a == pdg::n && pdg_b.is_nucleon() && pdg_c == pdg::p) ||
            (pdg_a == pdg::p && pdg_b == pdg::n && pdg_c.is_nucleon()) ||
            (pdg_a == pdg::n && pdg_b == pdg::p && pdg_c.is_nucleon())) {
          int symmetry_factor = 1;  // already true for N̅np → N̅d case

          ParticleList::iterator it =
              std::find_if(incoming_particles_.begin(),
                           incoming_particles_.end(), [](ParticleData x) {
                             return x.pdgcode().antiparticle_sign() == -1;
                           });
          if (it == incoming_particles_.end()) {
            /* Meaning no anti-N found by find_if,
             * therefore not N̅np → N̅d, but Nnp → Nd. */
            symmetry_factor = 2;  // for Nnp → Nd (2 factorial)
            // It is already clear here that we have a double of two N
            if (pdg_a == pdg_b) {
              it = incoming_particles_.begin();
            } else {
              // If a and b are not the double, then c has to be part of it
              it = incoming_particles_.begin() + 2;
            }
          }
          const ParticleType& type_N = it->type();

          const double spin_degn =
              react_degen_factor(spin_factor_inc, type_N.spin_degeneracy(),
                                 type_deuteron->spin_degeneracy());

          add_reaction(std::make_unique<CollisionBranch>(
              type_N, *type_deuteron,
              probability_three_to_two(type_N, *type_deuteron, dt, gcell_vol,
                                       symmetry_factor * spin_degn),
              ProcessType::MultiParticleThreeToTwo));
        }

        // Np̅n̅ → Nd̅, N̅p̅n̅ → N̅d̅
        if ((pdg_a.is_nucleon() && pdg_b == -pdg::p && pdg_c == -pdg::n) ||
            (pdg_a.is_nucleon() && pdg_b == -pdg::n && pdg_c == -pdg::p) ||
            (pdg_a == -pdg::p && pdg_b.is_nucleon() && pdg_c == -pdg::n) ||
            (pdg_a == -pdg::n && pdg_b.is_nucleon() && pdg_c == -pdg::p) ||
            (pdg_a == -pdg::p && pdg_b == -pdg::n && pdg_c.is_nucleon()) ||
            (pdg_a == -pdg::n && pdg_b == -pdg::p && pdg_c.is_nucleon())) {
          int symmetry_factor = 1;  // already true for Np̅n̅ → Nd̅ case

          ParticleList::iterator it =
              std::find_if(incoming_particles_.begin(),
                           incoming_particles_.end(), [](ParticleData x) {
                             return x.pdgcode().antiparticle_sign() == 1;
                           });
          if (it == incoming_particles_.end()) {
            /* Meaning no N found by find_if,
             * therefore not Np̅n̅ → Nd̅, but N̅p̅n̅ → N̅d̅. */
            symmetry_factor = 2;  // for N̅p̅n̅ → N̅d̅ (2 factorial)
            // It is already clear here that we have a double of two N̅
            if (pdg_a == pdg_b) {
              it = incoming_particles_.begin();
            } else {
              // If a and b are not the double, then c has to be part of it
              it = incoming_particles_.begin() + 2;
            }
          }
          const ParticleType& type_N = it->type();

          const double spin_degn =
              react_degen_factor(spin_factor_inc, type_N.spin_degeneracy(),
                                 type_anti_deuteron->spin_degeneracy());

          add_reaction(std::make_unique<CollisionBranch>(
              type_N, *type_anti_deuteron,
              probability_three_to_two(type_N, *type_anti_deuteron, dt,
                                       gcell_vol, symmetry_factor * spin_degn),
              ProcessType::MultiParticleThreeToTwo));
        }
      }
    }
  }
  // 4 -> 2
  if (incoming_particles_.size() == 4 &&
      incl_multi[IncludedMultiParticleReactions::A3_Nuclei_4to2] == 1) {
    std::map<PdgCode, int> c;  // counts incoming PdgCodes
    int spin_factor_inc = 1;
    for (const ParticleData& data : incoming_particles_) {
      c[data.pdgcode()]++;
      spin_factor_inc *= data.pdgcode().spin_degeneracy();
    }
    // Nucleons, antinucleons, and pions can catalyze
    const int n_possible_catalysts_incoming =
        c[pdg::n] + c[pdg::p] + c[-pdg::p] + c[-pdg::n] + c[pdg::pi_p] +
        c[pdg::pi_z] + c[pdg::pi_m];

    for (PdgCode pdg_nucleus :
         {pdg::triton, pdg::antitriton, pdg::he3, pdg::antihe3,
          pdg::hypertriton, pdg::antihypertriton}) {
      const ParticleTypePtr type_nucleus = ParticleType::try_find(pdg_nucleus);
      // Nucleus can be formed if and only if:
      // 1) Incoming particles contain enough components (like p, n, Lambda)
      // 2) In (incoming particles - components) there is still a catalyst
      // This is including the situation like nnpp. Can be that t(nnp) is formed
      // and p is catalyst, can be that he-3(ppn) is formed and n is catalyst.
      // Both reactions should be added.
      const int n_nucleus_components_that_can_be_catalysts =
          pdg_nucleus.nucleus_p() + pdg_nucleus.nucleus_ap() +
          pdg_nucleus.nucleus_n() + pdg_nucleus.nucleus_an();
      const bool incoming_contain_nucleus_components =
          c[pdg::p] >= pdg_nucleus.nucleus_p() &&
          c[-pdg::p] >= pdg_nucleus.nucleus_ap() &&
          c[pdg::n] >= pdg_nucleus.nucleus_n() &&
          c[-pdg::n] >= pdg_nucleus.nucleus_an() &&
          c[pdg::Lambda] >= pdg_nucleus.nucleus_La() &&
          c[-pdg::Lambda] >= pdg_nucleus.nucleus_aLa();
      const bool can_form_nucleus =
          type_nucleus && incoming_contain_nucleus_components &&
          n_possible_catalysts_incoming -
                  n_nucleus_components_that_can_be_catalysts ==
              1;

      if (!can_form_nucleus) {
        continue;
      }
      // Find the catalyst
      std::map<PdgCode, int> catalyst_count = c;
      catalyst_count[pdg::p] -= pdg_nucleus.nucleus_p();
      catalyst_count[-pdg::p] -= pdg_nucleus.nucleus_ap();
      catalyst_count[pdg::n] -= pdg_nucleus.nucleus_n();
      catalyst_count[-pdg::n] -= pdg_nucleus.nucleus_an();
      catalyst_count[pdg::Lambda] -= pdg_nucleus.nucleus_La();
      catalyst_count[-pdg::Lambda] -= pdg_nucleus.nucleus_aLa();
      PdgCode pdg_catalyst = PdgCode::invalid();
      for (const auto i : catalyst_count) {
        if (i.second == 1) {
          pdg_catalyst = i.first;
          break;
        }
      }
      if (pdg_catalyst == PdgCode::invalid()) {
        logg[LScatterActionMulti].error("Something went wrong while forming",
                                        pdg_nucleus, " from ",
                                        incoming_particles_);
      }
      const ParticleTypePtr type_catalyst =
          ParticleType::try_find(pdg_catalyst);
      const double spin_degn =
          react_degen_factor(spin_factor_inc, type_catalyst->spin_degeneracy(),
                             type_nucleus->spin_degeneracy());
      double symmetry_factor = 1.0;
      for (const auto i : c) {
        symmetry_factor *= (i.second == 3)   ? 6.0  // 3!
                           : (i.second == 2) ? 2.0  // 2!
                                             : 1.0;
        if (i.second > 3 || i.second < 0) {
          logg[LScatterActionMulti].error("4<->2 error, incoming particles ",
                                          incoming_particles_);
        }
      }

      add_reaction(std::make_unique<CollisionBranch>(
          *type_catalyst, *type_nucleus,
          probability_four_to_two(*type_catalyst, *type_nucleus, dt, gcell_vol,
                                  symmetry_factor * spin_degn),
          ProcessType::MultiParticleFourToTwo));
    }
  }
  // 5 -> 2
  if (incoming_particles_.size() == 5) {
    if (incl_multi[IncludedMultiParticleReactions::NNbar_5to2] == 1) {
      if (all_incoming_particles_are_pions_have_zero_charge_only_one_piz()) {
        const int spin_factor_inc =
            incoming_particles_[0].pdgcode().spin_degeneracy() *
            incoming_particles_[1].pdgcode().spin_degeneracy() *
            incoming_particles_[2].pdgcode().spin_degeneracy() *
            incoming_particles_[3].pdgcode().spin_degeneracy() *
            incoming_particles_[4].pdgcode().spin_degeneracy();
        const double symmetry_factor = 4.0;  // 2! * 2!

        const ParticleTypePtr type_p = ParticleType::try_find(pdg::p);
        const ParticleTypePtr type_anti_p = ParticleType::try_find(-pdg::p);
        const ParticleTypePtr type_n = ParticleType::try_find(pdg::n);
        const ParticleTypePtr type_anti_n = ParticleType::try_find(-pdg::n);

        const double spin_degn =
            react_degen_factor(spin_factor_inc, type_p->spin_degeneracy(),
                               type_anti_p->spin_degeneracy());

        if (type_p && type_n) {
          const double prob = probability_five_to_two(
              type_p->mass(), dt, gcell_vol,
              symmetry_factor * spin_degn);  // same for ppbar and nnbar
          add_reaction(std::make_unique<CollisionBranch>(
              *type_p, *type_anti_p, prob,
              ProcessType::MultiParticleFiveToTwo));
          add_reaction(std::make_unique<CollisionBranch>(
              *type_n, *type_anti_n, prob,
              ProcessType::MultiParticleFiveToTwo));
        }
      }
    }
  }
}

void ScatterActionMulti::generate_final_state() {
  logg[LScatterActionMulti].debug("Incoming particles for ScatterActionMulti: ",
                                  incoming_particles_);

  /* Decide for a particular final state. */
  const CollisionBranch* proc =
      choose_channel<CollisionBranch>(reaction_channels_, total_probability_);
  process_type_ = proc->get_type();
  outgoing_particles_ = proc->particle_list();
  partial_probability_ = proc->weight();

  logg[LScatterActionMulti].debug("Chosen channel: ", process_type_,
                                  outgoing_particles_);

  switch (process_type_) {
    case ProcessType::MultiParticleThreeMesonsToOne:
      /* n->1 annihilation */
      annihilation();
      break;
    case ProcessType::MultiParticleThreeToTwo:
    case ProcessType::MultiParticleFourToTwo:
    case ProcessType::MultiParticleFiveToTwo:
      /* n->2 scattering */
      n_to_two();
      break;
    default:
      throw InvalidScatterActionMulti(
          "ScatterActionMulti::generate_final_state: Invalid process type " +
          std::to_string(static_cast<int>(process_type_)) + " was requested.");
  }

  /* The production point of the new particles.  */
  FourVector middle_point = get_interaction_point();

  for (ParticleData& new_particle : outgoing_particles_) {
    // Boost to the computational frame
    new_particle.boost_momentum(
        -total_momentum_of_outgoing_particles().velocity());
    /* Set positions of the outgoing particles */
    new_particle.set_4position(middle_point);
  }
}

double ScatterActionMulti::calculate_I3(const double sqrts) const {
  const double m1 = incoming_particles_[0].type().mass();
  const double m2 = incoming_particles_[1].type().mass();
  const double m3 = incoming_particles_[2].type().mass();

  if (sqrts < m1 + m2 + m3) {
    return 0.0;
  }
  const double x1 = (m1 - m2) * (m1 - m2), x2 = (m1 + m2) * (m1 + m2),
               x3 = (sqrts - m3) * (sqrts - m3),
               x4 = (sqrts + m3) * (sqrts + m3);
  const double qmm = x3 - x1, qmp = x3 - x2, qpm = x4 - x1, qpp = x4 - x2;
  const double kappa = std::sqrt(qpm * qmp / (qpp * qmm));
  const double tmp = std::sqrt(qmm * qpp);
  const double c1 =
      4.0 * m1 * m2 * std::sqrt(qmm / qpp) * (x4 - m3 * sqrts + m1 * m2);
  const double c2 = 0.5 * (m1 * m1 + m2 * m2 + m3 * m3 + sqrts * sqrts) * tmp;
  const double c3 = 8 * m1 * m2 / tmp *
                    ((m1 * m1 + m2 * m2) * (m3 * m3 + sqrts * sqrts) -
                     2 * m1 * m1 * m2 * m2 - 2 * m3 * m3 * sqrts * sqrts);
  const double c4 =
      -8 * m1 * m2 / tmp * smash::pow_int(sqrts * sqrts - m3 * m3, 2);
  const double res =
      c1 * gsl_sf_ellint_Kcomp(kappa, GSL_PREC_DOUBLE) +
      c2 * gsl_sf_ellint_Ecomp(kappa, GSL_PREC_DOUBLE) +
      c3 * gsl_sf_ellint_Pcomp(kappa, -qmp / qmm, GSL_PREC_DOUBLE) +
      c4 * gsl_sf_ellint_Pcomp(kappa, -x1 * qmp / (x2 * qmm), GSL_PREC_DOUBLE);
  return res;
}

double ScatterActionMulti::probability_three_to_one(
    const ParticleType& type_out, double dt, const double gcell_vol,
    const int degen_sym_factor) const {
  const double e1 = incoming_particles_[0].momentum().x0();
  const double e2 = incoming_particles_[1].momentum().x0();
  const double e3 = incoming_particles_[2].momentum().x0();
  const double sqrts = sqrt_s();

  const double gamma_decay = type_out.get_partial_width(
      sqrts, {&incoming_particles_[0].type(), &incoming_particles_[1].type(),
              &incoming_particles_[2].type()});

  const double I_3 = calculate_I3(sqrts);
  const double ph_sp_3 =
      1. / (8 * M_PI * M_PI * M_PI) * 1. / (16 * sqrts * sqrts) * I_3;

  const double spec_f_val = type_out.spectral_function(sqrts);

  return dt / (gcell_vol * gcell_vol) * M_PI / (4. * e1 * e2 * e3) *
         gamma_decay / ph_sp_3 * spec_f_val * std::pow(hbarc, 5.0) *
         degen_sym_factor;
}

double ScatterActionMulti::probability_three_to_two(
    const ParticleType& type_out1, const ParticleType& type_out2, double dt,
    const double gcell_vol, const double degen_sym_factor) const {
  const double e1 = incoming_particles_[0].momentum().x0();
  const double e2 = incoming_particles_[1].momentum().x0();
  const double e3 = incoming_particles_[2].momentum().x0();
  const double m4 = type_out1.mass();
  const double m5 = type_out2.mass();

  const double sqrts = sqrt_s();
  const double xs =
      CrossSections::two_to_three_xs(type_out1, type_out2, sqrts) / gev2_mb;
  const double lamb = lambda_tilde(sqrts * sqrts, m4 * m4, m5 * m5);

  const double I_3 = calculate_I3(sqrts);
  const double ph_sp_3 =
      1. / (8 * M_PI * M_PI * M_PI) * 1. / (16 * sqrts * sqrts) * I_3;

  return dt / (gcell_vol * gcell_vol) * 1. / (4. * e1 * e2 * e3) * lamb /
         (ph_sp_3 * 8 * M_PI * sqrts * sqrts) * xs * std::pow(hbarc, 5.0) *
         degen_sym_factor;
}

double ScatterActionMulti::probability_four_to_two(
    const ParticleType& type_out1, const ParticleType& type_out2, double dt,
    const double gcell_vol, const double degen_sym_factor) const {
  const double e1 = incoming_particles_[0].momentum().x0();
  const double e2 = incoming_particles_[1].momentum().x0();
  const double e3 = incoming_particles_[2].momentum().x0();
  const double e4 = incoming_particles_[3].momentum().x0();
  const double m5 = type_out1.mass();
  const double m6 = type_out2.mass();

  const double man_s = sqrt_s() * sqrt_s();
  const double xs =
      CrossSections::two_to_four_xs(type_out1, type_out2, sqrt_s()) / gev2_mb;
  const double lamb = lambda_tilde(man_s, m5 * m5, m6 * m6);
  const double ph_sp_4 = parametrizaton_phi4(man_s);

  return dt / std::pow(gcell_vol, 3.0) * 1. / (16. * e1 * e2 * e3 * e4) * xs /
         (4. * M_PI * man_s) * lamb / ph_sp_4 * std::pow(hbarc, 8.0) *
         degen_sym_factor;
}

double ScatterActionMulti::parametrizaton_phi4(const double man_s) const {
  int n_nucleons = 0, n_pions = 0, n_lambdas = 0;
  double sum_m = 0.0, prod_m = 1.0;
  for (const ParticleData& data : incoming_particles_) {
    const PdgCode pdg = data.type().pdgcode();
    n_nucleons += pdg.is_nucleon();  // including anti-nucleons
    n_pions += pdg.is_pion();
    n_lambdas += pdg.is_Lambda();  // including anti-Lambda
    sum_m += data.type().mass();
    prod_m *= data.type().mass();
  }
  const double x = 1.0 - sum_m / std::sqrt(man_s);
  const double x2 = x * x;
  const double x3 = x2 * x;
  double g = -1.0;

  if (n_nucleons == 3 && n_pions == 1) {  // NNNpi
    g = (1.0 + 0.862432 * x - 3.4853 * x2 + 1.70259 * x3) /
        (1.0 + 0.387376 * x - 1.34128 * x2 + 0.154489 * x3);
  } else if (n_nucleons == 4) {  // NNNN
    g = (1.0 - 1.72285 * x + 0.728331 * x2) /
        (1.0 - 0.967146 * x - 0.0103633 * x2);
  } else if (n_nucleons == 2 && n_lambdas == 1 && n_pions == 1) {  // LaNNpi
    g = (1.0 + 0.937064 * x - 3.56864 * x2 + 1.721 * x3) /
        (1.0 + 0.365202 * x - 1.2854 * x2 + 0.138444 * x3);
  } else if (n_nucleons == 3 && n_lambdas == 1) {  // LaNNN
    g = (1.0 + 0.882401 * x - 3.4074 * x2 + 1.62454 * x3) /
        (1.0 + 1.61741 * x - 2.12543 * x2 - 0.0902067 * x3);
  }

  if (g > 0.0) {
    return (std::sqrt(prod_m) * sum_m * sum_m * std::pow(x, 3.5) * g) /
           (840. * std::sqrt(2) * std::pow(M_PI, 4.0) * std::pow(1 - x, 4.0));
  } else {
    logg[LScatterActionMulti].error("parametrizaton_phi4: no parametrization ",
                                    "available for ", incoming_particles_);
    return 0.0;
  }
}

double ScatterActionMulti::probability_five_to_two(
    const double mout, double dt, const double gcell_vol,
    const double degen_sym_factor) const {
  const double e1 = incoming_particles_[0].momentum().x0();
  const double e2 = incoming_particles_[1].momentum().x0();
  const double e3 = incoming_particles_[2].momentum().x0();
  const double e4 = incoming_particles_[3].momentum().x0();
  const double e5 = incoming_particles_[4].momentum().x0();

  const double man_s = sqrt_s() * sqrt_s();
  const double lamb = lambda_tilde(man_s, mout * mout, mout * mout);
  const double ph_sp_5 = parametrizaton_phi5_pions(man_s);

  // Matching the NNbar anihilation cross section defintion for 2-to-5
  const double xs =
      std::max(0., ppbar_total(man_s) - ppbar_elastic(man_s)) / gev2_mb;

  return dt / std::pow(gcell_vol, 4.0) * 1. / (32. * e1 * e2 * e3 * e4 * e5) *
         xs / (4. * M_PI * man_s) * lamb / ph_sp_5 * std::pow(hbarc, 11.0) *
         degen_sym_factor;
}

double ScatterActionMulti::parametrizaton_phi5_pions(const double man_s) const {
  // see function documentation for parameter naming
  const double s_zero = 25 * pion_mass * pion_mass;
  const double fit_a = 2.1018e-10;
  const double fit_alpha = 1.982;
  return fit_a * std::pow(man_s - s_zero, 5.0) *
         std::pow(1 + man_s / s_zero, -fit_alpha);
}

void ScatterActionMulti::annihilation() {
  if (outgoing_particles_.size() != 1) {
    std::string s =
        "Annihilation: "
        "Incorrect number of particles in final state: ";
    s += std::to_string(outgoing_particles_.size()) + ".";
    throw InvalidScatterActionMulti(s);
  }
  // Set the momentum of the formed particle in its rest frame.
  outgoing_particles_[0].set_4momentum(
      total_momentum_of_outgoing_particles().abs(), 0., 0., 0.);
  // Make sure to assign formation times before boost to the computational frame
  assign_formation_time_to_outgoing_particles();
  if (spin_interaction_type_ != SpinInteractionType::Off) {
    assign_unpolarized_spin_vector_to_outgoing_particles();
  }

  logg[LScatterActionMulti].debug("Momentum of the new particle: ",
                                  outgoing_particles_[0].momentum());
}

void ScatterActionMulti::n_to_two() {
  sample_2body_phasespace();
  // Make sure to assign formation times before boost to the computational frame
  assign_formation_time_to_outgoing_particles();
  if (spin_interaction_type_ != SpinInteractionType::Off) {
    assign_unpolarized_spin_vector_to_outgoing_particles();
  }
  logg[LScatterActionMulti].debug(incoming_particles_.size(),
                                  "->2 scattering:", incoming_particles_,
                                  " -> ", outgoing_particles_);
}

bool ScatterActionMulti::three_different_pions(
    const ParticleData& data_a, const ParticleData& data_b,
    const ParticleData& data_c) const {
  // We want a combination of pi+, pi- and pi0
  const PdgCode pdg_a = data_a.pdgcode();
  const PdgCode pdg_b = data_b.pdgcode();
  const PdgCode pdg_c = data_c.pdgcode();

  return (pdg_a.is_pion() && pdg_b.is_pion() && pdg_c.is_pion()) &&
         (pdg_a != pdg_b && pdg_b != pdg_c && pdg_c != pdg_a);
}

bool ScatterActionMulti::two_pions_eta(const ParticleData& data_a,
                                       const ParticleData& data_b,
                                       const ParticleData& data_c) const {
  // We want a combination of pi0, pi0 and eta or pi+, pi- and eta
  const PdgCode pdg_a = data_a.pdgcode();
  const PdgCode pdg_b = data_b.pdgcode();
  const PdgCode pdg_c = data_c.pdgcode();

  return (pdg_a == pdg::pi_z && pdg_b == pdg::pi_z && pdg_c == pdg::eta) ||
         (pdg_a == pdg::pi_z && pdg_b == pdg::eta && pdg_c == pdg::pi_z) ||
         (pdg_a == pdg::eta && pdg_b == pdg::pi_z && pdg_c == pdg::pi_z) ||

         (pdg_a == pdg::eta && pdg_b == pdg::pi_m && pdg_c == pdg::pi_p) ||
         (pdg_a == pdg::eta && pdg_b == pdg::pi_p && pdg_c == pdg::pi_m) ||
         (pdg_a == pdg::pi_m && pdg_b == pdg::pi_p && pdg_c == pdg::eta) ||
         (pdg_a == pdg::pi_m && pdg_b == pdg::eta && pdg_c == pdg::pi_p) ||
         (pdg_a == pdg::pi_p && pdg_b == pdg::pi_m && pdg_c == pdg::eta) ||
         (pdg_a == pdg::pi_p && pdg_b == pdg::eta && pdg_c == pdg::pi_m);
}

bool ScatterActionMulti::
    all_incoming_particles_are_pions_have_zero_charge_only_one_piz() const {
  const bool all_inc_pi =
      all_of(incoming_particles_.begin(), incoming_particles_.end(),
             [](const ParticleData& data) { return data.is_pion(); });
  const int no_of_piz = std::count_if(
      incoming_particles_.begin(), incoming_particles_.end(),
      [](const ParticleData& data) { return data.pdgcode() == pdg::pi_z; });

  int total_state_charge = 0;
  for (const ParticleData& part : incoming_particles_) {
    total_state_charge += part.pdgcode().charge();
  }

  return (all_inc_pi && total_state_charge == 0 && no_of_piz == 1);
}

void ScatterActionMulti::format_debug_output(std::ostream& out) const {
  out << "MultiParticleScatter of " << incoming_particles_;
  if (outgoing_particles_.empty()) {
    out << " (not performed)";
  } else {
    out << " to " << outgoing_particles_;
  }
}

}  // namespace smash
