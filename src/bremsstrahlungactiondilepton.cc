/*
 *
 *    Copyright (c) 2026
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/bremsstrahlungactiondilepton.h"

#include "smash/constants.h"
#include "smash/formfactors.h"
#include "smash/outputinterface.h"
#include "smash/parametrizations.h"
#include "smash/particledata.h"
#include "smash/particletype.h"
#include "smash/pdgcode.h"
#include "smash/pdgcode_constants.h"
#include "smash/pow.h"
#include "smash/random.h"

namespace smash {

namespace {

/**
 * Helper function for calculating R_2 as defined in \iref{Weil:2013mya},
 * eq. (42).
 *
 * \param[in] s Mandelstam variable s [GeV²]
 *
 * \return R_2(m_{inv}^2) (dimensionless)
 */
double R_2_helper(const double s) {
  const double m_pn = 2 * nucleon_mass;
  return (s < (m_pn * m_pn)) ? 0.0 : std::sqrt(1.0 - (m_pn * m_pn) / s);
}

}  // namespace

static constexpr int LScatterAction = LogArea::ScatterAction::id;

BremsstrahlungActionDilepton::BremsstrahlungActionDilepton(
    const ParticleList &in, const double time,
    const double hadronic_cross_section_input,
    DileptonBremsPionFormFactor ff_type)
    : ScatterAction(in[0], in[1], time),
      reaction_type_(dilepton_brems_reaction_type_(in)),
      hadronic_cross_section_(hadronic_cross_section_input),
      form_factor_type_(ff_type) {}

BremsstrahlungActionDilepton::ReactionType
BremsstrahlungActionDilepton::dilepton_brems_reaction_type_(
    const ParticleList &in) {
  if (in.size() != 2) {
    return ReactionType::no_reaction;
  }

  const PdgCode a = in[0].pdgcode();
  const PdgCode b = in[1].pdgcode();

  switch (pack(a.code(), b.code())) {
    case (pack(pdg::p, pdg::n)):
    case (pack(pdg::n, pdg::p)):
      return ReactionType::np;

    default:
      return ReactionType::no_reaction;
  }
}

void BremsstrahlungActionDilepton::add_dummy_hadronic_process(
    double reaction_cross_section) {
  CollisionBranchPtr dummy_process = std::make_unique<CollisionBranch>(
      incoming_particles_[0].type(), incoming_particles_[1].type(),
      reaction_cross_section, ProcessType::BremsstrahlungDilepton);

  add_collision(std::move(dummy_process));

  // Define all outgoing particles at the end of the reaction.
  static const ParticleTypePtr e_p_particle = &ParticleType::find(pdg::e_p);
  static const ParticleTypePtr e_m_particle = &ParticleType::find(pdg::e_m);
  static const ParticleTypePtr p_particle = &ParticleType::find(pdg::p);
  static const ParticleTypePtr n_particle = &ParticleType::find(pdg::n);

  CollisionBranchList final_state_list;

  assert(reaction_type_ != ReactionType::no_reaction);
  if (reaction_type_ == ReactionType::no_reaction) {
    logg[LScatterAction].fatal()
        << "Problem in " << __func__
        << ". Looks like an unknown reaction "
           "type for this dilepton bremsstrahlung process is present.";
    throw std::runtime_error(
        "Unreachable code was reached, "
        "please check logs for details.");
  }

  // For the 'np' reaction, the final state is 'pn e⁺e⁻' in this order.
  final_state_list.push_back(std::make_unique<CollisionBranch>(
      *p_particle, *n_particle, *e_p_particle, *e_m_particle,
      reaction_cross_section, ProcessType::BremsstrahlungDilepton));

  add_processes<CollisionBranch>(std::move(final_state_list),
                                 collision_processes_dilepton_bremsstrahlung_,
                                 cross_section_dilepton_bremsstrahlung_);
}

void BremsstrahlungActionDilepton::perform_dilepton_bremsstrahlung(
    const OutputsList &outputs) {
  // Only one photon is created per event.
  generate_final_state();
  for (const auto &output : outputs) {
    // we only care about the dilepton output, the function will take care
    if (output->is_dilepton_output()) {
      // we do not care about the local density
      output->at_interaction(*this, 0.0);
    }
  }
}

void BremsstrahlungActionDilepton::generate_final_state() {
  assert(collision_processes_dilepton_bremsstrahlung_.size() == 1);
  if (collision_processes_dilepton_bremsstrahlung_.size() != 1) {
    logg[LScatterAction].fatal()
        << "Problem in " << __func__
        << ". The function expects exactly "
           "one process branch for the dilepton bremsstrahlung process.";
    throw std::runtime_error(
        "Unreachable code was reached, "
        "please check logs for details.");
  }

  auto *proc = collision_processes_dilepton_bremsstrahlung_[0].get();

  outgoing_particles_ = proc->particle_list();
  process_type_ = proc->get_type();
  FourVector interaction_point = get_interaction_point();

  assert(outgoing_particles_.size() == 4);
  constexpr double m_p = nucleon_mass;
  constexpr double m_n = nucleon_mass;
  constexpr double m_e = electron_mass;

  const double M_min = 2.0 * m_e;
  const double M_max = sqrt_s() - m_p - m_n;
  // Check if it is possible to create a dilepton pair, i.e. M_min < M_max.
  if (M_max <= M_min) {
    weight_ = 0.0;
    return;
  } else {
    m_inv_ = random::uniform(M_min, M_max);
  }

  /* After fixing M, the momentum q depends on the CM energy sqrt_s.
   * There is no lower limit beside being positive, but the upper limit is given
   * by the kinematics of the 3-body final state.
   */
  const double q_min = 0.0;
  const double q_max = pCM(sqrt_s(), m_inv_, m_p + m_n);

  // Sample q_ uniformly in [q_min, q_max] if kinematically allowed.
  if (q_max > q_min) {
    q_ = random::uniform(q_min, q_max);
  } else {
    weight_ = 0.0;
    return;
  }

  // The virtual photon carries 4-momentum p_ll with direction (theta, phi)
  const double E_ll = std::sqrt(q_ * q_ + m_inv_ * m_inv_);

  Angles phitheta;
  phitheta.distribute_isotropically();
  const FourVector p_ll(E_ll, q_ * phitheta.threevec());

  // Calculate the recoil for the (pn)' subsystem after the collision.
  const FourVector p_recoil = FourVector(sqrt_s(), 0.0, 0.0, 0.0) - p_ll;

  /* If the invariant mass of (pn)' is smaller than the rest masses of p and n,
   * the sampling check fails and the function returns false, as the action is
   * energetically not possible. Hence, set the weight to 0 and return.
   */
  if (!sample_2body_isotropic_(p_recoil, outgoing_particles_[0],
                               outgoing_particles_[1])) {
    weight_ = 0.0;
    return;
  }

  /* Isotropic 2-body decay in the virtual photon rest frame, then boost to
   * pn-CM frame.
   */
  if (!sample_2body_isotropic_(p_ll, outgoing_particles_[2],
                               outgoing_particles_[3])) {
    weight_ = 0.0;
    return;
  }

  weight_ =
      diff_xs_pn_dilepton_(m_inv_, q_, sqrt_s()) / hadronic_cross_section_;

  weight_ *= incoming_particles_[0].xsec_scaling_factor() *
             incoming_particles_[1].xsec_scaling_factor();

  // Set positions and boost to computational frame
  for (auto &new_particle : outgoing_particles_) {
    new_particle.set_formation_time(time_of_execution_);
    new_particle.set_4position(interaction_point);
    new_particle.boost_momentum(
        -total_momentum_of_outgoing_particles().velocity());
  }
}

bool BremsstrahlungActionDilepton::sample_2body_isotropic_(
    const FourVector &p_parent, ParticleData &child_1, ParticleData &child_2) {
  const double M_parent = p_parent.abs();
  // Check whether the decay is energetically possible.
  if (M_parent < child_1.type().mass() + child_2.type().mass()) {
    return false;
  }

  const double pcm =
      pCM(M_parent, child_1.type().mass(), child_2.type().mass());

  Angles phitheta_children;
  phitheta_children.distribute_isotropically();
  child_1.set_4momentum(child_1.type().mass(),
                        pcm * phitheta_children.threevec());
  child_2.set_4momentum(child_2.type().mass(),
                        -pcm * phitheta_children.threevec());

  const ThreeVector beta = p_parent.velocity();
  child_1.boost_momentum(-beta);
  child_2.boost_momentum(-beta);

  return true;
}

double BremsstrahlungActionDilepton::diff_xs_pn_dilepton_(
    const double M, const double q, const double sqrts) const {
  const double s = sqrts * sqrts;
  constexpr double m_p = nucleon_mass;
  constexpr double m_n = nucleon_mass;
  const double m_pn = m_p + m_n;
  const double E = std::sqrt(q * q + M * M);

  // Check if the sampled point is kinematically allowed. If not, return 0.
  if (E <= 0.0 || M <= 0.0)
    return 0.0;

  const double sigma_bar =
      (s - (m_pn) * (m_pn)) / (2.0 * (m_p * m_p)) * np_elastic(s);
  const double R2_s = R_2_helper(s);
  if (R2_s <= 0.0)
    return 0.0;

  const double s2 = s + M * M - 2.0 * E * sqrts;
  const double R2_s2 = R_2_helper(s2);
  const double prefactor =
      fine_structure * fine_structure / (6.0 * M_PI * M_PI * M_PI);

  // Factor q²/(ME³) after dE->dq substitution in diff. cross section formula.
  return prefactor * (q * q) / (M * E * E * E) * sigma_bar * (R2_s2 / R2_s) *
         pion_em_form_factor_sq_(M);
}

double BremsstrahlungActionDilepton::pion_em_form_factor_sq_(
    const double M) const {
  const double rho_mass = ParticleType::find(pdg::rho_z).mass();
  const double rho_width = ParticleType::find(pdg::rho_z).total_width(M);

  switch (form_factor_type_) {
    case DileptonBremsPionFormFactor::FF1:
      return pion_em_form_factor_sqr_FF1(M * M, rho_mass, rho_width);
    case DileptonBremsPionFormFactor::FF2:
      return pion_em_form_factor_sqr_FF2(M * M, rho_mass, rho_width);
    case DileptonBremsPionFormFactor::Off:
      return 1.0;
    default:
      using namespace std::string_literals;  // NOLINT(build/namespaces)
      throw std::logic_error("Problem in "s + __func__ +
                             ". Unknown pion form factor.");
  }
}

}  // namespace smash
