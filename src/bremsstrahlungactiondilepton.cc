/*
 *
 *    Copyright (c) 2019-2020,2022,2024-2026
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/bremsstrahlungactiondilepton.h"

#include "smash/outputinterface.h"
#include "smash/random.h"

namespace smash {
static constexpr int LScatterAction = LogArea::ScatterAction::id;

/** The implementation of the BremsstrahlungActionDilepton class, which is a special
 * ScatterAction that produces a dilepton pair through bremsstrahlung. The final
 * state particles are not further propagated, only written to the dilepton output.
 */
// ── Constructor ──────────────────────────────────────────────────────────────
BremsstrahlungActionDilepton::BremsstrahlungActionDilepton(
    const ParticleList &in, const double time,
    const double hadronic_cross_section_input, FormFactorType ff_type)
    : ScatterAction(in[0], in[1], time),
      reac_(dilepton_brems_reaction_type(in)),
      hadronic_cross_section_(hadronic_cross_section_input),
      form_factor_type_(ff_type) {}

// ── Reaction type ─────────────────────────────────────────────────────────────      
BremsstrahlungActionDilepton::ReactionType
BremsstrahlungActionDilepton::dilepton_brems_reaction_type(const ParticleList &in) {
  /* Currently, only n+p dilepton bremsstrahlung is implemented. This function checks
   * if the incoming particles correspond to this process and returns the
   * corresponding enum. If the incoming particles do not correspond to any
   * implemented process, no_reaction is returned.
   */
  if (in.size() != 2) {
    return ReactionType::no_reaction;
  }

  PdgCode a = in[0].pdgcode();
  PdgCode b = in[1].pdgcode();

  switch (pack(a.code(), b.code())) {
    case (pack(pdg::p, pdg::n)):
    case (pack(pdg::n, pdg::p)):
      return ReactionType::np;

    default:
      return ReactionType::no_reaction;
  }
}

// ── Dummy hadronic process ────────────────────────────────────────────────────
// First, same pattern as BremsstrahlungAction::add_dummy_hadronic_process():
// The hadronic scatteraction (p+n elastic) is registered as a dummy
// to satisfy the ScatterAction machinery, while the actual dilepton
// emission is handled otherwise and does not rely on this machinery.
void BremsstrahlungActionDilepton::add_dummy_hadronic_process(
    double reaction_cross_section) {
  CollisionBranchPtr dummy_process = std::make_unique<CollisionBranch>(
      incoming_particles_[0].type(), incoming_particles_[1].type(),
      reaction_cross_section, ProcessType::BremsstrahlungDilepton);

  add_collision(std::move(dummy_process));

  // Second, define all outgoing particles at the end of the reaction,
  // i.e. the righthand side of the reaction np (-> np\gamma) -> npe⁺e⁻.
  // As this is the only process and there is no branching logic as
  // with the photon case, the add_single_process() function is not needed and 
  // the outgoing particles can be defined directly here.
  static const ParticleTypePtr e_p_particle = &ParticleType::find(pdg::e_p);
  static const ParticleTypePtr e_m_particle = &ParticleType::find(pdg::e_m);
  static const ParticleTypePtr p_particle = &ParticleType::find(pdg::p);
  static const ParticleTypePtr n_particle = &ParticleType::find(pdg::n);

  // Third, add the dilepton bremsstrahlung process branch to the final state list.
  CollisionBranchList final_state_list;

  // TODO: This part is commented out for now as reac_ is for sure np, 
  //       but it can be reworked in case we want to add more reactions in the future.
  //       In this case, the function needs to be changed as current return is void.
  /**
   * if (reac_ != ReactionType::np) {
   *     return final_state_list;
   * }
   */

  final_state_list.push_back(std::make_unique<CollisionBranch>(
      *p_particle, *n_particle,
      *e_p_particle, *e_m_particle,
      reaction_cross_section,
      ProcessType::BremsstrahlungDilepton));
  
  add_processes<CollisionBranch>(final_state_list, 
    collision_processes_dilepton_bremsstrahlung_,
    cross_section_dilepton_bremsstrahlung_);
}

// TODO: Go on from here next week. The following functions are still to be implemented:
// ── Perform dilepton bremsstrahlung ────────────────────────────────────────────────
/* The function perform_dilepton_bremsstrahlung is called to create the final state and write it
 * to the output. It first generates the final state by calling generate_final_state(), and then 
 * iterates over the list of outputs to find the dilepton output and writes the interaction to it.
*/
void BremsstrahlungActionDilepton::perform_dilepton_bremsstrahlung(const OutputsList &outputs) {
  // Compared to photon bremsstrahlung, only one photon is created. Hence, no loop anymore.
  generate_final_state();
  for (const auto &output : outputs) {
    // we only care about the dilepton output, the function will take care of this
    if (output->is_dilepton_output()) {
      // we do not care about the local density
      output->at_interaction(*this, 0.0);
    }
  }
}

// ── Main: generate_final_state ────────────────────────────────────────────────
void BremsstrahlungActionDilepton::generate_final_state() {
  //TODO: Rework completely for dileptons.
  // we have only one reaction per incoming particle pair
  if (collision_processes_dilepton_bremsstrahlung_.size() != 1) {
    logg[LScatterAction].fatal()
        << "Problem in BremsstrahlungActionDilepton::generate_final_state().\nThe "
           "process branch has "
        << collision_processes_dilepton_bremsstrahlung_.size()
        << " entries. It should however have 1.";
    throw std::runtime_error("");
  }

  auto *proc = collision_processes_dilepton_bremsstrahlung_[0].get();

  outgoing_particles_ = proc->particle_list();
  process_type_ = proc->get_type();
  FourVector interaction_point = get_interaction_point();

  // Sample k and theta:
  // minimum cutoff for k to be in accordance with cross section calculations
  double delta_k;  // k-range
  double k_min = 0.001;
  double k_max =
      (sqrt_s() * sqrt_s() - 2 * outgoing_particles_[0].type().mass() * 2 *
                                 outgoing_particles_[1].type().mass()) /
      (2 * sqrt_s());

  if ((k_max - k_min) < 0.0) {
    // Make sure it is kinematically even possible to create a photon that is
    // in accordance with the cross section cutoff
    k_ = 0.0;
    delta_k = 0.0;
  } else {
    k_ = random::uniform(k_min, k_max);
    delta_k = (k_max - k_min);
  }
  theta_ = random::uniform(0.0, M_PI);

  // Sample the phase space anisotropically in the local rest frame
  sample_3body_phasespace();

  // Get differential cross sections
  std::pair<double, double> diff_xs_pair = dilepton_brems_diff_cross_sections();
  double diff_xs_k = diff_xs_pair.first;
  double diff_xs_theta = diff_xs_pair.second;

  // Assign weighting factor
  const double W_theta = diff_xs_theta * (M_PI - 0.0);
  const double W_k = diff_xs_k * delta_k;
  weight_ = std::sqrt(W_theta * W_k) / hadronic_cross_section();

  // Set position and formation time and boost back to computational frame
  for (auto &new_particle : outgoing_particles_) {
    // assuming decaying particles are always fully formed
    new_particle.set_formation_time(time_of_execution_);
    new_particle.set_4position(interaction_point);
    new_particle.boost_momentum(
        -total_momentum_of_outgoing_particles().velocity());
  }

  // Not sure if the process here is similar to the photon process and also
  // not really part of the normal processes. A constant arbitrary number 
  // has been set in a similar fashion without interacting with the existing one.
  //TODO: Check with Ren how this is needed for the dilepton bremsstrahlung process 
  //      and if it can be set to a constant like this or we have to keep the 
  //      photon one due to the scatteractionphoton.cc implementation.
  const auto id_process = ID_PROCESS_DILEPTON_BREMS;
  Action::check_conservation(id_process);
}

// ── Stage 1: Sample k and theta, assign weight, set positions and boost ──────────
void BremsstrahlungActionDilepton::sample_3body_phasespace() {
  //TODO: Rework completely for dileptons.
  assert(outgoing_particles_.size() == 3);
  const double m_a = outgoing_particles_[0].type().mass(),
               m_b = outgoing_particles_[1].type().mass(),
               m_c = outgoing_particles_[2].type().mass();
  const double sqrts = sqrt_s();
  const double E_ab = sqrts - m_c - k_;  // Ekin of the pion pair in cm frame
  const double pcm = pCM(sqrts, E_ab, m_c);  // cm momentum of (π pair - photon)
  const double pcm_pions = pCM(E_ab, m_a, m_b);  // cm momentum within pion pair

  // Photon angle: Phi random, theta from theta_ sampled above
  const Angles phitheta_photon(random::uniform(0.0, twopi), std::cos(theta_));
  outgoing_particles_[2].set_4momentum(m_c, pcm * phitheta_photon.threevec());
  // Boost velocity to cm frame of the two pions
  const ThreeVector beta_cm_pion_pair_photon =
      pcm * phitheta_photon.threevec() / std::sqrt(pcm * pcm + E_ab * E_ab);

  // Sample pion pair isotropically
  Angles phitheta;
  phitheta.distribute_isotropically();
  outgoing_particles_[0].set_4momentum(m_a, pcm_pions * phitheta.threevec());
  outgoing_particles_[1].set_4momentum(m_b, -pcm_pions * phitheta.threevec());
  outgoing_particles_[0].boost_momentum(beta_cm_pion_pair_photon);
  outgoing_particles_[1].boost_momentum(beta_cm_pion_pair_photon);
}


// ── Stage 2: Sample the phase space anisotropically, get differential cross sections and assign weight ──────────
std::pair<double, double> BremsstrahlungActionDilepton::dilepton_brems_diff_cross_sections() {
  //TODO: Rework completely for dileptons.
  static const ParticleTypePtr pi_z_particle = &ParticleType::find(pdg::pi_z);
  const double collision_energy = sqrt_s();
  double dsigma_dk;
  double dsigma_dtheta;

  if (reac_ == ReactionType::pi_p_pi_m) {
    if (outgoing_particles_[0].type() != *pi_z_particle) {
      // pi+- + pi+-- -> pi+- + pi+- + gamma
      dsigma_dk =
          (*pipi_pipi_opp_dsigma_dk_interpolation)(k_, collision_energy);
      dsigma_dtheta = (*pipi_pipi_opp_dsigma_dtheta_interpolation)(
          theta_, collision_energy);
    } else {
      // pi+- + pi+-- -> pi0 + pi0 + gamma
      dsigma_dk = (*pipi_pi0pi0_dsigma_dk_interpolation)(k_, collision_energy);
      dsigma_dtheta =
          (*pipi_pi0pi0_dsigma_dtheta_interpolation)(theta_, collision_energy);
    }
  } else if (reac_ == ReactionType::pi_p_pi_p ||
             reac_ == ReactionType::pi_m_pi_m) {
    dsigma_dk = (*pipi_pipi_same_dsigma_dk_interpolation)(k_, collision_energy);
    dsigma_dtheta =
        (*pipi_pipi_same_dsigma_dtheta_interpolation)(theta_, collision_energy);
  } else if (reac_ == ReactionType::pi_z_pi_p ||
             reac_ == ReactionType::pi_z_pi_m) {
    dsigma_dk = (*pipi0_pipi0_dsigma_dk_interpolation)(k_, collision_energy);
    dsigma_dtheta =
        (*pipi0_pipi0_dsigma_dtheta_interpolation)(theta_, collision_energy);
  } else if (reac_ == ReactionType::pi_z_pi_z) {
    dsigma_dk = (*pi0pi0_pipi_dsigma_dk_interpolation)(k_, collision_energy);
    dsigma_dtheta =
        (*pi0pi0_pipi_dsigma_dtheta_interpolation)(theta_, collision_energy);
  } else {
    throw std::runtime_error(
        "Unkown channel when computing differential cross sections for "
        "bremsstrahlung processes.");
  }

  // Prevent negative cross sections due to numerics in interpolation
  dsigma_dk = (dsigma_dk < 0.0) ? really_small : dsigma_dk;
  dsigma_dtheta = (dsigma_dtheta < 0.0) ? really_small : dsigma_dtheta;

  // Combine differential cross sections to a pair
  std::pair<double, double> diff_x_sections = {dsigma_dk, dsigma_dtheta};

  return diff_x_sections;
}

}  // namespace smash
