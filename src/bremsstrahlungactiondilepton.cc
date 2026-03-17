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
BremsstrahlungActionDilepton::dilepton_brems_reaction_type(const ParticleList &in)
 {
  /* Currently, only n+p dilepton bremsstrahlung is implemented. This function 
   * checks if the incoming particles correspond to this process and returns the
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

  // TODO: This part is commented out for now as reac_ is for sure 'np', 
  //       but it can be reworked in case we want to add more reactions in the future.
  //       In this case, the function needs to be changed as current return is void.
  /**
   * if (reac_ != ReactionType::np) {
   *     return final_state_list;
   * }
   */

  // For the 'np' reaction, the final state is 'pn e⁺e⁻' in this order. 
  // The order of particles written in CollisionBranch is relevant for the usage
  // in generate_final_state() to identify the outgoing_particles_ correctly.
  final_state_list.push_back(std::make_unique<CollisionBranch>(
      *p_particle, *n_particle,
      *e_p_particle, *e_m_particle,
      reaction_cross_section,
      ProcessType::BremsstrahlungDilepton));
  
  add_processes<CollisionBranch>(final_state_list, 
    collision_processes_dilepton_bremsstrahlung_,
    cross_section_dilepton_bremsstrahlung_);
}

// ── Perform dilepton bremsstrahlung ───────────────────────────────────────────
/* The function perform_dilepton_bremsstrahlung is called to create the final 
 * state and write it to the output. It first generates the final state by calling
 * generate_final_state(), and then iterates over the list of outputs to find the 
 * dilepton output and writes the interaction to it.
*/
void BremsstrahlungActionDilepton::perform_dilepton_bremsstrahlung(const OutputsList &outputs) 
{
  // Compared to photon bremsstrahlung, only one photon is created.
  // Hence, no loop anymore.
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
  // Sanity check: exactly one process branch expected
  if (collision_processes_dilepton_bremsstrahlung_.size() != 1) {
    logg[LScatterAction].fatal()
        << "Problem in BremsstrahlungActionDilepton::generate_final_state().\nThe "
           "process branch has "
        << collision_processes_dilepton_bremsstrahlung_.size()
        << " entries. It should however have 1.";
    throw std::runtime_error("");
  }
  // Get the process branch and define the outgoing particles from it.
  auto *proc = collision_processes_dilepton_bremsstrahlung_[0].get();

  outgoing_particles_ = proc->particle_list();
  process_type_ = proc->get_type();
  // Get the interaction point that is needed much later
  FourVector interaction_point = get_interaction_point();

  // Get the nucleon masses once to not get confused later.
  const double m_p = outgoing_particles_[0].type().mass();  // proton
  const double m_n = outgoing_particles_[1].type().mass();  // neutron

  // DOCUMENT: Fixed electron mass [GeV] included in constants.h now.
  //           Not sure if this would immediately work the same way as above via
  //           ParticleType, but to be seen. 
  const double m_e = outgoing_particles_[2].type().mass();  // electron mass

  // ── Step 1: Sample invariant mass M of dilepton pair ────────────────────────
  //
  // M runs from the kinematic threshold 2m_e up to the maximum value
  // allowed by 3-body kinematics: M_max = sqrt(s) - m_p - m_n.
  // Uniform sampling; physical distribution enters via the weight below.
  // (Analogous to k_ = random::uniform(k_min, k_max) in BremsstrahlungAction)

  double delta_M;
  // Minimum mass is 2m_e, as the dilepton pair has to be created on-shell.
  const double M_min = 2.0 * m_e;
  // Maximum mass is given by the kinematics of the 3-body final state.
  const double M_max = sqrt_s() - m_p - m_n;
  // Check if it is kinematically possible to create a dilepton pair.
  // If not, set the weight to 0 and return.
  if (M_max <= M_min) {
    weight_ = 0.0;
    return;
  }
  // If it is kinematically possible, sample M uniformly in the allowed range
  // and set delta_M accordingly. 
  else {
    M_ = random::uniform(M_min, M_max);
    delta_M = M_max - M_min;
  }

  // ── Step 2: Sample dilepton 3-momentum q ────────────────────────
  //
  // After fixing M, the momentum q depends on the CM energy sqrt_s.
  // There is a lower limit stemming from the kinematics of the 3-body final state. 
  // Given the sampled mass (which impacts the available phase space), 
  // the total energy in the CM frame must fit as well to create a dilepton pair. 
  // This translates into a minimum q_min that can be calculated from the 
  // kinematics of the 3-body final state. 
  
  // The upper limit q_max follows from the allowed kinematic range q²= E²-M² >= 0, 
  // which can be expressed in terms of the the implemented Källén function named
  // Action::lambda_tilde used in the 3-body kinematics. 
  // It holds q² = lambda_tilde(s, M², (m_p+m_n)²)/4s.

  // Calculate the lower limit of q from the kinematics of the 3-body final state 
  // given the sampled M_.
  double q_min;
  // Corresponding minimum energy E_min of the dilepton pair in the pn-CM frame 
  // given the beam energy.
  const double E_min = (sqrt_s()*sqrt_s() + M_*M_ - (m_p + m_n)*(m_p + m_n)) / 
    (2.0 * sqrt_s());
  // The minimum q_min is then given by sqrt(E_min² - M_²). If this is negative, 
  // set q_min to 0.
  if (E_min > M_) {
    q_min = std::sqrt(E_min * E_min - M_ * M_);
  } else {
    q_min = 0.0;
  }
  
  // Calculate q_max from the built-in pCM function.
  const double q_max = pCM(sqrt_s(), M_, m_p + m_n);
  
  // Sample q_ uniformly in [q_min, q_max] if kinematically allowed, otherwise set
  // the weight to 0 and return.
  double delta_q;
  if (q_max > q_min) {
    q_ = random::uniform(q_min, q_max);
    delta_q = q_max - q_min;
  } else {
    weight_ = 0.0;
    return;
  }

  // ── Step 3: Sample polar and azimuthal angles of dilepton in pn-CM frame ────
  //
  // Analogous to theta_ in BremsstrahlungAction but phi here explicitly as well.
  // Uniform in [0, pi]; azimuthal angle phi uniform in [0, 2pi].
  theta_ = random::uniform(0.0, M_PI);
  const double phi = random::uniform(0.0, twopi);

  // ── Step 4 (Stage 1): Construct dilepton 4-momentum in pn-CM frame ──────────
  //
  // The virtual photon carries 4-momentum p_ll with:
  //   |p_ll| = q,  E_ll = sqrt(q²+M²),  direction (theta, phi)
  const double E_ll = std::sqrt(q_ * q_ + M_ * M_);

  // Photon angle: phi, theta_ from above as in BremsstrahlungAction.
  const Angles phitheta(phi, std::cos(theta_));

  // Construct the dilepton 4-momentum in the pn-CM frame. 
  // The dilepton pair is treated as a single particle with mass M_ and momentum q_ 
  // in the direction given by phitheta.
  const FourVector p_ll(E_ll, q_ * phitheta.threevec());

  // To better distinguish the notation, let (pn)' denote the recoil nucleon pair 
  // after the collision.
  // Calculate the recoil for the (pn)' subsystem to check whether there is enough 
  // energy to create the outgoing pn pair. This is introduced compared 
  // to the photon case to ensure that the reaction is physically possible. 
  const FourVector p_recoil = total_momentum_of_outgoing_particles() - p_ll;

  // ── Step 5 (Stage 2): "Decay" of subsystems into e⁺e⁻ ───────────────────────
  //
  // The invariant mass of the (pn)' subsystem is sqrt(p_recoil²), so basically 
  // sqrt_s for (pn)' only.
  // If the invariant mass of (pn)' is smaller than the rest masses of p and n,
  // set the weight to 0 and return, as the action is energetically not possible.
  // Momentum of the (pn)' subsystem in its CM frame, calculated from the built-in
  // pCM function.
  // Within (pn)' subsystem, sample the angles of p and n isotropically...
  // ... and set the 4-momenta of the outgoing p and n in the (pn)' CM frame 
  // accordingly.
  // Obtain the velocity of the (pn)' subsystem relative to the pn-CM frame.
  // By construction of p_recoil above, the velocity of the (pn)' subsystem is 
  // opposite to the velocity of the dilepton pair. 
  // Boost the outgoing p and n from the moving (pn)' subsystem back to the 
  // pn-CM frame.
  // Note that the dilepton pair is, by construction, already in the pn-CM frame 
  // and no boost is needed for it here.
  sample_2body_isotropic(p_recoil, outgoing_particles_[0], outgoing_particles_[1]);

  // Isotropic 2-body decay in the virtual photon rest frame, then boost to 
  // pn-CM frame. This step is independent of the SPA formula and the PEFF —
  // it is pure kinematic bookkeeping.
  sample_2body_isotropic(p_ll, outgoing_particles_[2], outgoing_particles_[3]);

  // TODO: Continue here tomorrow with step 6.




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

// ── Stage 2: "Decay" of subsystem ─────────────────────────────────────────────
// This function assumes a 2-body decay of a parent particle into two daughters, 
// where the parent particle can be a subsystem of the reaction,
// e.g. the virtual photon in the dilepton case.
// The reaction is AB -> (AB)' + C and C can be a dilepton pair.
void BremsstrahlungActionDilepton::sample_2body_isotropic(
    const FourVector &p_parent, ParticleData &daughter1, ParticleData &daughter2)
{
  // The invariant mass of the (AB)' subsystem is sqrt(p_parent²), so basically 
  // sqrt_s for (AB)' only.
  const double M_parent = p_parent.abs();

  // If the invariant mass of the (AB)' subsystem is smaller than the sum of the 
  // rest masses of the daughters, set the weight to 0 and return, as the action 
  // is energetically not possible.
  if (M_parent < daughter1.type().mass() + daughter2.type().mass()) {
    weight_ = 0.0;
    return;
  }

  // Momentum of the 2-particle subsystem in its CM frame,
  // calculated from the built-in pCM function.
  const double pcm = pCM(M_parent, daughter1.type().mass(), 
    daughter2.type().mass());

  // Within 2-particle subsystem, sample the angles of daughters isotropically...
  Angles phitheta_daughters;
  phitheta_daughters.distribute_isotropically();
  // ... and set the 4-momenta of the outgoing daughters in the (AB)' CM frame 
  // accordingly.
  daughter1.set_4momentum(daughter1.type().mass(),
    pcm * phitheta_daughters.threevec());
  daughter2.set_4momentum(daughter2.type().mass(),
   -pcm * phitheta_daughters.threevec());
  // Obtain the velocity of the (AB)' subsystem relative to the AB-CM frame.
  // By construction of p_recoil above, the velocity of the (AB)' subsystem is 
  // opposite to the velocity of the dilepton pair.
  const ThreeVector beta = p_parent.velocity();
  // Boost the outgoing A and B from the moving (AB)' subsystem back to the 
  // AB-CM frame.
  // Note that the dilepton pair is, by construction, already in the AB-CM frame 
  // and no boost is needed for it here.
  daughter1.boost_momentum(-beta);
  daughter2.boost_momentum(-beta);
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
