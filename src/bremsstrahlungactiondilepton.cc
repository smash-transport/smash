/*
 *
 *    Copyright (c) 2019-2020,2022,2024-2026
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

// scatteraction.h as well as forwarddeclarations.h are included via main header 
// bremsstrahlungactiondilepton.h 
#include "smash/bremsstrahlungactiondilepton.h" // corresponding main header
#include "smash/formfactors.h" // for the form factors
#include "smash/constants.h" // for the lambda² value in the form factors
#include "smash/parametrizations.h"   // np_elastic(s)
#include "smash/particledata.h"
#include "smash/particletype.h"
#include "smash/pdgcode.h" // for pdg codes of particles
#include "smash/pdgcode_constants.h" 
#include "smash/pow.h" // for efficient power calculations
#include "smash/outputinterface.h" // for writing the dilepton output
#include "smash/random.h" // for random sampling of kinematic variables

namespace smash {
static constexpr int LScatterAction = LogArea::ScatterAction::id;

/** The implementation of the BremsstrahlungActionDilepton class, which is a special
 * ScatterAction that produces a dilepton pair through bremsstrahlung. The final
 * state particles are not further propagated, only written to the dilepton output.
 */
// ── Constructor ──────────────────────────────────────────────────────────────
BremsstrahlungActionDilepton::BremsstrahlungActionDilepton(
    const ParticleList &in, const double time,
    const double hadronic_cross_section_input, DileptonBremsFormFactor ff_type)
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
  
  add_processes<CollisionBranch>(std::move(final_state_list), 
    collision_processes_dilepton_bremsstrahlung_,
    cross_section_dilepton_bremsstrahlung_);
}

// ── Perform dilepton bremsstrahlung ───────────────────────────────────────────
/* The function perform_dilepton_bremsstrahlung is called to create the final 
 * state and write it to the output. It first generates the final state by calling
 * generate_final_state(), and then iterates over the list of outputs to find the 
 * dilepton output and writes the interaction to it.
*/
void BremsstrahlungActionDilepton::perform_dilepton_bremsstrahlung(
  const OutputsList &outputs) 
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
  // Note: Maybe it would be more physical to sample M according to the expected 
  //       distribution, but this would require some non-trivial changes in the 
  //       code and the weight calculation. Hence, left for future work and 
  //       keeping it simple for now.
  else {
    M_ = random::uniform(M_min, M_max);
    delta_M = M_max - M_min;
  }

  // ── Step 2: Sample dilepton 3-momentum q ────────────────────────
  //
  // After fixing M, the momentum q depends on the CM energy sqrt_s.
  // There is no lower limit beside being positive, but the upper limit is given 
  // by the kinematics of the 3-body final state. 
  
  // The upper limit q_max follows from the allowed kinematic range q²= E²-M² >= 0, 
  // which can be expressed in terms of the the implemented Källén function named
  // Action::lambda_tilde used in the 3-body kinematics. 
  // It holds q² = lambda_tilde(s, M², (m_p+m_n)²)/4s.

  // Set minimal momentum to zero.
  double q_min = 0.0;
  
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
  // Calculate total momentum via momentum conservation as sum of incoming momenta.
  // This is needed to calculate the recoil momentum of the (pn)' subsystem below.
  const FourVector total_momentum = incoming_particles_[0].momentum() + 
    incoming_particles_[1].momentum();
  // Calculate the recoil for the (pn)' subsystem to check whether there is enough 
  // energy to create the outgoing pn pair. This is introduced compared 
  // to the photon case to ensure that the reaction is physically possible. 
  const FourVector p_recoil = total_momentum - p_ll;

  // ── Step 5 (Stage 2): "Decay" of subsystems into e⁺e⁻ ───────────────────────
  //
  // The invariant mass of the (pn)' subsystem is sqrt(p_recoil²), so basically 
  // sqrt_s for (pn)' only.
  // If the invariant mass of (pn)' is smaller than the rest masses of p and n,
  // the sampling check fails and the function returns false, as the action is 
  // energetically not possible. Hence, set the weight to 0 and return, completely.
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
  if (!sample_2body_isotropic(p_recoil, outgoing_particles_[0], 
    outgoing_particles_[1])) {weight_ = 0.0; return;}

  // Isotropic 2-body decay in the virtual photon rest frame, then boost to 
  // pn-CM frame. This step is independent of the SPA formula and the PEFF —
  // it is pure kinematic bookkeeping but along the same lines as above.
  if (!sample_2body_isotropic(p_ll, outgoing_particles_[2], 
    outgoing_particles_[3])) {weight_ = 0.0; return;}

  // ── Step 6: Compute weight ──────────────────────────────────────────────────
  //
  // Integrand is dsigma/(dM dq dtheta dphi) evaluated at the sampled point.
  // Note on solid angle treatment:
  // We sample phi uniformly -> degree of freedom averages to 2pi.
  // We sample theta uniformly in [0,pi] -> As we sample 2pi * pi we get 2*pi².
  const double dsigma_dM_dq_dOmega = diff_xs_pn_dilepton(M_, q_, sqrt_s());
  const double W_M     = dsigma_dM_dq_dOmega * delta_M;
  const double W_q     = delta_q;
  const double W_Omega = twopi * M_PI;

  weight_ = W_M * W_q * W_Omega / hadronic_cross_section_;

  // xsec_scaling_factor() not applied here by intention:
  // BremsstrahlungActionDilepton solely covers elementary pn collisions, 
  // for which the scaling factor would always be equal to 1.
  // This is different to the treatment in BremsstrahlungAction (photons) 
  // as in that case there is a heavy ion context with formation time scaling.
  
  // Set positions and boost to computational frame
  for (auto &new_particle : outgoing_particles_) {
    new_particle.set_formation_time(time_of_execution_);
    new_particle.set_4position(interaction_point);
    new_particle.boost_momentum(
      -total_momentum_of_outgoing_particles().velocity());
  }

}

// ── Stage 2: "Decay" of subsystem ─────────────────────────────────────────────
// This function assumes a 2-body decay of a parent particle into two daughters, 
// where the parent particle can be a subsystem of the reaction, e.g. the virtual 
// photon in the dilepton case. The reaction is AB -> (AB)' + C and C can be a 
// dilepton pair.
// Function type is bool as it returns false if the decay is not kinematically 
// possible, true otherwise.
bool BremsstrahlungActionDilepton::sample_2body_isotropic(
    const FourVector &p_parent, ParticleData &daughter1, ParticleData &daughter2)
{
  // The invariant mass of the (AB)' subsystem is sqrt(p_parent²), so basically 
  // sqrt_s for (AB)' only.
  const double M_parent = p_parent.abs();

  // If the invariant mass of the (AB)' subsystem is smaller than the sum of the 
  // rest masses of the daughters, return false, as the action 
  // is energetically not possible.
  if (M_parent < daughter1.type().mass() + daughter2.type().mass()) {
    return false;
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
  
  // If the function has not returned false until here, the decay is kinematically
  // possible and the function can return true.
  return true;
}

// ── Differential cross section: SPA + PEFF ────────────────────────────────────
//
double BremsstrahlungActionDilepton::diff_xs_pn_dilepton(
    const double M, const double q, const double sqrts) const {
  //
  // M: invariant mass of dilepton pair
  // q: 3-momentum of dilepton in pn-CM frame
  // E = sqrt(q²+M²): dilepton energy in pn-CM frame
  // s: Mandelstam-s of pn collision
  // m_p: proton mass
  // m_n: neutron mass
  const double s  = sqrts * sqrts;
  const double m_p = outgoing_particles_[0].type().mass();  // proton
  const double m_n = outgoing_particles_[1].type().mass();  // neutron
  const double m_pn = m_p + m_n;
  const double E  = std::sqrt(q * q + M * M);

  // Check if the sampled point is kinematically allowed. If not, return 0.
  if (E <= 0.0 || M <= 0.0) return 0.0;

  // sigma_bar(s): Eq. (41) in Weil (2013), p. 23
  const double sigma_bar   = (s - (m_pn)*(m_pn)) / (2.0 * (m_p * m_p)) 
                              * np_elastic(s);

  // R2(s): Eq. (42) in Weil (2013), p. 23 with truncation at 0 to prevent 
  // numerical issues just in case
  const double R2_s = R_2_helper(s);
  // Check for potential division by zero in the SPA formula.
  if (R2_s <= 0.0) return 0.0;

  // s2: Eq. (43) in Weil (2013), p. 23
  const double s2   = s + M * M - 2.0 * E * sqrts;

  // R2(s2): same factor evaluated at reduced energy
  const double R2_s2 = R_2_helper(s2);

  // Prefactor alpha²/(6pi³)
  const double prefactor = fine_structure * fine_structure / 
                            (6.0 * M_PI * M_PI * M_PI);

  // Factor q²/(ME³)  [after dE -> dq substitution]
  const double factor_2 = (q * q) / (M * E * E * E);

  // Total differential cross section: Eq. (40) in Weil (2013), p. 23, without
  // the PEFF correction, which is applied later as a multiplicative factor.
  double diff_xs = prefactor * factor_2 * sigma_bar * (R2_s2 / R2_s);

  // ── Pion electromagnetic form factor (Shyam & Mosel 2010) ───────────────────
  //
  // The PEFF modifies the QM bremsstrahlung contribution via the
  // internal charged pion propagator in diagram 1(c) of the pn-channel.
  // The differential cross-section goes dsigma/dM -> dsigma/dM * |F_pi(M²)|²
  diff_xs *= pion_em_form_factor_sq(M * M);

  return diff_xs;
}

// ── Helper function R_2(s) ────────────────────────────────────────────────────
//
double BremsstrahlungActionDilepton::R_2_helper(const double s) const {
  const double m_p = outgoing_particles_[0].type().mass();  // proton
  const double m_n = outgoing_particles_[1].type().mass();  // neutron
  const double m_pn = m_p + m_n;
  if (s <= (m_pn * m_pn)) {
    return 0.0;
  } else {
    return std::sqrt(1.0 - (m_pn * m_pn) / s);
  }
}

// ── Form factor |F_pi(M²)|² (Shyam & Mosel 2010) ──────────────────────────────
//
double BremsstrahlungActionDilepton::pion_em_form_factor_sq(
    const double M_sq) const {
  // Obtain rho pole mass from particles.txt
  const double m_rho = ParticleType::find(pdg::rho_z).mass();
  // Energy-dependent rho width
  const double Gamma = gamma_rho(M_sq);

  switch (form_factor_type_)
  {
  case DileptonBremsFormFactor::FF1:
    return pion_em_form_factor_sqr_FF1(M_sq, m_rho, Gamma);
  case DileptonBremsFormFactor::FF2:
    return pion_em_form_factor_sqr_FF2(M_sq, m_rho, Gamma);
  case DileptonBremsFormFactor::Off:
  // If no form factor is applied, return 1.
    return 1.0;
  default:
  // TODO: Change this to throw an error if the form factor type is not recognized,
  // as there should be no case where this happens.
    return 1.0;
  }
}

// ── Energy-dependent Gamma_rho ────────────────────────────────────────────────
//
double BremsstrahlungActionDilepton::gamma_rho(double M_sq) const {
  // The energy-dependent width Gamma_rho(M²) as provided in Brown et.al (1986)
  // is given by the formula:
  // Gamma_rho(M²) = gamma0_rho * m_rho³ / [M * (2 * m_rho² - M²)]
  //                * (M² - 4 * m_pi²)^(3/2) / (m_rho² - 4 * m_pi²)^(3/2)
  // with 
  // m_rho:       rho pole mass, 
  // gamma0_rho:  width at rho pole mass,
  // m_pi:        pion mass.

  // Get masses and width from particles.txt
  const double m_rho = ParticleType::find(pdg::rho_z).mass();
  const double m_rho_sq = ParticleType::find(pdg::rho_z).mass_sqr();
  const double m_pi_sq = ParticleType::find(pdg::pi_z).mass_sqr();
  const double gamma0_rho = ParticleType::find(pdg::rho_z).width_at_pole();
  
  // Check if M² is above the 2-pion threshold. If not, return the width at the 
  // rho pole mass. 
  const double num = M_sq - 4.0 * m_pi_sq;
  if (num <= 0) return gamma0_rho;

  // If M² is above the 2-pion threshold, calculate the energy-dependent width  
  const double M          = std::sqrt(M_sq);
  const double num_sqrt   = std::sqrt(num);
  const double denom_sqrt = std::sqrt(m_rho_sq - 4.0 * m_pi_sq);

  return gamma0_rho * (m_rho_sq * m_rho) / (M * (2.0 * m_rho_sq - M_sq))
         * pow_int(num_sqrt / denom_sqrt, 3);
}

}  // namespace smash
