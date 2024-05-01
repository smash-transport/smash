/*
 *
 *    Copyright (c) 2019-2020,2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/bremsstrahlungaction.h"

#include "smash/crosssectionsbrems.h"
#include "smash/outputinterface.h"
#include "smash/random.h"

namespace smash {
static constexpr int LScatterAction = LogArea::ScatterAction::id;

BremsstrahlungAction::BremsstrahlungAction(
    const ParticleList &in, const double time, const int n_frac_photons,
    const double hadronic_cross_section_input)
    : ScatterAction(in[0], in[1], time),
      reac_(bremsstrahlung_reaction_type(in)),
      number_of_fractional_photons_(n_frac_photons),
      hadronic_cross_section_(hadronic_cross_section_input) {}

BremsstrahlungAction::ReactionType
BremsstrahlungAction::bremsstrahlung_reaction_type(const ParticleList &in) {
  if (in.size() != 2) {
    return ReactionType::no_reaction;
  }

  PdgCode a = in[0].pdgcode();
  PdgCode b = in[1].pdgcode();

  switch (pack(a.code(), b.code())) {
    case (pack(pdg::pi_z, pdg::pi_m)):
    case (pack(pdg::pi_m, pdg::pi_z)):
      return ReactionType::pi_z_pi_m;

    case (pack(pdg::pi_z, pdg::pi_p)):
    case (pack(pdg::pi_p, pdg::pi_z)):
      return ReactionType::pi_z_pi_p;

    case (pack(pdg::pi_m, pdg::pi_p)):
    case (pack(pdg::pi_p, pdg::pi_m)):
      return ReactionType::pi_p_pi_m;

    case (pack(pdg::pi_m, pdg::pi_m)):
      return ReactionType::pi_m_pi_m;

    case (pack(pdg::pi_p, pdg::pi_p)):
      return ReactionType::pi_p_pi_p;

    case (pack(pdg::pi_z, pdg::pi_z)):
      return ReactionType::pi_z_pi_z;

    default:
      return ReactionType::no_reaction;
  }
}

void BremsstrahlungAction::perform_bremsstrahlung(const OutputsList &outputs) {
  for (int i = 0; i < number_of_fractional_photons_; i++) {
    generate_final_state();
    for (const auto &output : outputs) {
      if (output->is_photon_output()) {
        // we do not care about the local density
        output->at_interaction(*this, 0.0);
      }
    }
  }
}

void BremsstrahlungAction::generate_final_state() {
  // we have only one reaction per incoming particle pair
  if (collision_processes_bremsstrahlung_.size() != 1) {
    logg[LScatterAction].fatal()
        << "Problem in BremsstrahlungAction::generate_final_state().\nThe "
           "brocess branch has "
        << collision_processes_bremsstrahlung_.size()
        << " entries. It should however have 1.";
    throw std::runtime_error("");
  }

  auto *proc = collision_processes_bremsstrahlung_[0].get();

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
  std::pair<double, double> diff_xs_pair = brems_diff_cross_sections();
  double diff_xs_k = diff_xs_pair.first;
  double diff_xs_theta = diff_xs_pair.second;

  // Assign weighting factor
  const double W_theta = diff_xs_theta * (M_PI - 0.0);
  const double W_k = diff_xs_k * delta_k;
  weight_ = std::sqrt(W_theta * W_k) /
            (number_of_fractional_photons_ * hadronic_cross_section());

  // Scale weight by cross section scaling factor of incoming particles
  weight_ *= incoming_particles_[0].xsec_scaling_factor() *
             incoming_particles_[1].xsec_scaling_factor();

  // Set position and formation time and boost back to computational frame
  for (auto &new_particle : outgoing_particles_) {
    // assuming decaying particles are always fully formed
    new_particle.set_formation_time(time_of_execution_);
    new_particle.set_4position(interaction_point);
    new_particle.boost_momentum(
        -total_momentum_of_outgoing_particles().velocity());
  }

  // Set unpolarized spin vectors
  assign_unpolarized_spin_vector_to_outgoing_particles();

  // Photons are not really part of the normal processes, so we have to set a
  // constant arbitrary number.
  const auto id_process = ID_PROCESS_PHOTON;
  Action::check_conservation(id_process);
}

void BremsstrahlungAction::sample_3body_phasespace() {
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

void BremsstrahlungAction::add_dummy_hadronic_process(
    double reaction_cross_section) {
  CollisionBranchPtr dummy_process = std::make_unique<CollisionBranch>(
      incoming_particles_[0].type(), incoming_particles_[1].type(),
      reaction_cross_section, ProcessType::Bremsstrahlung);
  add_collision(std::move(dummy_process));
}

CollisionBranchList BremsstrahlungAction::brems_cross_sections() {
  CollisionBranchList process_list;
  // ParticleList final_state_particles;
  static const ParticleTypePtr photon_particle =
      &ParticleType::find(pdg::photon);
  static const ParticleTypePtr pi_z_particle = &ParticleType::find(pdg::pi_z);
  static const ParticleTypePtr pi_p_particle = &ParticleType::find(pdg::pi_p);
  static const ParticleTypePtr pi_m_particle = &ParticleType::find(pdg::pi_m);

  // Create interpolation objects, if not yet existent; only trigger for one
  // of them as either all or none is created
  if (pi0pi0_pipi_dsigma_dtheta_interpolation == nullptr) {
    create_interpolations();
  }

  // Find cross section corresponding to given sqrt(s)
  double sqrts = sqrt_s();
  double xsection;

  if (reac_ == ReactionType::pi_p_pi_m) {
    // Here the final state is determined by the the final state provided by the
    // sampled process using Monte Carlo techniqus

    // In the case of two oppositely charged pions as incoming particles,
    // there are two potential final states: pi+ + pi- and pi0 + pi0
    double xsection_pipi = (*pipi_pipi_opp_interpolation)(sqrts);
    double xsection_pi0pi0 = (*pipi_pi0pi0_interpolation)(sqrts);

    // Prevent negative cross sections due to numerics in interpolation
    xsection_pipi = (xsection_pipi <= 0.0) ? really_small : xsection_pipi;
    xsection_pi0pi0 = (xsection_pi0pi0 <= 0.0) ? really_small : xsection_pi0pi0;

    // Necessary only to decide for a final state with pi+ and pi- as incoming
    // particles.
    CollisionBranchList process_list_pipi;

    // Add both processes to the process_list
    process_list_pipi.push_back(std::make_unique<CollisionBranch>(
        incoming_particles_[0].type(), incoming_particles_[1].type(),
        *photon_particle, xsection_pipi, ProcessType::Bremsstrahlung));
    process_list_pipi.push_back(std::make_unique<CollisionBranch>(
        *pi_z_particle, *pi_z_particle, *photon_particle, xsection_pi0pi0,
        ProcessType::Bremsstrahlung));

    // Decide for one of the possible final states
    double total_cross_section = xsection_pipi + xsection_pi0pi0;
    const CollisionBranch *proc =
        choose_channel<CollisionBranch>(process_list_pipi, total_cross_section);

    xsection = proc->weight();

    process_list.push_back(std::make_unique<CollisionBranch>(
        proc->particle_list()[0].type(), proc->particle_list()[1].type(),
        *photon_particle, xsection, ProcessType::Bremsstrahlung));

  } else if (reac_ == ReactionType::pi_m_pi_m ||
             reac_ == ReactionType::pi_p_pi_p ||
             reac_ == ReactionType::pi_z_pi_m ||
             reac_ == ReactionType::pi_z_pi_p) {
    // Here the final state hadrons are identical to the initial state hadrons
    if (reac_ == ReactionType::pi_m_pi_m || reac_ == ReactionType::pi_p_pi_p) {
      xsection = (*pipi_pipi_same_interpolation)(sqrts);
    } else {
      // One pi0 in initial and final state
      xsection = (*pipi0_pipi0_interpolation)(sqrts);
    }

    // Prevent negative cross sections due to numerics in interpolation
    xsection = (xsection <= 0.0) ? really_small : xsection;

    process_list.push_back(std::make_unique<CollisionBranch>(
        incoming_particles_[0].type(), incoming_particles_[1].type(),
        *photon_particle, xsection, ProcessType::Bremsstrahlung));

  } else if (reac_ == ReactionType::pi_z_pi_z) {
    // Here we have a hard-coded final state that differs from the initial
    // state, namely: pi0 + pi0 -> pi+- + pi-+ + gamma
    xsection = (*pi0pi0_pipi_interpolation)(sqrts);

    // Prevent negative cross sections due to numerics in interpolation
    xsection = (xsection <= 0.0) ? really_small : xsection;

    process_list.push_back(std::make_unique<CollisionBranch>(
        *pi_p_particle, *pi_m_particle, *photon_particle, xsection,
        ProcessType::Bremsstrahlung));
  } else {
    throw std::runtime_error("Unknown ReactionType in BremsstrahlungAction.");
  }

  return process_list;
}

std::pair<double, double> BremsstrahlungAction::brems_diff_cross_sections() {
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

void BremsstrahlungAction::create_interpolations() {
  // Read in tabularized values for sqrt(s), k and theta
  std::vector<double> sqrts = BREMS_SQRTS;
  std::vector<double> photon_momentum = BREMS_K;
  std::vector<double> photon_angle = BREMS_THETA;

  // Read in tabularized total cross sections
  std::vector<double> sigma_pipi_pipi_opp = BREMS_PIPI_PIPI_OPP_SIG;
  std::vector<double> sigma_pipi_pipi_same = BREMS_PIPI_PIPI_SAME_SIG;
  std::vector<double> sigma_pipi0_pipi0 = BREMS_PIPI0_PIPI0_SIG;
  std::vector<double> sigma_pipi_pi0pi0 = BREMS_PIPI_PI0PI0_SIG;
  std::vector<double> sigma_pi0pi0_pipi = BREMS_PI0PI0_PIPI_SIG;

  // Read in tabularized differential cross sections dSigma/dk
  std::vector<double> dsigma_dk_pipi_pipi_opp = BREMS_PIPI_PIPI_OPP_DIFF_SIG_K;
  std::vector<double> dsigma_dk_pipi_pipi_same =
      BREMS_PIPI_PIPI_SAME_DIFF_SIG_K;
  std::vector<double> dsigma_dk_pipi0_pipi0 = BREMS_PIPI0_PIPI0_DIFF_SIG_K;
  std::vector<double> dsigma_dk_pipi_pi0pi0 = BREMS_PIPI_PI0PI0_DIFF_SIG_K;
  std::vector<double> dsigma_dk_pi0pi0_pipi = BREMS_PI0PI0_PIPI_DIFF_SIG_K;

  // Read in tabularized differential cross sections dSigma/dtheta
  std::vector<double> dsigma_dtheta_pipi_pipi_opp =
      BREMS_PIPI_PIPI_OPP_DIFF_SIG_THETA;
  std::vector<double> dsigma_dtheta_pipi_pipi_same =
      BREMS_PIPI_PIPI_SAME_DIFF_SIG_THETA;
  std::vector<double> dsigma_dtheta_pipi0_pipi0 =
      BREMS_PIPI0_PIPI0_DIFF_SIG_THETA;
  std::vector<double> dsigma_dtheta_pipi_pi0pi0 =
      BREMS_PIPI_PI0PI0_DIFF_SIG_THETA;
  std::vector<double> dsigma_dtheta_pi0pi0_pipi =
      BREMS_PI0PI0_PIPI_DIFF_SIG_THETA;

  // Create interpolation objects containing linear interpolations for
  // total cross sections
  pipi_pipi_opp_interpolation = std::make_unique<InterpolateDataLinear<double>>(
      sqrts, sigma_pipi_pipi_opp);
  pipi_pipi_same_interpolation =
      std::make_unique<InterpolateDataLinear<double>>(sqrts,
                                                      sigma_pipi_pipi_same);
  pipi0_pipi0_interpolation =
      std::make_unique<InterpolateDataLinear<double>>(sqrts, sigma_pipi0_pipi0);
  pipi_pi0pi0_interpolation =
      std::make_unique<InterpolateDataLinear<double>>(sqrts, sigma_pipi_pi0pi0);
  pi0pi0_pipi_interpolation =
      std::make_unique<InterpolateDataLinear<double>>(sqrts, sigma_pi0pi0_pipi);

  // Create interpolation objects containing bicubic interpolations for
  // differential dSigma/dk
  pipi_pipi_opp_dsigma_dk_interpolation =
      std::make_unique<InterpolateData2DSpline>(photon_momentum, sqrts,
                                                dsigma_dk_pipi_pipi_opp);
  pipi_pipi_same_dsigma_dk_interpolation =
      std::make_unique<InterpolateData2DSpline>(photon_momentum, sqrts,
                                                dsigma_dk_pipi_pipi_same);
  pipi0_pipi0_dsigma_dk_interpolation =
      std::make_unique<InterpolateData2DSpline>(photon_momentum, sqrts,
                                                dsigma_dk_pipi0_pipi0);
  pipi_pi0pi0_dsigma_dk_interpolation =
      std::make_unique<InterpolateData2DSpline>(photon_momentum, sqrts,
                                                dsigma_dk_pipi_pi0pi0);
  pi0pi0_pipi_dsigma_dk_interpolation =
      std::make_unique<InterpolateData2DSpline>(photon_momentum, sqrts,
                                                dsigma_dk_pi0pi0_pipi);

  // Create interpolation objects containing bicubic interpolations for
  // differential dSigma/dtheta
  pipi_pipi_opp_dsigma_dtheta_interpolation =
      std::make_unique<InterpolateData2DSpline>(photon_angle, sqrts,
                                                dsigma_dtheta_pipi_pipi_opp);
  pipi_pipi_same_dsigma_dtheta_interpolation =
      std::make_unique<InterpolateData2DSpline>(photon_angle, sqrts,
                                                dsigma_dtheta_pipi_pipi_same);
  pipi0_pipi0_dsigma_dtheta_interpolation =
      std::make_unique<InterpolateData2DSpline>(photon_angle, sqrts,
                                                dsigma_dtheta_pipi0_pipi0);
  pipi_pi0pi0_dsigma_dtheta_interpolation =
      std::make_unique<InterpolateData2DSpline>(photon_angle, sqrts,
                                                dsigma_dtheta_pipi_pi0pi0);
  pi0pi0_pipi_dsigma_dtheta_interpolation =
      std::make_unique<InterpolateData2DSpline>(photon_angle, sqrts,
                                                dsigma_dtheta_pi0pi0_pipi);
}
}  // namespace smash
