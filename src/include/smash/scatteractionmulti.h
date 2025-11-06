/*
 *
 *    Copyright (c) 2020-2022,2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_SCATTERACTIONMULTI_H_
#define SRC_INCLUDE_SMASH_SCATTERACTIONMULTI_H_

#include "action.h"

namespace smash {

/**
 * \ingroup action
 * ScatterActionMulti is a special action which takes any number of incoming
 * particles and performs a scattering with the use of the stochastic criterion,
 * producing one or more final-state particles.
 */
class ScatterActionMulti : public Action {
 public:
  /**
   * Construct a ScatterActionMulti object.
   *
   * \param[in] in_plist List of reaction partners
   * \param[in] time Time at which the action is supposed to take place
   * \param[in] spin_interaction_type Type of spin interaction
   */
  ScatterActionMulti(const ParticleList& in_plist, double time,
                     const SpinInteractionType spin_interaction_type =
                         SpinInteractionType::Off);

  /**
   * Generate the final-state of the multi-particle scattering process.
   * Assign position and momenta to outgoing particles.
   *
   * \throw InvalidScatterActionMulti
   */
  void generate_final_state() override;

  /**
   * Get the total probability for the reaction (scaled with the cross section
   * scaling factors of the incoming particles).
   *
   * \return total probability.
   */
  double get_total_weight() const override;

  /**
   * Get the partial probability for the chosen channel (scaled with the cross
   * section scaling factors of the incoming particles).
   *
   * \return partial probability.
   */
  double get_partial_weight() const override;

  /**
   * Add all possible multi-particle reactions for the given incoming particles.
   *
   * \param[in] dt timestep size
   * \param[in] gcell_vol gcell_vol grid cell volume
   * \param[in] incl_multi Which multi-particle reactions are enabled?
   */
  void add_possible_reactions(double dt, const double gcell_vol,
                              const MultiParticleReactionsBitSet incl_multi);

  /**
   * Get list of possible reaction channels.
   *
   * \return list of possible reaction channels.
   */
  const CollisionBranchList& reaction_channels() { return reaction_channels_; }

  /**
   * \ingroup exception
   * Thrown when ScatterActionMulti is called to perform with unknown
   * combination of incoming and outgoing number of particles or unknown process
   * type.
   */
  class InvalidScatterActionMulti : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };

  /**
   * Calculate the integration necessary for the three-body phase space. The
   * defintion for the integral is given by
   * \f[I_3 = \int dm^2_{23}dm^2_{12} =
   * \int^{(M-m_3)^2}_{(m_1+m_2)^2}[m^2_{23, max}- m^2_{23, min}]dm^2_{12}\f]
   * see PDG book (chapter Kinematics) for defintions of variables. The numbered
   * masses reference the incoming particles and \f$M\f$ the mass of the
   * outgoing particles in this case, since we are looking at the backreaction
   * to the 1-to-3 decay.
   *
   * \param[in] sqrts center of mass energy of incoming particles
   *                  (= mass of outgoing particle)
   * \return result of integral
   */
  double calculate_I3(const double sqrts) const;

  /**
   * Calculate the parametrized 4-body relativistic phase space integral.
   *
   * The 4-body phase space is a n = 4 case of general integral over
   * \f[d\Phi_n = (2\pi)^4 \prod_{i=1}^n \frac{d^3p_i}{2E_i (2\pi)^3}
   *    \times \delta(E_{tot} -  \sum E_i) \delta^{(3)}(p_{tot} - \sum p_i)\f].
   * This is a Lorentz-invariant quantity, so the result of the integration
   * depends only on \f[ s = E_{tot}^2 - p_{tot}^2\f] and masses of the
   * particles \f[ m_i^2 = E_i^2 - p_i^2 \f]. The dimension in general case is
   * \f[\mathrm{GeV}^{2n-4}\f]. The \f[\hbar = c = 1 \f] convention is used
   * here, so when this integral is used in the acceptance probability for
   * collision, one has to restore \f[\hbar\f] to obtain correct dimensionless
   * probability. More on phase space integrals can be found, for example, in
   * CERN-68-15 report. For developers, I (oliiny) have compiled a document with
   * properties and parametrizations of many-body phase space integrals here:
   * github.com/smash-transport/smash-devel/files/7791360/n_body_relativistic_phase_space.pdf
   *
   * \param[in] man_s mandelstam s of reaction
   * \return phase space integral value for 4 bodies [GeV^4]
   */
  double parametrizaton_phi4(const double man_s) const;

 protected:
  /*
   * \ingroup logging
   * Writes information about this action to the \p out stream.
   */
  void format_debug_output(std::ostream& out) const override;

 private:
  /**
   * Add a new reaction channel.
   *
   * \param[in] p Channel to be added.
   */
  void add_reaction(CollisionBranchPtr p);

  /**
   * Add several new reaction channels at once.
   *
   * \param[in] pv list of channels to be added.
   */
  void add_reactions(CollisionBranchList pv);

  /**
   * Perform a n->1 annihilation process.
   * \throw InvalidScatterActionMulti
   */
  void annihilation();

  /// Perform a n->2 process.
  void n_to_two();

  /**
   * Calculate the probability for a 3m-to-1 reaction according to the
   * stochastic collision criterion as given in \iref{Staudenmaier:2021lrg}.
   *
   * The formula for the probablilty is not taken from a reference, but derived
   * following the same idea as specified e.g. in the paper above.
   *
   * \f[ P_{3\rightarrow 1} = \frac{\Delta t}{(\Delta^3x)^2}
   * \frac{\pi}{4E_1E_2E_3}\frac{\Gamma_{1\rightarrow3}}{\Phi_3}
   * \mathcal{A}(\sqrt{s}),\f]
   *
   * where \f$\Phi_3\f$ represents the 3-body phase space:
   * \f[\Phi_3 = \frac{1}{(2\pi)^3)}\frac{1}{16M^2}I_3.\f]
   *
   * The defintion for the integral \f$I_3\f$ is given in the documentation of
   * calculate_I3(). Degeneracy and symmetry factors are neglected in the
   * formula, since they are treated as input for the function.
   *
   * \param[in] type_out type of outgoing particle
   * \param[in] dt timestep size
   * \param[in] gcell_vol grid cell volume
   * \param[in] degen_sym_factor degeneracy factor for reaction
   *                             (including symmetry factors)
   * \return probabilty for 3-to-1 reaction
   */
  double probability_three_to_one(const ParticleType& type_out, double dt,
                                  const double gcell_vol,
                                  const int degen_sym_factor = 1) const;

  /**
   * Calculate the probability for a 3-to-2 reaction according to the
   * stochastic collision criterion as given in \iref{Staudenmaier:2021lrg}.
   *
   * \f[ P_{3 \rightarrow 2} = \frac{1}{4E_1E_2E_3} \frac{\Delta t}{(\Delta^3
   * x)^2} \frac{\tilde{\lambda}}{\Phi_38\pi s}\sigma_{2 \rightarrow 3},\f]
   *
   * where \f$\Phi_3\f$ represents the 3-body phase space:
   * \f[\Phi_3 = \frac{1}{(2\pi)^3)}\frac{1}{16M^2}I_3.\f]
   *
   * The defintion for the integral \f$I_3\f$ is given in the documentation of
   * calculate_I3(). Degeneracy and symmetry factors are neglected in the
   * formula, since they are treated as input for the function.
   *
   * \param[in] type_out1 type of outgoing particle 1
   * \param[in] type_out2 type of outgoing particle 2
   * \param[in] dt timestep size
   * \param[in] gcell_vol grid cell volume
   * \param[in] degen_sym_factor degeneracy factor for reaction
   *                             (including symmetry factors)
   * \return probabilty for 3-to-2 reaction
   */
  double probability_three_to_two(const ParticleType& type_out1,
                                  const ParticleType& type_out2, double dt,
                                  const double gcell_vol,
                                  const double degen_sym_factor = 1.0) const;
  /**
   * Calculate the probability for a 4-to-2 reaction according to the
   * stochastic collision criterion as given in \iref{Staudenmaier:2021lrg}.
   *
   * \f[ P_{4 \rightarrow 2} = \frac{1}{16E_1E_2E_3E_4}
   * \frac{\Delta t}{(\Delta^3 x)^3} \frac{\tilde{\lambda}}{\Phi_44\pi s}
   * \sigma_{2 \rightarrow 4},\f]
   *
   * where \f$\Phi_4\f$ represents the 4-body phase space. Degeneracy and
   * symmetry factors are neglected in the formula, since they are treated as
   * input for the function.
   *
   * \param[in] type_out1 type of outgoing particle 1
   * \param[in] type_out2 type of outgoing particle 2
   * \param[in] dt timestep size
   * \param[in] gcell_vol grid cell volume
   * \param[in] degen_sym_factor degeneracy factor for reaction
   *                             (including symmetry factors)
   * \return probabilty for 4-to-2 reaction
   */
  double probability_four_to_two(const ParticleType& type_out1,
                                 const ParticleType& type_out2, double dt,
                                 const double gcell_vol,
                                 const double degen_sym_factor = 1.0) const;
  /**
   * Calculate the probability for a 5-to-2 reaction according to the
   * stochastic collision criterion as given in \iref{Garcia-Montero:2021haa}.
   *
   * \f[ P_{5 \rightarrow 2} = \frac{1}{32E_1E_2E_3E_4E_5} \frac{\Delta
   * t}{(\Delta^3 x)^4} \frac{\tilde{\lambda}}{\Phi_54\pi s}\sigma_{2
   * \rightarrow 5},\f]
   *
   * where the defintion is given without the necessary symmetry and spin
   * degeneracy factors, which are input to the function and \f$\Phi_5\f$
   * represents the 5-body phase space, which is paramaterized for the relevent
   * 5 pion state here, see documentation of parametrizaton_phi5_pions().
   *
   *
   * \param[in] m_out mass of outgoing particle types (assumes equal masses)
   * \param[in] dt timestep size
   * \param[in] gcell_vol grid cell volume
   * \param[in] degen_sym_factor degeneracy factor for reaction
   *                             (including symmetry factors)
   * \return probabilty for 5-to-2 reaction
   */
  double probability_five_to_two(const double m_out, double dt,
                                 const double gcell_vol,
                                 const double degen_sym_factor = 1.0) const;

  /**
   * Calculate the parametrized 5-pion phase space. The
   * defintion for the paramterization is given by
   * \f[\Phi_5^{param.} = A(s-s_0)^5(1+\frac{s}{s_0})^{-\alpha}\f]
   * with \f$s_0 = 25 m_{\pi}^2\f$. \f$A\f$ and \f$\alpha\f$ are fitted to
   * reproduce the phase space distribution.
   *
   *
   * \param[in] man_s mandelstam s of reaction
   * \return phase space value for 5 pions
   */
  double parametrizaton_phi5_pions(const double man_s) const;

  /**
   * Determine the spin degeneracy factor (\f$D_{spin}\f$) for the N->2
   * reaction.
   *
   * \f[D_{spin} = \frac{(2J_{out1}+1)(2J_{out2}+1)}
   * {(2J_{in1}+1)(2J_{in2}+1)(2J_{in3}+1)...(2J_{inN}+1)}\f]
   *
   * \param[in] spin_factor_inc product of incoming spin degeneracy
   *                            (denominator in above expression)
   * \param[in] spin_degen_out1 degeneracy factor of outgoing particle 1
   * \param[in] spin_degen_out2 degeneracy factor of outgoing particle 2
   * \return spin degeneracy factor
   */
  double react_degen_factor(const int spin_factor_inc,
                            const int spin_degen_out1,
                            const int spin_degen_out2) const {
    return static_cast<double>(spin_degen_out1 * spin_degen_out2) /
           static_cast<double>(spin_factor_inc);
  }

  /**
   * Check wether the three incoming particles are π⁺,π⁻,π⁰ in any order.
   * Wrapper for unwieldy if statment.
   *
   * \param[in] data_a data for first incoming particle
   * \param[in] data_b data for second incoming particle
   * \param[in] data_c data for third incoming particle
   * \return true if combination of π⁺,π⁻,π⁰
   */
  bool three_different_pions(const ParticleData& data_a,
                             const ParticleData& data_b,
                             const ParticleData& data_c) const;

  /**
   * Check wether the three incoming particles are π⁺,π⁻,η or π⁰,π⁰,η in any
   * order. Wrapper for unwieldy if statment.
   *
   * \param[in] data_a data for first incoming particle
   * \param[in] data_b data for second incoming particle
   * \param[in] data_c data for third incoming particle
   * \return true if combination of π⁺,π⁻,η or π⁰,π⁰,η
   */
  bool two_pions_eta(const ParticleData& data_a, const ParticleData& data_b,
                     const ParticleData& data_c) const;

  /**
   * Check if 5 incoming particles match intial pion state for 5-to-2, which is
   * pi+ pi- pi+ pi- pi0 in order to match the NNbar resonance treatment.
   */
  bool all_incoming_particles_are_pions_have_zero_charge_only_one_piz() const;

  /// Total probability of reaction
  double total_probability_;

  /// Partial probability of the chosen outgoing channel
  double partial_probability_;

  /// Spin interaction type
  SpinInteractionType spin_interaction_type_;

  /// List of possible collisions
  CollisionBranchList reaction_channels_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_SCATTERACTIONMULTI_H_
