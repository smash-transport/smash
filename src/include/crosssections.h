/*
 *
 *    Copyright (c) 2015-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_CROSSSECTIONS_H_
#define SRC_INCLUDE_CROSSSECTIONS_H_

#include "forwarddeclarations.h"
#include "isoparticletype.h"
#include "particles.h"

namespace smash {

class cross_sections {
 public:

  cross_sections(const ParticleList &scat_particles, const double sqrt_s);

  /**
   * Determine the elastic cross section for this collision. If elastic_par is
   * given (and positive), we just use a constant cross section of that size,
   * otherwise a parametrization of the elastic cross section is used
   * (if available).
   *
   * \param[in] elast_par Elastic cross section parameter from the input file.
   *
   * \return A ProcessBranch object containing the cross section and
   * final-state IDs.
   */
  CollisionBranchPtr elastic(double elast_par);

  /**
   * Determine the (parametrized) elastic cross section for a
   * nucleon-nucleon collision.
   */
  double nn_el();

  /**
   * Determine the elastic cross section for a nucleon-kaon collision.
   * It is given by a parametrization of experimental data.
   */
  double nk_el();

  /**
   * Determine the elastic cross section for a nucleon-pion collision.
   * It is given by a parametrization of experimental data.
   */
  double npi_el();

  /**
  * Find all resonances that can be produced in a 2->1 collision of the two
  * input particles and the production cross sections of these resonances.
  *
  * Given the data and type information of two colliding particles,
  * create a list of possible resonance production processes
  * and their cross sections.
  *
  * \return A list of processes with resonance in the final state.
  * Each element in the list contains the type of the final-state particle
  * and the cross section for that particular process.
  */
  CollisionBranchList two_to_one();

  /**
   * Return the 2-to-1 resonance production cross section for a given resonance.
   *
   * \param[in] type_resonance Type information for the resonance to be
   * produced.
   * \param[in] cm_momentum_sqr Square of the center-of-mass momentum of the
   * two initial particles.
   *
   * \return The cross section for the process
   * [initial particle a] + [initial particle b] -> resonance.
   */
  double formation(const ParticleType &type_resonance,
                                             double cm_momentum_sqr);


  /** Find all inelastic 2->2 processes for the given scattering.
   * This function calls the different, more specific functions for
   * the different scatterings.
   */
  CollisionBranchList two_to_two();

  /** Find all inelastic 2->2 processes for Baryon-Baryon Scattering. */
  CollisionBranchList bb_xx();

  /** Find all inelastic 2->2 processes for Nucelon-Nucelon Scattering.
   * Calculate cross sections for resonance production from
   * nucleon-nucleon collisions (i.e. N N -> N R, N N -> Delta R).
   *
   * Checks are processed in the following order:
   * 1. Charge conservation
   * 2. Isospin factors (Clebsch-Gordan)
   * 3. Enough energy for all decay channels to be available for the resonance
   *
   * \return List of resonance production processes possible in the collision
   * of the two nucleons. Each element in the list contains the type(s) of the
   * final state particle(s) and the cross section for that particular process.
   */
  CollisionBranchList nn_xx();

  /** Find all inelastic 2->2 processes for Nucelon-Kaon Scattering. */
  CollisionBranchList nk_xx();

  /** Find all inelastic 2->2 processes for Delta-Kaon Scattering. */
  CollisionBranchList deltak_xx();

  /** Find all inelastic 2->2 processes for Hyperon-Pion Scattering. */
  CollisionBranchList ypi_xx();

  /** Determine the parametrized total cross section at high energies
   * for the given collision, which is non-zero for Baryon-Baryon and
   * Nucleon-Pion scatterings currently.
   */
  CollisionBranchPtr high_energy();

  /**
  * Calculate cross sections for resonance absorption
  * (i.e. NR->NN and ΔR->NN).
  *
  * \param[in] is_anti_particles Whether the colliding particles are
  * antiparticles
  *
  * \return List of possible resonance absorption processes. Each element of the
  * list contains the types of the final-state particles and the cross section
  * for that particular process.
  */
  CollisionBranchList bar_bar_to_nuc_nuc(const bool is_anti_particles);


  /**
   * Scattering matrix amplitude squared (divided by 16π) for resonance
   * production processes like NN → NR and NN → ΔR, where R is a baryon
   * resonance (Δ, N*, Δ*). Includes no spin or isospin factors.
   *
   * \param[in] type_a Type information for the first final-state particle.
   * \param[in] type_b Type information for the second final-state particle.
   * \param[in] twoI Twice the total isospin of the involved state.
   *
   * \return Matrix amplitude squared \f$ |\mathcal{M}(\sqrt{s})|^2/16\pi \f$.
   */
   // TODO WHY was this static before
  double nn_to_resonance_matrix_element(const ParticleType &type_a,
                                               const ParticleType &type_b,
                                               const int twoI) const;

  /**
   * Utility function to avoid code replication in nn_xx().
   */
  template <class IntegrationMethod>
  CollisionBranchList find_nn_xsection_from_type(
      const ParticleTypePtrList &type_res_1,
      const ParticleTypePtrList &type_res_2,
      const IntegrationMethod integrator);

  /** Determine the momenta of the incoming particles in the
   * center-of-mass system.
   */
  double cm_momentum() const;

  /**
   * Add a 2-to-2 channel to a collision branch list given a cross section.
   *
   * The cross section is only calculated if there is enough energy
   * for the process. If the cross section is small, the branch is not added.
   */
  template <typename F>
  inline void add_channel(CollisionBranchList& process_list, F get_xsection,
                          double sqrts, const ParticleType& type_a,
                          const ParticleType& type_b) {
    const double sqrt_s_min =
        type_a.min_mass_spectral() + type_b.min_mass_spectral();
    if (sqrts <= sqrt_s_min) {
      return;
    }
    const auto xsection = get_xsection();
    if (xsection > really_small) {
      process_list.push_back(make_unique<CollisionBranch>(
          type_a, type_b, xsection, ProcessType::TwoToTwo));
    }
  }


  /**
   * Calculate the detailed balance factor R such that
   * \f[ R = \sigma(AB \to CD) / \sigma(CD \to AB) \f]
   * where $A, B, C, D$ are stable.
   */
  inline double detailed_balance_factor_stable(double s, const ParticleType& a,
                                               const ParticleType& b,
                                               const ParticleType& c,
                                               const ParticleType& d) {
    double spin_factor = (c.spin() + 1) * (d.spin() + 1);
    spin_factor /= (a.spin() + 1) * (b.spin() + 1);
    double symmetry_factor = (1 + (a == b));
    symmetry_factor /= (1 + (c == d));
    const double momentum_factor = pCM_sqr_from_s(s, c.mass(), d.mass()) /
                                   pCM_sqr_from_s(s, a.mass(), b.mass());
    return spin_factor * symmetry_factor * momentum_factor;
  }

  /**
   * Calculate the detailed balance factor R such that
   * \f[ R = \sigma(AB \to CD) / \sigma(CD \to AB) \f]
   * where $A$ is unstable, $B$ is a kaon and $C, D$ are stable.
   */
  inline double detailed_balance_factor_RK(double sqrts, double pcm,
                                           const ParticleType& a,
                                           const ParticleType& b,
                                           const ParticleType& c,
                                           const ParticleType& d) {
    assert(!a.is_stable());
    assert(b.pdgcode().is_kaon());
    double spin_factor = (c.spin() + 1) * (d.spin() + 1);
    spin_factor /= (a.spin() + 1) * (b.spin() + 1);
    double symmetry_factor = (1 + (a == b));
    symmetry_factor /= (1 + (c == d));
    const double momentum_factor =
        pCM_sqr(sqrts, c.mass(), d.mass()) /
        (pcm * a.iso_multiplet()->get_integral_RK(sqrts));
    return spin_factor * symmetry_factor * momentum_factor;
  }

  /**
   * Calculate the detailed balance factor R such that
   * \f[ R = \sigma(AB \to CD) / \sigma(CD \to AB) \f]
   * where $A$ and $B$ are unstable, and $C$ and $D$ are stable.
   */
  inline double detailed_balance_factor_RR(double sqrts, double pcm,
                                           const ParticleType& particle_a,
                                           const ParticleType& particle_b,
                                           const ParticleType& particle_c,
                                           const ParticleType& particle_d) {
    assert(!particle_a.is_stable());
    assert(!particle_b.is_stable());
    double spin_factor = (particle_c.spin() + 1) * (particle_d.spin() + 1);
    spin_factor /= (particle_a.spin() + 1) * (particle_b.spin() + 1);
    double symmetry_factor = (1 + (particle_a == particle_b));
    symmetry_factor /= (1 + (particle_c == particle_d));
    const double momentum_factor =
        pCM_sqr(sqrts, particle_c.mass(), particle_d.mass()) /
        (pcm * particle_a.iso_multiplet()->get_integral_RR(particle_b, sqrts));
    return spin_factor * symmetry_factor * momentum_factor;
  }


  // CollisionBranchList NNbar_annihilation(); // TODO
  // CollisionBranchList NNbar_creation();  // TODO
  // CollisionBranchList call_correct_xs();
  // bool decide_string();



 private:

  /** List with data of scattering particles.  */
  ParticleList scattering_particles_;

  /** total energy in the center-of-mass frame. */
  double sqrt_s_;


};

}  // namespace smash

#endif  // SRC_INCLUDE_CROSSSECTIONS_H_
