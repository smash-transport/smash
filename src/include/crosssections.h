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
   * Determine the elastic cross section for a nucleon-kaon collision.
   * It is given by a parametrization of experimental data.
   */
  double nk_el();
  /**
   * Determine the (parametrized) elastic cross section for a
   * nucleon-nucleon collision.
   */
  double nn_el();
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
  CollisionBranchList two_to_one();  // exclude NN

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


  // CollisionBranchList call_correct_xs();
  //
  // CollisionBranchList two_to_two();
  //
  // CollisionBranchList bb_xx();
  // CollisionBranchList nn_xx();
  // CollisionBranchList nk_xx();
  // CollisionBranchList ypi_xx();
  // CollisionBranchList deltak_xx();
  // /** NNbar annihilation thru NNbar → ρh₁(1170); combined with the decays
  //  *  ρ → ππ and h₁(1170) → πρ, this gives a final state of 5 pions.
  //  *  Only use in cases when detailed balance MUST happen, i.e. in a box! */
  // CollisionBranchList NNbar_annihilation(); // TODO NNbar treatment flag
  // CollisionBranchList NNbar_creation();  // TODO implmntation will not work anymore
  //
  //





  // CollisionBranchPtr high_energy();
  //
  // CollisionBranchList bb_he(); //TODO
  // CollisionBranchList bm_he(); //TODO
  //
  //
  //
  // bool decide_string();


 private:

  // helper functions
  // CollisionBranchList bar_bar_to_nuc_nuc(const bool is_anti_particles);
  // double nn_to_resonance_matrix_element(const int twoI);

  /** List with data of scattering particles.  */
  ParticleList scattering_particles_;

  /** total energy in the center-of-mass frame. */
  double sqrt_s_;


};

}  // namespace smash

#endif  // SRC_INCLUDE_CROSSSECTIONS_H_
