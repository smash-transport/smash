/*
 *    Copyright (c) 2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_SMASH_ALPHACLUSTEREDNUCLEUS_H_
#define SRC_INCLUDE_SMASH_ALPHACLUSTEREDNUCLEUS_H_

#include <map>
#include <vector>

#include "angles.h"
#include "configuration.h"
#include "forwarddeclarations.h"
#include "input_keys.h"
#include "nucleus.h"
#include "threevector.h"

namespace smash {

/**
 * Child of \c Nucleus for alpha clustered nuclei.
 *
 * All options from the nucleus will still apply. The alpha-clustered nucleus
 * adds new or updated features which are outlined below.
 */
class AlphaClusteredNucleus : public Nucleus {
 public:
  /**
   * Constructor which takes a configuration and the number of test particles.
   *
   * \param[in] config the input configuration object
   * \param[in] n_test number of test particles
   * \param[in] automatic whether or not the alpha-clustering parameters should
   *                      be set automatically
   */
  AlphaClusteredNucleus(Configuration &config, int n_test, bool automatic);

  /**
   * Alpha-clustering sampling routine. This routine samples the nucleon using
   * the Woods-Saxon routine with parameters for Helium and shifts them towards
   * one of the vertices of the tetrahedron
   *
   * \return Resulting spatial position of the nucleon
   */
  ThreeVector distribute_nucleon() override;

  /**
   * Scales the tetrahedron vertex positions to have the specified side length.
   *
   * \param[in] side_length Side length to be used to scale the tetrahedron
   */
  void scale_tetrahedron_vertex_positions(double side_length);

 private:
  /**
   * Side length of the tetrahedron used for alpha-clustering.
   *
   * \note The default is set to the default of the projectile key, although
   * this class is instantiated for the target, too. However, this should be
   * fine since the keys defaults are forced to be the same in the database.
   */
  double tetrahedron_side_length_ =
      InputKeys::modi_collider_projectile_alphaClustered_sideLength
          .default_value();
  /**
   * Positions of the vertices of the regular tetrahedron with center at (0,0,0)
   * used for alpha-clustering.
   */
  std::vector<ThreeVector> tetrahedron_vertex_positions_ = {
      {1, 0.0, 0.0},
      {-1.0 / 3, std::sqrt(8) / 3, 0.0},
      {-1.0 / 3, -std::sqrt(8) / 6, std::sqrt(24) / 6},
      {-1.0 / 3, -std::sqrt(8) / 6, -std::sqrt(24) / 6}};

  /**
   * An index to iterate through the vertices of the tetrahedron. Used in the
   * alpha-clustering sampling routine to shift every 4th nucleon to the same
   * vertex.
   */
  int tetrahedron_vertex_index_ = 0;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_ALPHACLUSTEREDNUCLEUS_H_
