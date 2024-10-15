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
   * Constructor which takes a particle list and the number of test particles.
   * This constructor is only used for testing purposes.
   *
   * \param[in] particle_list Map with PDG code and number of particles which
   *                          make up the nucleus
   * \param[in] nTest number of test particles
   */
  AlphaClusteredNucleus(const std::map<PdgCode, int> &particle_list, int nTest);

  /**
   * Constructor which takes a configuration and the number of test particles.
   *
   * \param[in] config the input configuration object
   * \param[in] nTest number of test particles
   * \param[in] auto_alphaclustering whether or not the alpha-clustering
   *                                 parameters should be set automatically
   */
  AlphaClusteredNucleus(Configuration &config, int nTest,
                        bool auto_alphaclustering);

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
   * \param[in] sidelength Side length to be used to scale the tetrahedron
   */
  void scale_tetrahedron_vertex_positions(double sidelength);

 private:
  /// Side length of the tetrahedron used for alpha-clustering.
  double tetrahedron_sidelength_ = 3.42;
  /**
   * Positions of the vertices of the regular tetrahedron with center at (0,0,0)
   * used for alpha-clustering.
   */
  std::vector<ThreeVector> tetrahedron_vertex_positions_ = {
      {1, 0.0, 0.0},
      {-1.0 / 3, sqrt(8) / 3, 0.0},
      {-1.0 / 3, -sqrt(8) / 6, sqrt(24) / 6},
      {-1.0 / 3, -sqrt(8) / 6, -sqrt(24) / 6}};

  /**
   * An index to iterate through the vertices of the tetrahedron. Used in the
   * alpha-clustering sampling routine to shift every 4th nucleon to the same
   * vertex.
   */
  int tetrahedron_vertex_index_ = 0;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_ALPHACLUSTEREDNUCLEUS_H_
