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
 * AlphaClusteredNucleus: Child of nucleus for alpha clustered nuclei.
 *
 * All options from the nucleus will still apply. The alpha-clustered nucleus
 * adds new or updated features which are outlined below.
 */
class AlphaClusteredNucleus : public Nucleus {
 public:
  /**
   * Constructor for AlphaClusteredNuclei which takes a particle list and the
   * number of testparticles. This constructor is only used for testing
   * purposes. \param[in] particle_list Map with PDGCode and number of particles
   * which make up the nucleus \param[in] nTest number of testparticles
   */
  AlphaClusteredNucleus(const std::map<PdgCode, int> &particle_list, int nTest);
  /**
   * Constructor for AlphaClusteredNucleus, that needs the configuration
   * parameters from the inputfile and the number of testparticles \param[in]
   * config contains the parameters from the inputfile on the numbers of
   * particles with a certain PDG code \param[in] nTest number of testparticles
   * \param[in] auto_alphaclustering whether or not the alpha-clustering
   * parameters should be set automatically
   */
  AlphaClusteredNucleus(Configuration &config, int nTest,
                        bool auto_alphaclustering);
  /**
   * Alpha-Clustering sampling routine. This routine samples the nucleon using
   * the woods saxon routine with parameters for Helium and shifts them towards
   * one of the vertices of the tetrahedron
   *
   * \return Spatial position of nucleon from sampling the woods saxon
   * distribution and a shift towards a vertex of the tetrahedron
   */
  ThreeVector distribute_nucleon() override;
  /**
   * Scales the tetrahedron vertex positions to have the sidelength that is
   * specified in the config file. By default, this scales the vertex positions
   * so that the tetrahedron has sidelength 3.42fm. This value was taken from
   * \iref{Li:2020vrg}. \param[in] sidelength Sidelength for the tetrahedron
   * that was specified in the config file.
   */
  void scale_tetrahedron_vertex_positions(double sidelength);

 private:
  /// Sidelength of the tetrahedron used for alpha-clustering.
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
