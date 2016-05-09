/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SCATTERACTIONFINDERPHOTON_H_
#define SRC_INCLUDE_SCATTERACTIONFINDERPHOTON_H_

#include "scatteractionsfinder.h"

namespace Smash {

/**
 * The Scatter Action Finder for Photon producing reactions
*/

class ScatterActionsFinderPhoton : public ScatterActionsFinder {
 public:
    /** Initialize the finder with the given parameters. */
  ScatterActionsFinderPhoton(Configuration config,
                       const ExperimentParameters &parameters,
                       bool two_to_one, bool two_to_two, int nofp):
    ScatterActionsFinder(config, parameters, two_to_one, two_to_two),
    number_of_fractional_photons(nofp) {}
  /** Constructor for testing purposes. */
  ScatterActionsFinderPhoton(float elastic_parameter, int testparticles):
    ScatterActionsFinder(elastic_parameter, testparticles){}

  /* At the moment no final scattering actions are implemented in SMASH */
  // ActionList find_final_actions(const Particles &search_list) const override;

  /// Number of fractional photons produced per single reaction
  int number_of_fractional_photons = 100;

 private:
  /* Construct only ScatterActionPhoton objects. */
  ScatterActionPtr construct_scatter_action(
      const ParticleData &data_a, const ParticleData &data_b,
      float time_until_collision) const override;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTIONFINDERPHOTON_H_
