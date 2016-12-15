/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SCATTERACTIONSFINDERPHOTON_H_
#define SRC_INCLUDE_SCATTERACTIONSFINDERPHOTON_H_

#include <vector>

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
                       bool two_to_one, bool two_to_two, double low_snn_cut,
                       bool strings_switch, const std::vector<bool> &nucleus_id,
                       int N_tot, int N_proj, int nofp)
      : ScatterActionsFinder(config, parameters, two_to_one, two_to_two,
                             low_snn_cut, strings_switch, nucleus_id,
                             N_tot, N_proj),
        number_of_fractional_photons(nofp) {}

  /** Constructor for testing purposes. */
  ScatterActionsFinderPhoton(float elastic_parameter, int testparticles,
                              const std::vector<bool> &nucleus_id)
      : ScatterActionsFinder(elastic_parameter, testparticles, nucleus_id) {}

  /// Number of fractional photons produced per single reaction
  int number_of_fractional_photons;

 private:
  /* Construct only ScatterActionPhoton objects. */
  ScatterActionPtr construct_scatter_action(
      const ParticleData &data_a, const ParticleData &data_b,
      float time_until_collision) const override;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTIONSFINDERPHOTON_H_
