#ifndef SRC_INCLUDE_SCATTERACTIONFINDERPHOTON_H_
#define SRC_INCLUDE_SCATTERACTIONFINDERPHOTON_H_

#include "scatteractionsfinder.h"

namespace Smash {

/**
 * The Scatter Action Finder for Photon producing reactions
*/

class ScatterActionsFinderPhoton : public ScatterActionsFinder {
  public:
    /** Initialize the finder */
    ScatterActionsFinderPhoton() {} 

    /** At the moment no final scattering actions are implemented in SMASH */
    // ActionList find_final_actions(const Particles &search_list) const override;
    
  private:
    /* Construct only ScatterActionPhoton objects. */
    ScatterActionPtr construct_scatter_action(const ParticleData &data_a, const ParticleData &data_b, float time_until_collision) const override;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTIONFINDERPHOTON_H_

