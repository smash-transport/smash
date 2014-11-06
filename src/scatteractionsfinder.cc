/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteractionsfinder.h"

#include "include/action.h"
#include "include/configuration.h"
#include "include/constants.h"
#include "include/cxx14compat.h"
#include "include/experimentparameters.h"
#include "include/logging.h"
#include "include/macros.h"
#include "include/particles.h"
#include "include/resonances.h"

namespace Smash {
/*!\Userguide
* \page input_collision_term_ Collision_Term
* \key Sigma (float, optional, default = 0.0 [mb]) \n
* Elastic cross section parameter
*/

ScatterActionsFinder::ScatterActionsFinder(
    Configuration config, const ExperimentParameters &parameters)
    : ActionFinderInterface(parameters.timestep_duration()) {
/*read in parameter for elastic cross section */ 
  if (config.has_value({"Collision_Term", "Sigma"})) {
	elastic_parameter_ =  config.take({"Collision_Term", "Sigma"});
  } 
}
  
double ScatterActionsFinder::collision_time(const ParticleData &p1,
                                            const ParticleData &p2) {
  const auto &log = logger<LogArea::FindScatter>();
  /* UrQMD collision time
   * arXiv:1203.4418 (5.15): in computational frame
   * position of particle a: x_a
   * position of particle b: x_b
   * momentum of particle a: p_a
   * momentum of particle b: p_b
   * t_{coll} = - (x_a - x_b) . (p_a - p_b) / (p_a - p_b)^2
   */
  ThreeVector pos_diff = p1.position().threevec() - p2.position().threevec();
  ThreeVector velo_diff = p1.velocity() - p2.velocity();
  log.trace(source_location, "\n"
            "Scatter ", p1, "\n"
            "    <-> ", p2, "\n"
            "=> position difference: ", pos_diff, " [fm]",
            ", velocity difference: ", velo_diff, " [GeV]");
  /* Zero momentum leads to infite distance, particles are not approaching. */
  if (fabs(velo_diff.sqr()) < really_small) {
    return -1.0;
  } else {
    return -pos_diff * velo_diff / velo_diff.sqr();
  }
}

ActionPtr ScatterActionsFinder::check_collision(
    const ParticleData &data_a, const ParticleData &data_b) const {
  const auto &log = logger<LogArea::FindScatter>();

  /* just collided with this particle */
  if (data_a.id_process() >= 0 && data_a.id_process() == data_b.id_process()) {
    log.debug("Skipping collided particles at time ", data_a.position().x0(),
              " due to process ", data_a.id_process(),
              "\n    ", data_a,
              "\n<-> ", data_b);
    return nullptr;
  }

  /* check according timestep: positive and smaller */
  const float time_until_collision = collision_time(data_a, data_b);
  if (time_until_collision < 0.f || time_until_collision >= dt_) {
    return nullptr;
  }

  /* Create ScatterAction object. */
  std::unique_ptr<ScatterAction> act;
  if (data_a.is_baryon() && data_b.is_baryon()) {
    act = make_unique<ScatterActionBaryonBaryon>(data_a, data_b,
                                                 time_until_collision);
  } else if (data_a.is_baryon() || data_b.is_baryon()) {
    act = make_unique<ScatterActionBaryonMeson>(data_a, data_b,
                                                time_until_collision);
  } else {
    act = make_unique<ScatterActionMesonMeson>(data_a, data_b,
                                               time_until_collision);
  }

  /* Add various subprocesses.  */
  /* (1) elastic */
  act->add_process(act->elastic_cross_section(elastic_parameter_));
  /* (2) resonance formation (2->1) */
  act->add_processes(act->resonance_cross_sections());
  /* (3) 2->2 (inelastic) */
  act->add_processes(act->two_to_two_cross_sections());

  {
    /* distance criteria according to cross_section */
    const double distance_squared = act->particle_distance();
    if (distance_squared >= act->weight() * fm2_mb * M_1_PI) {
      return nullptr;
    }
    log.debug("particle distance squared: ", distance_squared,
              "\n    ", data_a,
              "\n<-> ", data_b);
  }

  return std::move(act);
}

template <typename F>
inline void iterate_all_pairs(
    const ParticleList &search_list,
    const std::vector<const ParticleList *> &neighbors_list, F &&closure) {
  const auto end0 = search_list.end();
  const auto end1 = neighbors_list.end();
  for (auto it0 = search_list.begin(); it0 != end0; ++it0) {
    for (auto it1 = std::next(it0); it1 != end0; ++it1) {
      if (it0->id() < it1->id()) {
        closure(*it0, *it1);
      } else {
        closure(*it1, *it0);
      }
    }
    for (auto it1 = neighbors_list.begin(); it1 != end1; ++it1) {
      const ParticleList &inner_neighbors_list = **it1;
      const auto end2 = inner_neighbors_list.end();
      for (auto it2 = inner_neighbors_list.begin(); it2 != end2; ++it2) {
        if (it0->id() < it2->id()) {
          closure(*it0, *it2);
        } else {
          closure(*it2, *it0);
        }
      }
    }
  }
}

std::vector<ActionPtr> ScatterActionsFinder::find_possible_actions(
    const ParticleList &search_list,
    const std::vector<const ParticleList *> &neighbors_list) const {
  std::vector<ActionPtr> actions;

  iterate_all_pairs(search_list, neighbors_list,
                    [&](const ParticleData &p1, const ParticleData &p2) {
    /* Check if collision is possible. */
    ActionPtr act = check_collision(p1, p2);

    /* Add to collision list. */
    if (act != nullptr) {
      actions.push_back(std::move(act));
    }
  });
  return std::move(actions);
}

}  // namespace Smash
