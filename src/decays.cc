/*
 *
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */

#include "include/decays.h"

#include "include/Parameters.h"
#include "include/ParticleData.h"
#include "include/ParticleType.h"
#include "include/constants.h"

/* does_decay - does a resonance decay on this timestep? */
bool does_decay(std::vector<ParticleData> *particle,
  std::vector<ParticleType> *particle_type, std::map<int, int> *map_type,
  std::list<int> *collision_list, const Parameters &parameters, int id_res) {
  /* local rest frame velocity */
  FourVector velocity_lrf;
  velocity_lrf.set_x0(1.0);
  velocity_lrf.set_x1((*particle)[id_res].momentum().x1()
                      / (*particle)[id_res].momentum().x0());
  velocity_lrf.set_x2((*particle)[id_res].momentum().x2()
                      / (*particle)[id_res].momentum().x0());
  velocity_lrf.set_x3((*particle)[id_res].momentum().x3()
                      / (*particle)[id_res].momentum().x0());

  /* The clock goes slower in the rest frame of the resonance */
  double inverse_gamma = velocity_lrf.Dot(velocity_lrf);
  double resonance_frame_timestep = parameters.eps() * inverse_gamma;

  /* Exponential decay. Average lifetime t_avr = 1 / width
   * t / t_avr = width * t (remember GeV-fm conversion)
   * P(decay at Delta_t) = width * Delta_t
   * P(alive after n steps) = (1 - width * Delta_t)^n
   * = (1 - width * Delta_t)^(t / Delta_t)
   * -> exp(-width * t) when Delta_t -> 0
   */
  if (drand48() < resonance_frame_timestep
      * (*particle_type)[(*map_type)[id_res]].width() / hbarc) {
    /* Time is up! Set the particle to decay at this timestep */
    (*particle)[id_res].set_collision(2, 0.0, -1);
    collision_list->push_back(id_res);
    return true;
  } 
  return false;
}
