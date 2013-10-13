/*
 *
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <list>
#include <map>

#include "include/decays.h"

#include "include/Laboratory.h"
#include "include/ParticleData.h"
#include "include/ParticleType.h"
#include "include/Particles.h"
#include "include/constants.h"
#include "include/FourVector.h"
#include "include/outputroutines.h"
#include "include/resonances.h"

/* check_decays - does a resonance decay on this timestep? */
void check_decays(Particles *particles, std::list<int> *decay_list,
  const Laboratory &parameters) {
  FourVector velocity_lrf;
  velocity_lrf.set_x0(1.0);

  for (std::map<int, ParticleData>::iterator i = particles->begin();
       i != particles->end(); ++i) {
    /* particle doesn't decay */
    if (particles->type(i->first).width() < 0.0)
      continue;
    /* local rest frame velocity */
    velocity_lrf.set_x1(i->second.momentum().x1() / i->second.momentum().x0());
    velocity_lrf.set_x2(i->second.momentum().x2() / i->second.momentum().x0());
    velocity_lrf.set_x3(i->second.momentum().x3() / i->second.momentum().x0());

    /* The clock goes slower in the rest frame of the resonance */
    double inverse_gamma = sqrt(velocity_lrf.Dot(velocity_lrf));
    double resonance_frame_timestep = parameters.eps() * inverse_gamma;

    /* Exponential decay. Average lifetime t_avr = 1 / width
     * t / t_avr = width * t (remember GeV-fm conversion)
     * P(decay at Delta_t) = width * Delta_t
     * P(alive after n steps) = (1 - width * Delta_t)^n
     * = (1 - width * Delta_t)^(t / Delta_t)
     * -> exp(-width * t) when Delta_t -> 0
     */
    if (drand48() < resonance_frame_timestep
                    * particles->type(i->first).width() / hbarc) {
      /* Time is up! Set the particle to decay at this timestep */
      i->second.set_collision(2, 0.0, -1);
      decay_list->push_back(i->first);
    }
  }
}



/* 1->2 and 1->3 decay processes */
size_t decay_particles(Particles *particles, std::list<int> *decay_list,
  size_t id_process) {
  FourVector velocity_CM;

  for (std::list<int>::iterator id = decay_list->begin();
    id != decay_list->end(); ++id) {
    /* relevant particle id's for the collision */
    int id_a = *id;
    int interaction_type = particles->data(id_a).process_type();

    if (interaction_type != 2)
      printf("Decays warning: ID %i (%s) has process type %i.\n",
           id_a, particles->type(id_a).name().c_str(), interaction_type);

    FourVector initial_momentum(particles->data(id_a).momentum());

    /* 1->2 resonance decay */
    printd("Process: Resonance decay. ");
    write_oscar(particles->data(id_a), particles->type(id_a), 1, 2);
    printd_momenta("Resonance momenta before decay", particles->data(id_a));

    /* boost to rest frame */
    velocity_CM.set_x0(1.0);
    velocity_CM.set_x1(particles->data(id_a).momentum().x1()
                       / particles->data(id_a).momentum().x0());
    velocity_CM.set_x2(particles->data(id_a).momentum().x2()
                       / particles->data(id_a).momentum().x0());
    velocity_CM.set_x3(particles->data(id_a).momentum().x3()
                       / particles->data(id_a).momentum().x0());
    particles->data_pointer(id_a)->set_momentum(
        particles->data(id_a).momentum().LorentzBoost(velocity_CM));
    particles->data_pointer(id_a)->set_position(
        particles->data(id_a).position().LorentzBoost(velocity_CM));

    printd_momenta("Boosted resonance momenta before decay",
                   particles->data(id_a));

    size_t old_max_id = particles->id_max();
    size_t id_new_a = resonance_decay(particles, id_a);
    size_t id_new_b = id_new_a + 1;

    printd_momenta("particle 1 momenta in lrf", particles->data(id_new_a));
    printd_momenta("particle 2 momenta in lrf", particles->data(id_new_b));

    boost_back_CM(particles->data_pointer(id_new_a),
                  particles->data_pointer(id_new_b), &velocity_CM);

    write_oscar(particles->data(id_new_a), particles->type(id_new_a));
    write_oscar(particles->data(id_new_b), particles->type(id_new_b));

    size_t new_particles = particles->id_max() - old_max_id;
    int id_new_c = -1;
    if (new_particles == 3) {
      id_new_c = id_new_b + 1;
      FourVector velocity = velocity_CM;
      velocity *= -1;
      velocity.set_x0(1.0);
      FourVector momentum_c(particles->data(id_new_c).momentum());
      FourVector position_c(particles->data(id_new_c).position());
      /* Boost the momenta back to lab frame */
      momentum_c = momentum_c.LorentzBoost(velocity);
      /* Boost the positions back to lab frame */
      position_c = position_c.LorentzBoost(velocity);

      particles->data_pointer(id_new_c)->set_momentum(momentum_c);
      particles->data_pointer(id_new_c)->set_position(position_c);
      write_oscar(particles->data(id_new_c), particles->type(id_new_c));
    }

    printd_momenta("particle 1 momenta in comp", particles->data(id_new_a));
    printd_momenta("particle 2 momenta in comp", particles->data(id_new_b));

    FourVector final_momentum(particles->data(id_new_a).momentum()
      + particles->data(id_new_b).momentum());

    /* unset collision time for both particles + keep id + unset partner */
    particles->data_pointer(id_new_a)->set_collision_past(id_process);
    particles->data_pointer(id_new_b)->set_collision_past(id_process);
    if (new_particles == 3) {
      final_momentum += particles->data(id_new_c).momentum();
      particles->data_pointer(id_new_c)->set_collision_past(id_process);
    }
    printd("Particle map has now %zu elements. \n", particles->size());

    id_process++;

    FourVector momentum_difference;
    momentum_difference += initial_momentum;
    momentum_difference -= final_momentum;
    if (fabs(momentum_difference.x0()) > really_small) {
      printf("Process %zu type %i particle %s decay to %s and %s ",
        id_process, interaction_type, particles->type(id_a).name().c_str(),
             particles->type(id_new_a).name().c_str(),
             particles->type(id_new_b).name().c_str());
      if (new_particles == 3) {
        printf("and %s ", particles->type(id_new_c).name().c_str());
      }
      printf("time %g\n", particles->data(id_a).position().x0());
      printf("Warning: Interaction type %i E conservation violation %g\n",
             interaction_type, momentum_difference.x0());
    }
    if (fabs(momentum_difference.x1()) > really_small)
      printf("Warning: Interaction type %i px conservation violation %g\n",
             interaction_type, momentum_difference.x1());
    if (fabs(momentum_difference.x2()) > really_small)
      printf("Warning: Interaction type %i py conservation violation %g\n",
             interaction_type, momentum_difference.x2());
    if (fabs(momentum_difference.x3()) > really_small)
      printf("Warning: Interaction type %i pz conservation violation %g\n",
             interaction_type, momentum_difference.x3());

    particles->remove(id_a);
    printd("ID %i has decayed and removed from the list.\n", id_a);
  } /* end for std::list<int>::iterator id = decay_list->begin() */
  /* empty the decay table */
  decay_list->clear();
  printd("Decay list done.\n");

  /* return how many processes we have handled so far*/
  return id_process;
}
