/*
 *
 *    Copyright (c) 2015-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "setup.h"

#include <algorithm>

#include "../include/smash/actions.h"
#include "../include/smash/decayaction.h"

using namespace smash;

TEST(construct_and_insert) {
  Test::create_smashon_particletypes();

  // use different times for different actions
  constexpr double time_1 = 1.;
  constexpr double time_2 = 2.;
  constexpr double time_3 = 3.;
  constexpr double time_4 = 4.;
  constexpr double time_5 = 5.;
  constexpr double time_6 = 6.;

  constexpr double current_time = 10.5;

  // create arbitrary particle
  ParticleData testparticle =
      Test::smashon(Test::Momentum{0.2, 0., .1, 0.},
                    Test::Position{current_time, 1., .9, 1.});

  // add actions to list
  ActionList action_vec;
  action_vec.push_back(make_unique<DecayAction>(testparticle, time_4));
  action_vec.push_back(make_unique<DecayAction>(testparticle, time_1));
  action_vec.push_back(make_unique<DecayAction>(testparticle, time_6));

  // construct the Actions object
  Actions actions(std::move(action_vec));
  VERIFY(!actions.is_empty());

  // create new actions that are then inserted into the Actions object
  ActionList new_actions;
  new_actions.push_back(make_unique<DecayAction>(testparticle, time_5));
  new_actions.push_back(make_unique<DecayAction>(testparticle, time_2));
  new_actions.push_back(make_unique<DecayAction>(testparticle, time_3));

  // insert actions
  actions.insert(std::move(new_actions));

  // verify that the actions are in the right order
  COMPARE(actions.pop()->time_of_execution(), current_time + time_1);
  COMPARE(actions.pop()->time_of_execution(), current_time + time_2);
  COMPARE(actions.pop()->time_of_execution(), current_time + time_3);
  COMPARE(actions.pop()->time_of_execution(), current_time + time_4);
  COMPARE(actions.pop()->time_of_execution(), current_time + time_5);
  COMPARE(actions.pop()->time_of_execution(), current_time + time_6);

  VERIFY(actions.is_empty());
}
