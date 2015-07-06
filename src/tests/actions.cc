/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "setup.h"

#include <algorithm>

#include "../include/actions.h"
#include "../include/decayaction.h"

using namespace Smash;

TEST(construct_and_insert) {
  Test::create_smashon_particletypes();

  // use different times for different actions
  constexpr float time_1 = 1.f;
  constexpr float time_2 = 2.f;
  constexpr float time_3 = 3.f;
  constexpr float time_4 = 4.f;
  constexpr float time_5 = 5.f;
  constexpr float time_6 = 6.f;

  constexpr float current_time = 10.5f;

  // create arbitrary particle
  ParticleData testparticle = Test::smashon(Test::Momentum{0.2, 0., .1, 0.},
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
