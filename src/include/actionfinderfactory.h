/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_ACTIONFINDERFACTORY_H_
#define SRC_INCLUDE_ACTIONFINDERFACTORY_H_

namespace Smash {

class ActionFinderFactory {
 public:
  virtual std::vector<ActionPtr> find_possible_actions(const ParticleList &)
      const = 0;

 private:
};

}  // namespace Smash

#endif  // SRC_INCLUDE_ACTIONFINDERFACTORY_H_
