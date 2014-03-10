/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_REGION_H_
#define SRC_INCLUDE_REGION_H_

namespace Smash {

class Region {
 public:
  Region(const ParticleList &);

 private:
  const ParticleList particles_;
  std::vector<ActionPtr> all_actions_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_REGION_H_
