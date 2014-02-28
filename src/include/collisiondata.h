/*
 *    Copyright (c) 2012-2013
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_COLLISIONDATA_H_
#define SRC_INCLUDE_COLLISIONDATA_H_

#include <cmath>
#include <vector>

#include "include/fourvector.h"

namespace Smash {

class CollisionData {
 public:
  /* Use improbable values for default constructor */
  CollisionData() : process_type_(-1), collision_time_(0.0) {}
  int process_type(void) const;
  double collision_time(void) const;
  void inline set_collision_time(const double &collision_time);
  void inline set_collision(int id_a, int id_b, int collision_type,
                            const double &collision_time);
  void inline set_collision_past(void);
  int id_partner(int i) const;

 private:
  /* Type of interaction. 0: 2->2, 2: 1->2 >99 ( = PDG code): 2->1 */
  int process_type_;
  /* collision time */
  double collision_time_;
  /* Particles id's we'd collide against */
  std::vector<int> id_partner_;
};

/* look up the process type */
int inline CollisionData::process_type(void) const {
  return process_type_;
}

/* look up the collision time */
double inline CollisionData::collision_time(void) const {
  return collision_time_;
}

/* set the collision time */
void inline CollisionData::set_collision_time(const double &collision_t) {
  collision_time_ = collision_t;
}

/* look up the id of the collision partner */
int inline CollisionData::id_partner(int i) const {
  return id_partner_[i];
}

/* set possible collision data */
void inline CollisionData::set_collision(int id_a, int id_b,
  int proc_type, const double &collision_t) {
  id_partner_.push_back(id_a);
  id_partner_.push_back(id_b);
  process_type_ = proc_type;
  collision_time_ = collision_t;
}

/* set happened collision data */
void inline CollisionData::set_collision_past() {
  collision_time_ = 0.0;
  id_partner_.clear();
  process_type_ = -1;
}

}  // namespace Smash

#endif  // SRC_INCLUDE_COLLISIONDATA_H_
