/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_OUTPUTINTERFACE_H_
#define SRC_INCLUDE_OUTPUTINTERFACE_H_

namespace Smash {
class Particles;

class OutputInterface {
 public:
  virtual ~OutputInterface() = default;

  virtual void at_runstart() = 0;
  virtual void at_eventstart(const Particles &particles, const int evt_num) = 0;
  virtual void at_eventend(const Particles &particles, const int evt_num) = 0;
  //virtual void at_collision(const Collisions &collisions) = 0;
  virtual void at_outtime(const Particles &particles, const int evt_num, const int timestep) = 0;
  virtual void at_runend() = 0;
  virtual void at_crash() = 0;

};

}  // namespace Smash

#endif  // SRC_INCLUDE_OUTPUTINTERFACE_H_
