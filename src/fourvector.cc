/*
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#include <cmath>

#include "include/FourVector.h"

/* LorentzBoost - This is equivalent to -velocity Boost from ROOT */
FourVector FourVector::LorentzBoost(const FourVector &velocity) const {
  double velocity_squared = velocity.Dot(velocity);
  /* Lorentz gamma. 4-velocity squared = 1-v^2 */
  double gamma = velocity_squared > 0 ? 1.0 / std::sqrt(velocity_squared) : 0;

  return FourVector(gamma * this->Dot(velocity),
      this->x1() - gamma * velocity.x1()
    * (gamma / (gamma + 1) * this->Dot(velocity) + this->x0() / (gamma + 1)),
  this->x2() - gamma * velocity.x2()
    * (gamma / (gamma + 1) * this->Dot(velocity) + this->x0() / (gamma + 1)),
  this->x3() - gamma * velocity.x3()
    * (gamma / (gamma + 1) * this->Dot(velocity) + this->x0() / (gamma + 1)));
}
