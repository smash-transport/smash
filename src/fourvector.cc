/*
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#include <cmath>

#include "include/FourVector.h"

/* LorentzBoost - This is equivalent to -velocity Boost from ROOT */
FourVector FourVector::LorentzBoost(FourVector origin, FourVector velocity) {
  double velocity_squared = velocity.Dot(velocity);
  /* Lorentz gamma. 4-velocity squared = 1-v^2 */
  double gamma = velocity_squared > 0 ? 1.0 / std::sqrt(velocity_squared) : 0;
  FourVector boosted_vector;

  boosted_vector.x0_ = gamma * origin.Dot(velocity);
  boosted_vector.x1_ = origin.x1() - gamma * velocity.x1()
    * (gamma / (gamma + 1) * origin.Dot(velocity) + origin.x0() / (gamma + 1));
  boosted_vector.x2_ = origin.x2() - gamma * velocity.x2()
    * (gamma / (gamma + 1) * origin.Dot(velocity) + origin.x0() / (gamma + 1));
  boosted_vector.x3_ = origin.x3() - gamma * velocity.x3()
    * (gamma / (gamma + 1) * origin.Dot(velocity) + origin.x0() / (gamma + 1));
  return boosted_vector;
}
