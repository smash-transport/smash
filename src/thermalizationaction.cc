/*
 *    Copyright (c) 2017-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "smash/thermalizationaction.h"

namespace smash {

ThermalizationAction::ThermalizationAction(const GrandCanThermalizer &gct,
                                           double absolute_execution_time)
    : Action(gct.particles_to_remove(), gct.particles_to_insert(),
             absolute_execution_time, ProcessType::Thermalization) {}
}  // namespace smash
