/*
 *
 *    Copyright (c) 2020-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_CROSSSECTIONSBREMS_H_
#define SRC_INCLUDE_CROSSSECTIONSBREMS_H_

#include <initializer_list>

#include "interpolation.h"

namespace smash {

/// Interpolation object containing Xsections for π + π -> π + π + γ
static std::unique_ptr<InterpolateDataLinear<double>> pipi_interpolation =
    nullptr;

/// Center-of-mass energy.
const std::initializer_list<double> BREMS_PIPI_SQRTS = {0.0, 0.5, 1.0, 1.5,
                                                        2.0, 2.5, 3.0};
/// π + π -> π + π + γ cross section
const std::initializer_list<double> BREMS_PIPI_SIG = {0.5, 0.8, 1.1, 1.7,
                                                      2.6, 2.8, 3.6};

/// Interpolation object containing Xsections for π0 + π -> π0 + π + γ
static std::unique_ptr<InterpolateDataLinear<double>> pi0pi_interpolation =
    nullptr;

/// Center-of-mass energy.
const std::initializer_list<double> BREMS_PI0PI_SQRTS = {0.0, 0.5, 1.0, 1.5,
                                                         2.0, 2.5, 3.0};
/// π0 + π -> π0 + π + γ cross section
const std::initializer_list<double> BREMS_PI0PI_SIG = {0.5, 0.6, 1.4, 1.7,
                                                       2.2, 2.8, 3.6};

}  // namespace smash

#endif  // SRC_INCLUDE_CROSSSECTIONSBREMS_H_
