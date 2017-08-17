/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/random.h"
#include <random>

namespace Smash {
/*thread_local*/ Random::Engine Random::engine;
}  // namespace Smash
