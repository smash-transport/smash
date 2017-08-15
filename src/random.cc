/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <random>
#include "include/random.h"

namespace Smash {
/*thread_local (see #3075)*/ Random::Engine Random::engine;
}  // namespace Smash
