/*
 *
 *    Copyright (c) 2014-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/random.h"
#include <random>

namespace smash {
/*thread_local (see #3075)*/ Random::Engine Random::engine;
}  // namespace smash
