/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/random.h"
#include <random>

namespace smash {
/*thread_local (see #3075)*/ random::Engine random::engine;
}  // namespace smash
