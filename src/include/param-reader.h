/*
 *    Copyright (c) 2012
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_PARAM_READER_H_
#define SRC_INCLUDE_PARAM_READER_H_

#include <list>

#include "include/parameters.h"

extern const char *sep;

/* read params file parameters */
void process_config(std::list<Parameters> *configuration, char *path);

#endif  // SRC_INCLUDE_PARAM_READER_H_
