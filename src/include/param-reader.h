/*
 *    Copyright (c) 2012 maximilian attems <attems@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_PARAM_READER_H_
#define SRC_INCLUDE_PARAM_READER_H_

#include <list>

#include "../include/Parameters.h"

/* forward declaration */
class Box;
class Laboratory;

extern const char *sep;

/* read params file parameters */
void process_params(char *paramfile, std::list<Parameters> *configuration);
void assign_params(std::list<Parameters> *configuration, Laboratory *lab);
void assign_params(std::list<Parameters> *configuration, Box *box);

#endif  // SRC_INCLUDE_PARAM_READER_H_
