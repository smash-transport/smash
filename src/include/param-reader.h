/*
 *    Copyright (c) 2012 maximilian attems <attems@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_PARAM_READER_H_
#define SRC_INCLUDE_PARAM_READER_H_

#include <vector>

/* forward declaration */
class Box;
class Laboratory;
class Parameters;

extern const char *sep;

/* read params file parameters */
void process_params(char *paramfile, std::vector<Parameters> *configuration);
void assign_params(std::vector<Parameters> *configuration,
  Laboratory *lab);
void assign_params(std::vector<Parameters> *configuration, Box *box);

#endif  // SRC_INCLUDE_PARAM_READER_H_
