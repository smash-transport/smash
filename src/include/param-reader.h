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
class Sphere;

extern const char *sep;

/* read params file parameters */
void process_laboratory_config(Laboratory *parameters, char *path);
void process_box_config(Box *cube, char *path);
void process_sphere_config(Sphere *ball, char *path);

#endif  // SRC_INCLUDE_PARAM_READER_H_
