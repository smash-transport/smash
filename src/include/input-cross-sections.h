/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_INPUT_CROSS_SECTIONS_H_
#define SRC_INCLUDE_INPUT_CROSS_SECTIONS_H_

/* forward declarations */
class Particles;

/* read cross section parametrizations */
void input_cross_sections(CrossSections *cross_sections, char *path);

#endif  // SRC_INCLUDE_INPUT_CROSS_SECTIONS_H_
