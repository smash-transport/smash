/*
 *    Copyright (c) 2013
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_INPUT_CROSS_SECTIONS_H_
#define SRC_INCLUDE_INPUT_CROSS_SECTIONS_H_

namespace Smash {

/* forward declarations */
class Particles;

/* read cross section parametrizations */
void input_cross_sections(CrossSections *cross_sections, char *path);

}  // namespace Smash

#endif  // SRC_INCLUDE_INPUT_CROSS_SECTIONS_H_
