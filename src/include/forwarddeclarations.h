/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_FORWARDDECLARATIONS_H_
#define SRC_INCLUDE_FORWARDDECLARATIONS_H_

namespace std {
template <typename T>
class allocator;
template <typename T, typename A>
class vector;
}  // namespace std

namespace Smash {

class ParticleData;
class ParticleType;
class Particles;
class ProcessBranch;

using ParticleList = std::vector<ParticleData, std::allocator<ParticleData>>;

}  // namespace Smash

#endif  // SRC_INCLUDE_FORWARDDECLARATIONS_H_
