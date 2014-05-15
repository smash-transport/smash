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

// the forward declarations should not appear in doxygen output
#ifndef DOXYGEN

namespace std {
template <typename T>
class allocator;
template <typename T, typename A>
class vector;

template <typename T>
struct default_delete;
template <typename T, typename Deleter>
class unique_ptr;
}  // namespace std

namespace boost {
namespace filesystem {
class path;
}  // namespace filesystem
}  // namespace boost

namespace Smash {

template <typename T>
using build_unique_ptr_ = std::unique_ptr<T, std::default_delete<T>>;
template <typename T>
using build_vector_ = std::vector<T, std::allocator<T>>;

class BoxModus;
class Clock;
class Configuration;
class CrossSections;
class FourVector;
class ModusDefault;
class OutputInterface;
class ParticleData;
class Particles;
class ParticleType;
class ProcessBranch;
struct ExperimentParameters;

using OutputsList = build_vector_<build_unique_ptr_<OutputInterface>>;
using ParticleList = build_vector_<ParticleData>;

namespace bf = boost::filesystem;

}  // namespace Smash

#endif  // DOXYGEN
#endif  // SRC_INCLUDE_FORWARDDECLARATIONS_H_
