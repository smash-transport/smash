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

#include <iosfwd>

// the forward declarations should not appear in doxygen output
#ifndef DOXYGEN

#ifdef _LIBCPP_BEGIN_NAMESPACE_STD
_LIBCPP_BEGIN_NAMESPACE_STD
#else
namespace std {
#endif

template <typename T>
class allocator;
template <typename T, typename A>
class vector;

template <typename T>
struct default_delete;
template <typename T, typename Deleter>
class unique_ptr;

#ifdef _LIBCPP_END_NAMESPACE_STD
_LIBCPP_END_NAMESPACE_STD
#else
}  // namespace std
#endif

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

class Action;
class ScatterAction;
class BoxModus;
class Clock;
class Configuration;
class CrossSections;
class DecayModes;
class FourVector;
class ThreeVector;
class ModusDefault;
class OutputInterface;
class ParticleData;
class Particles;
class ParticleType;
class ParticleTypePtr;
class PdgCode;
class DecayBranch;
class CollisionBranch;
struct ExperimentParameters;

using ActionPtr = build_unique_ptr_<Action>;
using ScatterActionPtr = build_unique_ptr_<ScatterAction>;
using ActionList = build_vector_<ActionPtr>;

using OutputsList = build_vector_<build_unique_ptr_<OutputInterface>>;
using ParticleList = build_vector_<ParticleData>;

using ParticleTypeList = build_vector_<ParticleType>;
using ParticleTypePtrList = build_vector_<ParticleTypePtr>;

template<typename T>
using ProcessBranchPtr = build_unique_ptr_<T>;
template<typename T>
using ProcessBranchList = build_vector_<ProcessBranchPtr<T>>;
using DecayBranchPtr = build_unique_ptr_<DecayBranch>;
using DecayBranchList = build_vector_<DecayBranchPtr>;
using CollisionBranchPtr = build_unique_ptr_<CollisionBranch>;
using CollisionBranchList = build_vector_<CollisionBranchPtr>;

namespace bf = boost::filesystem;

}  // namespace Smash

#endif  // DOXYGEN
#endif  // SRC_INCLUDE_FORWARDDECLARATIONS_H_
