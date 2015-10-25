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
class IsoParticleType;
class PdgCode;
class DecayBranch;
class CollisionBranch;
class Tabulation;
struct ExperimentParameters;
enum class CalculationFrame {
    CenterOfVelocity,
    CenterOfMass,
    FixedTarget,
};

/// Possible methods of impact parameter sampling.
enum class Sampling {
  /// Sample from uniform distribution.
  Uniform,
  /// Sample from areal / quadratic distribution.
  Quadratic,
  /// Sample from custom, user-defined distribution.
  Custom,
};

/** Allows to choose which kind of density to calculate.
*  The baryon density is necessary for the Skyrme potential.
*  For the symmetry potential one needs to know the isospin density.
*/
enum class DensityType {
  None = 0,
  Hadron = 1,
  Baryon = 2,
  BaryonicIsospin = 3,
  Pion = 4,
};

/** The time step mode.
 */
enum class TimeStepMode : char {
  /// Don't use time steps; propagate from action to action.
  None,
  /// Use fixed time step.
  Fixed,
  /// Use time step that adapts to the state of the system.
  Adaptive,
};

/** Initial condition for a particle in a box.
*
* If PeakedMomenta is used, all particles have the same momentum
* \f$p = 3 \cdot T\f$ with T the temperature.
*
* Else, a thermalized ensemble is generated (the momenta are sampled
* from a Maxwell-Boltzmann distribution).
*
* In either case, the positions in space are chosen randomly.
*/
enum class BoxInitialCondition {
  ThermalMomenta,
  PeakedMomenta,
};


using ActionPtr = build_unique_ptr_<Action>;
using ScatterActionPtr = build_unique_ptr_<ScatterAction>;
using ActionList = build_vector_<ActionPtr>;

using OutputsList = build_vector_<build_unique_ptr_<OutputInterface>>;
using ParticleList = build_vector_<ParticleData>;

using ParticleTypeList = build_vector_<ParticleType>;
using ParticleTypePtrList = build_vector_<ParticleTypePtr>;
using IsoParticleTypeList = build_vector_<IsoParticleType>;
using IsoParticleTypeRefList = build_vector_<const IsoParticleType&>;

template<typename T>
using ProcessBranchPtr = build_unique_ptr_<T>;
template<typename T>
using ProcessBranchList = build_vector_<ProcessBranchPtr<T>>;
using DecayBranchPtr = build_unique_ptr_<DecayBranch>;
using DecayBranchList = build_vector_<DecayBranchPtr>;
using CollisionBranchPtr = build_unique_ptr_<CollisionBranch>;
using CollisionBranchList = build_vector_<CollisionBranchPtr>;

using TabulationPtr = build_unique_ptr_<Tabulation>;

namespace bf = boost::filesystem;


}  // namespace Smash

#endif  // DOXYGEN
#endif  // SRC_INCLUDE_FORWARDDECLARATIONS_H_
