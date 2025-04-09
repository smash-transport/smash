/*
 *
 *    Copyright (c) 2014-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_FORWARDDECLARATIONS_H_
#define SRC_INCLUDE_SMASH_FORWARDDECLARATIONS_H_

/// @cond
// exclude most content here from documentation

#include <bitset>
#include <iosfwd>
#include <memory>
#include <vector>

namespace smash {

class Action;
class ScatterAction;
class ScatterActionMulti;
class BoxModus;
class Clock;
class CollisionBranch;
class Configuration;
class CrossSections;
class DecayBranch;
class DecayModes;
class DecayType;
class ExperimentBase;
class FourVector;
class IsoParticleType;
template <typename T>
class Key;
class ModusDefault;
class OutputInterface;
class ParticleData;
class Particles;
class ParticleType;
class ParticleTypePtr;
class PdgCode;
class ScatterActionsFinderParameters;
class Tabulation;
class ThreeVector;

struct ExperimentParameters;
struct InitialConditionParameters;
struct StringTransitionParameters;
struct Nucleoncorr;

/// @endcond

/// The calculation frame
enum class CalculationFrame {
  CenterOfVelocity,
  CenterOfMass,
  FixedTarget,
};

/// Option to use Fermi Motion
enum class FermiMotion {
  /// Don't use fermi motion.
  Off,
  /// Use fermi motion in combination with potentials.
  On,
  /// Use fermi motion without potentials.
  Frozen,
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

/// Modes of calculating the gradients
enum class DerivativesMode {
  CovariantGaussian,
  FiniteDifference,
  Off,
};

/**
 * This enum is here only to serve InputKeys class, but it is unused and
 * referring to a removed SMASH input key.
 */
enum class RestFrameDensityDerivativesMode {
  On,
  Off,
};

/**
 * Modes of calculating the field gradients: chain rule or direct. The modes
 * only make sense for the VDF potentials.
 */
enum class FieldDerivativesMode {
  ChainRule,
  Direct,
};

/// Modes of smearing
enum class SmearingMode {
  CovariantGaussian,
  Discrete,
  Triangular,
};

/// The time step mode.
enum class TimeStepMode : char {
  /// Don't use time steps; propagate from action to action.
  None,
  /// Use fixed time step.
  Fixed,
};

/**
 * Initial condition for a particle in a box.
 *
 * In all cases, the positions in space are chosen randomly.
 */
enum class BoxInitialCondition {
  /// A thermalized ensemble is generated, with momenta sampled from a
  /// Maxwell-Boltzmann distribution
  ThermalMomentaBoltzmann,
  /// A thermalized ensemble is generated, with momenta of baryons(mesons)
  /// sampled from a Fermi(Bose) distribution
  ThermalMomentaQuantum,
  /// All particles have the same momentum \f$p = 3 \cdot T\f$ with T being the
  /// temperature.
  PeakedMomenta,
};

/// Initial condition for a particle in a sphere
enum class SphereInitialCondition {
  /// A thermalized ensemble is generated, with momenta sampled from a
  /// Maxwell-Boltzmann distribution
  ThermalMomentaBoltzmann,
  /// A thermalized ensemble is generated, with momenta of baryons(mesons)
  /// sampled from a Fermi(Bose) distribution
  ThermalMomentaQuantum,
  /// Off-equilibrium distribution used in massless comparisons of SMASH to the
  /// extended universe metric. See eq. (76) in \iref{Bazow:2016oky}
  IC_ES,
  /// Off-equilibrium distribution used in massless comparisons of SMASH to the
  /// extended universe metric. See eq. (77) in \iref{Bazow:2016oky}
  IC_1M,
  /// Off-equilibrium distribution used in massless comparisons of SMASH to the
  /// extended universe metric. See eq. (78) in \iref{Bazow:2016oky}
  IC_2M,
  /// A generalization of IC_ES for the non-zero mass case; note that there is
  /// currently no analytical comparison possible with this distribution.
  IC_Massive,
};

/**
 * Defines properties of expansion for the metric (e.g. FRW)
 *
 * If anything else than NoExpansion is used, then a non-zero
 * Hubble parameter is computed and corrections are brought to the
 * propagation of all particles according to selected expanding
 * metric.
 */
enum class ExpansionMode {
  NoExpansion,
  MasslessFRW,
  MassiveFRW,
  Exponential,
};

/// Treatment of N Nbar Annihilation
enum class NNbarTreatment {
  /// No Annihilation
  NoAnnihilation,
  /// Use intermediate Resonances
  Resonances,
  /// Directly create 5 pions, use with multi-particle reactions
  TwoToFive,
  /// Use string fragmentation
  Strings,
};

/**
 * Represents thermodynamic quantities that can be printed out
 * \see_key{key_output_thermo_type_}
 */
enum class ThermodynamicQuantity : char {
  /// Density in the Eckart frame
  EckartDensity,
  /// Energy-momentum tensor in lab frame
  Tmn,
  /// Energy-momentum tensor in Landau rest frame
  TmnLandau,
  /// Velocity of the Landau rest frame
  LandauVelocity,
  /// Electric (Q), baryonic (B) and strange (S) currents
  j_QBS
};

/// Criteria used to check collisions
enum class CollisionCriterion {
  /// Geometric criterion.
  Geometric,
  /// Stochastic Criteiron.
  Stochastic,
  /// Covariant Criterion
  Covariant
};

/// Whether and when only final state particles should be printed.
enum class OutputOnlyFinal {
  /// Print only final-state particles.
  Yes,
  /// Print initial, intermediate and final-state particles.
  No,
  /// Print only final-state particles, and those only if the event is not
  /// empty.
  IfNotEmpty,
};

/// The different groups of 2 to 2 reactions that one can include
// Because std::bitset does not handle enum classes, this is a simple enum.
enum IncludedReactions {
  All = 50,
  Elastic = 0,
  NN_to_NR = 1,
  NN_to_DR = 2,
  KN_to_KN = 3,
  KN_to_KDelta = 4,
  Strangeness_exchange = 5,
  NNbar = 6,
  PiDeuteron_to_NN = 7,
  PiDeuteron_to_pidprime = 8,
  NDeuteron_to_Ndprime = 9,
};

/// Container for the 2 to 2 reactions in the code
typedef std::bitset<10> ReactionsBitSet;

/// The different groups of multi-particle reactions that one can include
// Because std::bitset does not handle enum classes, this is a simple enum.
enum IncludedMultiParticleReactions {
  Meson_3to1 = 0,
  Deuteron_3to2 = 1,
  NNbar_5to2 = 2,
  A3_Nuclei_4to2 = 3,
};

/// Container for the n to m reactions in the code
typedef std::bitset<4> MultiParticleReactionsBitSet;

/// Possible spin interaction types
enum class SpinInteractionType {
  /// All spin interactions
  On,
  /// No spin interactions
  Off,
  /// Spin flips in elastic collisions only
  Elastic,
};

/**
 * Defines the algorithm used for the forced thermalization.
 *  For the description of algorithms see \iref{Oliinychenko:2016vkg}.
 *  All of them intend to conserve the net baryon number, strangeness
 *  and electric charge, as well as energy. Mode sampling is the fastest,
 *  but least theoretically robust, unbiased BF is the slowest
 *  (even hangs completely from time to time), but it is also the most
 *  theoretically robust.
 */
enum class ThermalizationAlgorithm {
  ModeSampling,
  BiasedBF,
  UnbiasedBF,
};

/// Defines how the number of events is determined.
enum class EventCounting {
  /// The desired number of events is simulated disregarding of whether an
  /// interaction took place.
  FixedNumber,
  /// Events are simulated until there are at least a given number of ensembles
  /// in which an interaction took place
  MinimumNonEmpty,
  /// Unused, only in the code for internal logic
  Invalid,
};

/// Determine how total cross sections for collision finding should be computed.
enum class TotalCrossSectionStrategy {
  /// Sum the existing partial contributions
  BottomUp,
  /// Use parametrizations based on existing data, rescaling with AQM for
  /// unmeasured processes
  TopDown,
  /// Mix the two above, using the parametrizations only for measured processes,
  /// and summing up partials for unmeasured interactions
  TopDownMeasured,
};

/**
 *  Which pseudo-resonance fills the inelastic gap in the transition to string
 * region of cross sections. \see_key{key_CT_pseudoresonance_}
 */
enum class PseudoResonance {
  /// No pseudo-resonance is created
  None,
  /// Resonance of largest mass for all processes
  Largest,
  /// Resonance with the pole mass closest from the invariant mass of incoming
  /// particles for all processes
  Closest,
  /// Heaviest possible resonance from processes with at least one resonance in
  /// the incoming particles
  LargestFromUnstable,
  /// Closest resonance for a given mass from processes with at least one
  /// resonance in the incoming particles
  ClosestFromUnstable,
};

/// Possible methods to convert SMASH particle into fluid cells.
/// \see_key{key_MC_IC_type_}
enum class FluidizationType {
  /// Hypersurface crossed at a fixed proper time
  ConstantTau,
  /// Dynamic fluidization based on local densities
  Dynamic,
};

/// The different processes from where fluidizable particles are produced.
/// \see_key{key_MC_IC_fluidizable_processes}
// Because std::bitset does not handle enum classes, this is a simple enum.
enum IncludedFluidizableProcesses {
  From_Elastic = 0,
  From_Decay = 1,
  From_Inelastic = 2,
  From_SoftString = 3,
  From_HardString = 4,
};

typedef std::bitset<5> FluidizableProcessesBitSet;

/**
 * Allows to choose which kind of density to calculate.
 * The baryon density is necessary for the Skyrme potential.
 * For the symmetry potential one needs to know the isospin density.
 */
enum class DensityType {
  None = 0,
  Hadron = 1,
  Baryon = 2,
  BaryonicIsospin = 3,
  Pion = 4,
  Isospin3_tot = 5,
  Charge = 6,
  Strangeness = 7,
};

/// @cond
template <typename T>
using build_unique_ptr_ = std::unique_ptr<T, std::default_delete<T>>;
template <typename T>
using build_vector_ = std::vector<T, std::allocator<T>>;

using ActionPtr = build_unique_ptr_<Action>;
using ScatterActionPtr = build_unique_ptr_<ScatterAction>;
using ScatterActionMultiPtr = build_unique_ptr_<ScatterActionMulti>;
using ActionList = build_vector_<ActionPtr>;

using OutputPtr = build_unique_ptr_<OutputInterface>;
using OutputsList = build_vector_<OutputPtr>;

using ParticleList = build_vector_<ParticleData>;
using ParticleTypeList = build_vector_<ParticleType>;
using ParticleTypePtrList = build_vector_<ParticleTypePtr>;
using IsoParticleTypeList = build_vector_<IsoParticleType>;

template <typename T>
using ProcessBranchPtr = build_unique_ptr_<T>;
template <typename T>
using ProcessBranchList = build_vector_<ProcessBranchPtr<T>>;
using DecayBranchPtr = build_unique_ptr_<DecayBranch>;
using DecayBranchList = build_vector_<DecayBranchPtr>;
using CollisionBranchPtr = build_unique_ptr_<CollisionBranch>;
using CollisionBranchList = build_vector_<CollisionBranchPtr>;

using TabulationPtr = build_unique_ptr_<Tabulation>;
using ExperimentPtr = build_unique_ptr_<ExperimentBase>;
using DecayTypePtr = build_unique_ptr_<DecayType>;

/// @endcond

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_FORWARDDECLARATIONS_H_
