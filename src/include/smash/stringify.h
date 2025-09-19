/*
 *
 *    Copyright (c) 2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_STRINGIFY_H_
#define SRC_INCLUDE_SMASH_STRINGIFY_H_

#include <string>
#include <vector>

#include "einhard.hpp"

#include "forwarddeclarations.h"

namespace smash {

/**
 * Convert a ThermodynamicQuantity enum value to its corresponding string.
 *
 * \param[in] quantity The ThermodynamicQuantity enum value to convert.
 *
 * \return std::string Corresponding string representation.
 * \throws std::invalid_argument If the enum value is unhandled.
 */
std::string to_string(ThermodynamicQuantity quantity);

/**
 * Convert a CalculationFrame enum value to its corresponding string.
 *
 * \param[in] frame The CalculationFrame enum value to convert.
 *
 * \return std::string Corresponding string representation.
 * \throws std::invalid_argument If the enum value is unhandled.
 */
std::string to_string(CalculationFrame frame);

/**
 * Convert a FermiMotion enum value to its corresponding string.
 *
 * \param[in] motion The FermiMotion enum value to convert.
 *
 * \return std::string Corresponding string representation.
 * \throws std::invalid_argument If the enum value is unhandled.
 */
std::string to_string(FermiMotion motion);

/**
 * Convert a DensityType enum value to its corresponding string.
 *
 * \param[in] type The DensityType enum value to convert.
 *
 * \return std::string Corresponding string representation.
 * \throws std::invalid_argument If the enum value is unhandled.
 */
std::string to_string(DensityType type);

/**
 * Convert an ExpansionMode enum value to its corresponding string.
 *
 * \param[in] mode The ExpansionMode enum value to convert.
 *
 * \return std::string Corresponding string representation.
 * \throws std::invalid_argument If the enum value is unhandled.
 */
std::string to_string(ExpansionMode mode);

/**
 * Convert a DerivativesMode enum value to its corresponding string.
 *
 * \param[in] mode The DerivativesMode enum value to convert.
 *
 * \return std::string Corresponding string representation.
 * \throws std::invalid_argument If the enum value is unhandled.
 */
std::string to_string(DerivativesMode mode);

/**
 * Convert a FieldDerivativesMode enum value to its corresponding string.
 *
 * \param[in] mode The FieldDerivativesMode enum value to convert.
 *
 * \return std::string Corresponding string representation.
 * \throws std::invalid_argument If the enum value is unhandled.
 */
std::string to_string(FieldDerivativesMode mode);

/**
 * Convert a SmearingMode enum value to its corresponding string.
 *
 * \param[in] mode The SmearingMode enum value to convert.
 *
 * \return std::string Corresponding string representation.
 * \throws std::invalid_argument If the enum value is unhandled.
 */
std::string to_string(SmearingMode mode);

/**
 * Convert a TimeStepMode enum value to its corresponding string.
 *
 * \param[in] mode The TimeStepMode enum value to convert.
 *
 * \return std::string Corresponding string representation.
 * \throws std::invalid_argument If the enum value is unhandled.
 */
std::string to_string(TimeStepMode mode);

/**
 * Convert a BoxInitialCondition enum value to its corresponding string.
 *
 * \param[in] cond The BoxInitialCondition enum value to convert.
 *
 * \return std::string Corresponding string representation.
 * \throws std::invalid_argument If the enum value is unhandled.
 */
std::string to_string(BoxInitialCondition cond);

/**
 * Convert a SphereInitialCondition enum value to its corresponding string.
 *
 * \param[in] cond The SphereInitialCondition enum value to convert.
 *
 * \return std::string Corresponding string representation.
 * \throws std::invalid_argument If the enum value is unhandled.
 */
std::string to_string(SphereInitialCondition cond);

/**
 * Convert a NNbarTreatment enum value to its corresponding string.
 *
 * \param[in] t The NNbarTreatment enum value to convert.
 *
 * \return std::string Corresponding string representation.
 * \throws std::invalid_argument If the enum value is unhandled.
 */
std::string to_string(NNbarTreatment t);

/**
 * Convert a Sampling enum value to its corresponding string.
 *
 * \param[in] s The Sampling enum value to convert.
 *
 * \return std::string Corresponding string representation.
 * \throws std::invalid_argument If the enum value is unhandled.
 */
std::string to_string(Sampling s);

/**
 * Convert a ThermalizationAlgorithm enum value to its corresponding string.
 *
 * \param[in] algo The ThermalizationAlgorithm enum value to convert.
 *
 * \return std::string Corresponding string representation.
 * \throws std::invalid_argument If the enum value is unhandled.
 */
std::string to_string(ThermalizationAlgorithm algo);

/**
 * Convert a CollisionCriterion enum value to its corresponding string.
 *
 * \param[in] c The CollisionCriterion enum value to convert.
 *
 * \return std::string Corresponding string representation.
 * \throws std::invalid_argument If the enum value is unhandled.
 */
std::string to_string(CollisionCriterion c);

/**
 * Convert a TotalCrossSectionStrategy enum value to its corresponding string.
 *
 * \param[in] s The TotalCrossSectionStrategy enum value to convert.
 *
 * \return std::string Corresponding string representation.
 * \throws std::invalid_argument If the enum value is unhandled.
 */
std::string to_string(TotalCrossSectionStrategy s);

/**
 * Convert a PseudoResonance enum value to its corresponding string.
 *
 * \param[in] p The PseudoResonance enum value to convert.
 *
 * \return std::string Corresponding string representation.
 * \throws std::invalid_argument If the enum value is unhandled.
 */
std::string to_string(PseudoResonance p);

/**
 * Convert a FluidizationType enum value to its corresponding string.
 *
 * \param[in] f The FluidizationType enum value to convert.
 *
 * \return std::string Corresponding string representation.
 * \throws std::invalid_argument If the enum value is unhandled.
 */
std::string to_string(FluidizationType f);

/**
 * Convert an OutputOnlyFinal enum value to its corresponding string.
 *
 * \param[in] o The OutputOnlyFinal enum value to convert.
 *
 * \return std::string Corresponding string representation.
 * \throws std::invalid_argument If the enum value is unhandled.
 */
std::string to_string(OutputOnlyFinal o);

/**
 * Convert a LogLevel enum value to its corresponding string.
 *
 * \param[in] level The LogLevel enum value to convert.
 *
 * \return std::string Corresponding string representation.
 * \throws std::invalid_argument If the enum value is unhandled.
 */
std::string to_string(einhard::LogLevel level);

/**
 * Convert a ReactionsBitSet to a vector of strings for all set reactions.
 *
 * \param[in] s The ReactionsBitSet to convert.
 *
 * \return std::vector<std::string> Vector of all set reaction names.
 */
std::vector<std::string> to_string(const ReactionsBitSet &s);

/**
 * Convert a MultiParticleReactionsBitSet to a vector of strings for all set
 * reactions.
 *
 * \param[in] s The MultiParticleReactionsBitSet to convert.
 *
 * \return std::vector<std::string> Vector of all set reaction names.
 */
std::vector<std::string> to_string(const MultiParticleReactionsBitSet &s);

/**
 * Convert a FluidizableProcessesBitSet to a vector of strings for all set
 * processes.
 *
 * \param[in] s The FluidizableProcessesBitSet to convert.
 *
 * \return std::vector<std::string> Vector of all set process names.
 */
std::vector<std::string> to_string(const FluidizableProcessesBitSet &s);

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_STRINGIFY_H_
