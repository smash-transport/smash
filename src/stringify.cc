/*
 *
 *    Copyright (c) 2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/stringify.h"

#include <stdexcept>
#include <string_view>

namespace smash {

// Throw because an enum value was not handled in the \c to_string function.
static void throw_unhandled_enum(std::string_view enum_name, int value) {
  throw std::invalid_argument("Unhandled " + std::string(enum_name) +
                              " enum value " + std::to_string(value) +
                              " passed to conversion function to_string().");
}

std::string to_string(IncludedReactions reaction) {
  switch (reaction) {
    case IncludedReactions::All:
      return "All";
    case IncludedReactions::Elastic:
      return "Elastic";
    case IncludedReactions::NN_to_NR:
      return "NN_to_NR";
    case IncludedReactions::NN_to_DR:
      return "NN_to_DR";
    case IncludedReactions::KN_to_KN:
      return "KN_to_KN";
    case IncludedReactions::KN_to_KDelta:
      return "KN_to_KDelta";
    case IncludedReactions::Strangeness_exchange:
      return "Strangeness_exchange";
    case IncludedReactions::NNbar:
      return "NNbar";
    case IncludedReactions::PiDeuteron_to_NN:
      return "PiDeuteron_to_NN";
    case IncludedReactions::PiDeuteron_to_pidprime:
      return "PiDeuteron_to_pidprime";
    case IncludedReactions::NDeuteron_to_Ndprime:
      return "NDeuteron_to_Ndprime";
  }
  throw_unhandled_enum("IncludedReactions", static_cast<int>(reaction));
}

std::string to_string(IncludedMultiParticleReactions reaction) {
  switch (reaction) {
    case IncludedMultiParticleReactions::Meson_3to1:
      return "Meson_3to1";
    case IncludedMultiParticleReactions::Deuteron_3to2:
      return "Deuteron_3to2";
    case IncludedMultiParticleReactions::NNbar_5to2:
      return "NNbar_5to2";
    case IncludedMultiParticleReactions::A3_Nuclei_4to2:
      return "A3_Nuclei_4to2";
  }
  throw_unhandled_enum("IncludedMultiParticleReactions",
                       static_cast<int>(reaction));
}

std::string to_string(ThermodynamicQuantity quantity) {
  switch (quantity) {
    case ThermodynamicQuantity::EckartDensity:
      return "rho_eckart";
    case ThermodynamicQuantity::Tmn:
      return "tmn";
    case ThermodynamicQuantity::TmnLandau:
      return "tmn_landau";
    case ThermodynamicQuantity::LandauVelocity:
      return "landau_velocity";
    case ThermodynamicQuantity::j_QBS:
      return "j_QBS";
  }
  throw_unhandled_enum("ThermodynamicQuantity", static_cast<int>(quantity));
}

std::string to_string(CalculationFrame frame) {
  switch (frame) {
    case CalculationFrame::CenterOfVelocity:
      return "center of velocity";
    case CalculationFrame::CenterOfMass:
      return "center of mass";
    case CalculationFrame::FixedTarget:
      return "fixed target";
  }
  throw_unhandled_enum("CalculationFrame", static_cast<int>(frame));
}

std::string to_string(FermiMotion motion) {
  switch (motion) {
    case FermiMotion::Off:
      return "off";
    case FermiMotion::On:
      return "on";
    case FermiMotion::Frozen:
      return "frozen";
  }
  throw_unhandled_enum("FermiMotion", static_cast<int>(motion));
}

std::string to_string(DensityType type) {
  switch (type) {
    case DensityType::Hadron:
      return "hadron";
    case DensityType::Baryon:
      return "baryon";
    case DensityType::BaryonicIsospin:
      return "baryonic isospin";
    case DensityType::Pion:
      return "pion";
    case DensityType::Isospin3_tot:
      return "total isospin";
    case DensityType::None:
      return "none";
    case DensityType::Charge:
      return "charge";
    case DensityType::Strangeness:
      return "strangeness";
  }
  throw_unhandled_enum("DensityType", static_cast<int>(type));
}

std::string to_string(ExpansionMode mode) {
  switch (mode) {
    case ExpansionMode::NoExpansion:
      return "NoExpansion";
    case ExpansionMode::MasslessFRW:
      return "MasslessFRW";
    case ExpansionMode::MassiveFRW:
      return "MassiveFRW";
    case ExpansionMode::Exponential:
      return "Exponential";
  }
  throw_unhandled_enum("ExpansionMode", static_cast<int>(mode));
}

std::string to_string(DerivativesMode mode) {
  switch (mode) {
    case DerivativesMode::CovariantGaussian:
      return "Covariant Gaussian";
    case DerivativesMode::FiniteDifference:
      return "Finite difference";
    case DerivativesMode::Off:
      return "Off";
  }
  throw_unhandled_enum("DerivativesMode", static_cast<int>(mode));
}

std::string to_string(FieldDerivativesMode mode) {
  switch (mode) {
    case FieldDerivativesMode::ChainRule:
      return "Chain Rule";
    case FieldDerivativesMode::Direct:
      return "Direct";
  }
  throw_unhandled_enum("FieldDerivativesMode", static_cast<int>(mode));
}

std::string to_string(SmearingMode mode) {
  switch (mode) {
    case SmearingMode::CovariantGaussian:
      return "Covariant Gaussian";
    case SmearingMode::Discrete:
      return "Discrete";
    case SmearingMode::Triangular:
      return "Triangular";
  }
  throw_unhandled_enum("SmearingMode", static_cast<int>(mode));
}

std::string to_string(TimeStepMode mode) {
  switch (mode) {
    case TimeStepMode::None:
      return "None";
    case TimeStepMode::Fixed:
      return "Fixed";
  }
  throw_unhandled_enum("TimeStepMode", static_cast<int>(mode));
}

std::string to_string(BoxInitialCondition cond) {
  switch (cond) {
    case BoxInitialCondition::ThermalMomentaBoltzmann:
      return "thermal momenta";
    case BoxInitialCondition::ThermalMomentaQuantum:
      return "thermal momenta quantum";
    case BoxInitialCondition::PeakedMomenta:
      return "peaked momenta";
  }
  throw_unhandled_enum("BoxInitialCondition", static_cast<int>(cond));
}

std::string to_string(SphereInitialCondition cond) {
  switch (cond) {
    case SphereInitialCondition::ThermalMomentaBoltzmann:
      return "thermal momenta";
    case SphereInitialCondition::ThermalMomentaQuantum:
      return "thermal momenta quantum";
    case SphereInitialCondition::IC_ES:
      return "IC_ES";
    case SphereInitialCondition::IC_1M:
      return "IC_1M";
    case SphereInitialCondition::IC_2M:
      return "IC_2M";
    case SphereInitialCondition::IC_Massive:
      return "IC_Massive";
  }
  throw_unhandled_enum("SphereInitialCondition", static_cast<int>(cond));
}

std::string to_string(NNbarTreatment t) {
  switch (t) {
    case NNbarTreatment::NoAnnihilation:
      return "no annihilation";
    case NNbarTreatment::Resonances:
      return "resonances";
    case NNbarTreatment::TwoToFive:
      return "two to five";
    case NNbarTreatment::Strings:
      return "strings";
  }
  throw_unhandled_enum("NNbarTreatment", static_cast<int>(t));
}

std::string to_string(Sampling s) {
  switch (s) {
    case Sampling::Quadratic:
      return "quadratic";
    case Sampling::Custom:
      return "custom";
    case Sampling::Uniform:
      return "uniform";
  }
  throw_unhandled_enum("Sampling", static_cast<int>(s));
}

std::string to_string(ThermalizationAlgorithm algo) {
  switch (algo) {
    case ThermalizationAlgorithm::ModeSampling:
      return "mode sampling";
    case ThermalizationAlgorithm::BiasedBF:
      return "biased BF";
    case ThermalizationAlgorithm::UnbiasedBF:
      return "unbiased BF";
  }
  throw_unhandled_enum("ThermalizationAlgorithm", static_cast<int>(algo));
}

std::string to_string(CollisionCriterion c) {
  switch (c) {
    case CollisionCriterion::Geometric:
      return "Geometric";
    case CollisionCriterion::Stochastic:
      return "Stochastic";
    case CollisionCriterion::Covariant:
      return "Covariant";
  }
  throw_unhandled_enum("CollisionCriterion", static_cast<int>(c));
}

std::string to_string(TotalCrossSectionStrategy s) {
  switch (s) {
    case TotalCrossSectionStrategy::BottomUp:
      return "BottomUp";
    case TotalCrossSectionStrategy::TopDown:
      return "TopDown";
    case TotalCrossSectionStrategy::TopDownMeasured:
      return "TopDownMeasured";
  }
  throw_unhandled_enum("TotalCrossSectionStrategy", static_cast<int>(s));
}

std::string to_string(PseudoResonance p) {
  switch (p) {
    case PseudoResonance::None:
      return "None";
    case PseudoResonance::Largest:
      return "Largest";
    case PseudoResonance::Closest:
      return "Closest";
    case PseudoResonance::LargestFromUnstable:
      return "LargestFromUnstable";
    case PseudoResonance::ClosestFromUnstable:
      return "ClosestFromUnstable";
  }
  throw_unhandled_enum("PseudoResonance", static_cast<int>(p));
}

std::string to_string(FluidizationType f) {
  switch (f) {
    case FluidizationType::ConstantTau:
      return "Constant_Tau";
    case FluidizationType::Dynamic:
      return "Dynamic";
  }
  throw_unhandled_enum("FluidizationType", static_cast<int>(f));
}

std::string to_string(OutputOnlyFinal o) {
  switch (o) {
    case OutputOnlyFinal::Yes:
      return "Yes";
    case OutputOnlyFinal::No:
      return "No";
    case OutputOnlyFinal::IfNotEmpty:
      return "IfNotEmpty";
  }
  throw_unhandled_enum("OutputOnlyFinal", static_cast<int>(o));
}

std::string to_string(IncludedFluidizableProcesses p) {
  switch (p) {
    case IncludedFluidizableProcesses::From_Elastic:
      return "Elastic";
    case IncludedFluidizableProcesses::From_Decay:
      return "Decay";
    case IncludedFluidizableProcesses::From_Inelastic:
      return "Inelastic";
    case IncludedFluidizableProcesses::From_SoftString:
      return "SoftString";
    case IncludedFluidizableProcesses::From_HardString:
      return "HardString";
  }
  throw_unhandled_enum("IncludedFluidizableProcesses", static_cast<int>(p));
}

std::string to_string(einhard::LogLevel level) {
  switch (level) {
    case einhard::LogLevel::ALL:
      return "ALL";
    case einhard::LogLevel::TRACE:
      return "TRACE";
    case einhard::LogLevel::DEBUG:
      return "DEBUG";
    case einhard::LogLevel::INFO:
      return "INFO";
    case einhard::LogLevel::WARN:
      return "WARN";
    case einhard::LogLevel::ERROR:
      return "ERROR";
    case einhard::LogLevel::FATAL:
      return "FATAL";
    case einhard::LogLevel::OFF:
      return "OFF";
  }
  throw_unhandled_enum("IncludedFluidizableProcesses", static_cast<int>(level));
}

}  // namespace smash
