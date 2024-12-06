/*
 *
 *    Copyright (c) 2016-2021,2023
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_PDGCODE_CONSTANTS_H_
#define SRC_INCLUDE_SMASH_PDGCODE_CONSTANTS_H_

namespace smash {
/**
 * Constants representing PDG codes.
 *
 * '_p' is short for '+', '_pp' for '++', '_m' for '-' and '_z' for '0'.
 */
namespace pdg {

/// Invalid particle.
constexpr int invalid = 0x0;

/// Photon.
constexpr int photon = 0x22;

/// Proton.
constexpr int p = 0x2212;
/// Neutron.
constexpr int n = 0x2112;

/// N(1520)⁺.
constexpr int N1520_p = 0x2124;
/// N(1520)⁰.
constexpr int N1520_z = 0x1214;

/// N(1535)⁺.
constexpr int N1535_p = 0x22212;
/// N(1535)⁰.
constexpr int N1535_z = 0x22112;

/// Δ⁺⁺.
constexpr int Delta_pp = 0x2224;
/// Δ⁺.
constexpr int Delta_p = 0x2214;
/// Δ⁰.
constexpr int Delta_z = 0x2114;
/// Δ⁻.
constexpr int Delta_m = 0x1114;

/// Λ.
constexpr int Lambda = 0x3122;
/// Σ⁺.
constexpr int Sigma_p = 0x3222;
/// Σ⁰.
constexpr int Sigma_z = 0x3212;
/// Σ⁻.
constexpr int Sigma_m = 0x3112;
/// Ξ⁰.
constexpr int Xi_z = 0x3322;
/// Ξ⁻.
constexpr int Xi_m = 0x3312;
/// Ω⁻.
constexpr int Omega_m = 0x3334;

/// π⁺.
constexpr int pi_p = 0x211;
/// π⁰.
constexpr int pi_z = 0x111;
/// π⁻.
constexpr int pi_m = -0x211;

/// K⁺.
constexpr int K_p = 0x321;
/// K⁰.
constexpr int K_z = 0x311;
/// K̄⁰.
constexpr int Kbar_z = -0x311;
/// K̄⁻.
constexpr int K_m = -0x321;

/// η.
constexpr int eta = 0x221;
/// ω.
constexpr int omega = 0x223;

/// ρ⁺.
constexpr int rho_p = 0x213;
/// ρ⁰.
constexpr int rho_z = 0x113;
/// ρ⁻.
constexpr int rho_m = -0x213;

/// h₁(1170).
constexpr int h1 = 0x10223;

/*
 * Constants representing PDG codes of nuclei.
 */

/// Deuteron.
constexpr int64_t deuteron = 0x1000010020;
/// Anti-deuteron in decimal digits.
constexpr int64_t antideuteron = -0x1000010020;
/// Deuteron-prime resonance.
constexpr int64_t dprime = 0x1000010021;
/// Triton.
constexpr int64_t triton = 0x1000010030;
/// Anti-triton.
constexpr int64_t antitriton = -0x1000010030;
/// He-3
constexpr int64_t he3 = 0x1000020030;
/// Anti-He-3
constexpr int64_t antihe3 = -0x1000020030;
/// Hypertriton
constexpr int64_t hypertriton = 0x1010010030;
/// Anti-Hypertriton
constexpr int64_t antihypertriton = -0x1010010030;

}  // namespace pdg

/**
 * Pack two int32_t into an uint64_t.
 * This is useful for switch statements on pairs.
 *
 * \param x First integer to be packed.
 * \param y Second integer to be packed.
 * \return Combined integer.
 */
constexpr uint64_t pack(int32_t x, int32_t y) {
  return (static_cast<uint64_t>(static_cast<uint32_t>(x)) << 32) |
         static_cast<uint64_t>(static_cast<uint32_t>(y));
  //^ Casting to an intermediate 32-bit integer is important!
}

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_PDGCODE_CONSTANTS_H_
