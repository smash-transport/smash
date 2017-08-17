/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_PDGCODE_CONSTANTS_H_
#define SRC_INCLUDE_PDGCODE_CONSTANTS_H_

namespace Smash {
/**
 * Constants representing PDG codes.
 *
 * '_p' is short for '+', '_pp' for '++', '_m' for '-' and '_z' for '0'.
 */
namespace pdg {

constexpr int invalid = 0x0;

constexpr int photon = 0x22;

constexpr int p = 0x2212;
constexpr int n = 0x2112;

constexpr int N1535_p = 0x22212;
constexpr int N1535_z = 0x22112;

constexpr int Delta_pp = 0x2224;
constexpr int Delta_p = 0x2214;
constexpr int Delta_z = 0x2114;
constexpr int Delta_m = 0x1114;

constexpr int Lambda = 0x3122;
constexpr int Sigma_p = 0x3222;
constexpr int Sigma_z = 0x3212;
constexpr int Sigma_m = 0x3112;
constexpr int Xi_z = 0x3322;
constexpr int Xi_m = 0x3312;
constexpr int Omega_m = 0x3334;

constexpr int pi_p = 0x211;
constexpr int pi_z = 0x111;
constexpr int pi_m = -0x211;

constexpr int K_p = 0x321;
constexpr int K_z = 0x311;
constexpr int Kbar_z = -0x311;
constexpr int K_m = -0x321;

constexpr int eta = 0x221;
constexpr int omega = 0x223;

constexpr int rho_p = 0x213;
constexpr int rho_z = 0x113;
constexpr int rho_m = -0x213;

constexpr int h1 = 0x10223;

}  // namespace pdg

/// Pack two int32_t into an uint64_t.
/// This is useful for switch statements on pairs.
constexpr uint64_t pack(int32_t x, int32_t y) {
  return (static_cast<uint64_t>(static_cast<uint32_t>(x)) << 32) |
         static_cast<uint64_t>(static_cast<uint32_t>(y));
  //^ Casting to an intermediate 32-bit integer is important!
}

}  // namespace Smash

#endif  // SRC_INCLUDE_PDGCODE_CONSTANTS_H_
