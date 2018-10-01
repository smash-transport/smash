/*
 *
 *    Copyright (c) 2017-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"

#include "setup.h"

#include "../include/smash/formfactors.h"
#include "../include/smash/isoparticletype.h"
#include "../include/smash/kinematics.h"
#include "../include/smash/particletype.h"

using namespace smash;

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
      "π       0.138   7.7e-9  -   111     211\n"
      "ρ       0.776   0.149   -   113     213\n"
      "f₂      1.275   0.185   +   225\n"
      "N       0.938   0       +  2112    2212\n"
      "Δ       1.232   0.117   +  1114    2114    2214    2224\n");
}

TEST(init_decay_channels) {
  const std::string decays_input(
      "ρ \n"
      "1.0\t1\tπ π\n \n"
      "f₂ \n"
      "1.0\t0\tρ ρ\n"
      "Δ \n"
      "1.0\t1\tN π\n");
  DecayModes::load_decaymodes(decays_input);
  ParticleType::check_consistency();
}

TEST(pp_DeltaDelta_integral) {
  const ParticleType &proton = ParticleType::find(0x2212);
  const ParticleType &delta_plus = ParticleType::find(0x2214);
  const ParticleType &pi = ParticleType::find(0x211);
  const double mpi = pi.mass();
  const double mN = proton.mass();
  const double mD0 = delta_plus.mass();
  const double GD0 = delta_plus.width_at_pole();
  FUZZY_COMPARE(mpi, 0.138);
  FUZZY_COMPARE(mN, 0.938);
  FUZZY_COMPARE(mD0, 1.232);
  FUZZY_COMPARE(GD0, 0.117);
  COMPARE_ABSOLUTE_ERROR(pCM(mD0, mN, mpi), 0.228161951, 1.e-9);
  COMPARE_ABSOLUTE_ERROR(delta_plus.total_width(2.0), 0.400452404, 1.e-9);
  COMPARE_ABSOLUTE_ERROR(delta_plus.spectral_function_no_norm(1.3), 2.343254132,
                         1.e-9);
  COMPARE_ABSOLUTE_ERROR(delta_plus.spectral_function(1.3), 2.044684068, 1.e-9);
  std::vector<double> sqrts_list = {2.662, 3.162, 3.662, 4.162,
                                    4.662, 5.162, 5.662, 50.};
  std::vector<double> integrals2d = {0.180796, 0.539887, 0.858989, 1.15759,
                                     1.44408,  1.72266,  1.99572,  24.3628};
  for (size_t i = 0; i < sqrts_list.size(); i++) {
    const double integral_from_SMASH =
        delta_plus.iso_multiplet()->get_integral_RR(delta_plus, sqrts_list[i]);
    double rel_err = 1.e-4;
    if (i == sqrts_list.size() - 2) {
      // the error gets larger because the SMASH tabulation is extrapolated
      rel_err = 2.e-3;
    }
    if (i == sqrts_list.size() - 1) {
      // the error is even larger, when extrapolating from 2-5 GeV to 50 GeV,
      // but in that region 10% is still ok.
      rel_err = 0.1;
    }
    COMPARE_RELATIVE_ERROR(integrals2d[i], integral_from_SMASH, rel_err)
        << " at sqrts = " << sqrts_list[i];
  }
  /* Integrals were independently computed with the following Mathematica code:
     pCM[m_, m1_, m2_] := Sqrt[(m^2 - (m1 + m2)^2)*(m^2 - (m1 - m2)^2)]/(2 m)
     mN = 0.938; mpi = 0.138; mD0 = 1.232; GD0 = 0.117; La = 0.197327053;
     G[m_] := GD0* mD0/m*(pCM[m, mN, mpi]/ pCM[mD0, mN, mpi])^3 *
       (pCM[mD0, mN, mpi]^2 + La^2)/(pCM[m, mN, mpi]^2 + La^2)
     Anonorm[m_] := 2/Pi*m^2 G[m]/((m^2 - mD0^2)^2 + (m*G[m])^2)
     norm = NIntegrate[Anonorm[m], {m, mN + mpi, Infinity}]
     A[m_] := Anonorm[m]/norm
     UnstableIntegrand[m_, m1_, m2_] :=
       A[m1]*A[m2]*HeavisideTheta[m - m1 - m2] *pCM[m, m1, m2]
     UnstableIntegral[m_] := NIntegrate[
       UnstableIntegrand[m, m1, m2], {m1, mN + mpi, m - (mN + mpi)},
         {m2, mN + mpi, m - (mN + mpi)}]
   */
}

TEST(f2_width) {
  const ParticleType &f2 = ParticleType::find(0x225);
  const ParticleType &rho = ParticleType::find(0x113);
  const ParticleType &pi = ParticleType::find(0x111);
  const double mpi = pi.mass();
  FUZZY_COMPARE(pi.mass(), 0.138);
  FUZZY_COMPARE(rho.mass(), 0.776);
  FUZZY_COMPARE(rho.width_at_pole(), 0.149);
  FUZZY_COMPARE(f2.mass(), 1.275);
  FUZZY_COMPARE(f2.width_at_pole(), 0.185);
  COMPARE_RELATIVE_ERROR(rho.iso_multiplet()->get_integral_RR(rho, 2.0),
                         0.417261, 5e-4);
  COMPARE_ABSOLUTE_ERROR(post_ff_sqr(2.0, f2.mass(), mpi * 4.0, 0.6),
                         0.00366995, 1e-8);

  std::vector<double> f2mass_list = {0.862, 1.162, 1.462, 1.762, 2.062, 2.362,
                                     2.662, 2.962, 3.262, 3.562, 3.862, 4.162,
                                     4.462, 4.762, 5.062, 5.362, 5.662};
  std::vector<double> f2_width_list = {
      0.00193368,   0.158777,    0.223674,     0.155464,    0.0480971,
      0.0156635,    0.00577278,  0.00237419,   0.00106874,  0.000518157,
      0.000267202,  0.000145138, 0.0000824081, 0.000048615, 0.0000296518,
      0.0000186238, 0.0000120054};
  for (size_t i = 0; i < f2mass_list.size(); i++) {
    const double f2width_m = f2.total_width(f2mass_list[i]);
    double rel_err = 5.e-3;
    // Extrapolation leads to larger errors
    if (i > 10) {
      rel_err = 5.e-2;
    }
    if (i > 13) {
      rel_err = 1.e-1;
    }
    COMPARE_RELATIVE_ERROR(f2_width_list[i], f2width_m, rel_err)
        << " at mass = " << f2mass_list[i];
  }
  /* Integrals were independently computed with the following Mathematica code:
     mpi = 0.138; mrho0 = 0.776; Grho0 = 0.149; m0f2 = 1.275; G0f2 = 0.185;
     La = 0.197327053; Lcut = 0.6;
     pCM[m_, m1_, m2_] := Sqrt[(m^2 - (m1 + m2)^2)*(m^2 - (m1 - m2)^2)]/(2 m)
     Grho[m_] := Grho0*mrho0/ m*(pCM[m, mpi, mpi]/ pCM[mrho0, mpi, mpi])^3 *
       (pCM[mrho0, mpi, mpi]^2 + La^2)/(pCM[m, mpi, mpi]^2 + La^2)
     Anonormrho[m_] := 2/Pi*m^2 Grho[m]/((m^2 - mrho0^2)^2 + (m*Grho[m])^2)
     normrho = NIntegrate[Anonormrho[m], {m, 2*mpi, Infinity}]
     Arho[m_] := Anonormrho[m]/normrho
     Integralf2[m_] := NIntegrate[ Arho[m1] Arho[m2] pCM[m, m1, m2]
       HeavisideTheta[m - m1 - m2],
       {m1, 2 mpi, m - 2 mpi}, {m2, 2 mpi, m - 2 mpi}]
     ffPost[m_] := (Lcut^4 + (16 mpi^2 - m0f2^2)^2/ 4)/
       (Lcut^4 + (m^2 - (16 mpi^2 + m0f2^2)/2)^2)
     Gf2[m_] := G0f2 *m0f2/m*(ffPost[m]/ffPost[m0f2])^2 *
       Integralf2[m]/Integralf2[m0f2]
   */
}
