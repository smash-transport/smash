/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "setup.h"

#include "../include/integrate.h"

using namespace Smash;

TEST(init_particle_types) {
  Test::create_actual_particletypes();
  Test::create_actual_decaymodes();
}

TEST(width_Delta) {
  const ParticleType &t = ParticleType::find(0x2214);
  for (int i = 0; i < 100; i++) {
    float m = 1. + i*0.01;
    float w = t.total_width (m);
    printf("%7.3f %7.3f \n", m, w);
  }
}

TEST(width_Roper) {
  const ParticleType &t = ParticleType::find(0x202212);
  printf("%7.4f \n",t.minimum_mass());
  for (int i = 0; i < 100; i++) {
    float m = 1. + i*0.01;
    float wtot = t.total_width (m);
    printf("%7.4f %7.4f ", m, wtot);
    /* Print all partial widths. */
    ProcessBranchList<DecayBranch> partial = t.get_partial_widths(m);
    for (const auto &branch : partial) {
      printf("%7.4f ", branch->weight());
    }
    printf("\n");
  }
}


/* Compare the out-width vs the integrated in-width,
 * according to equ. (2.60) in Effenberger's thesis,
 * for a given resonance type, decay branch and resonance mass. */
static void compare_in_vs_out_width(const ParticleType &t,
                                    const DecayBranch &the_branch,
                                    const float m_R) {
  // calculate out-width
  const float Gam_out = t.partial_width(m_R, &the_branch);

  // integrate in-width
  const ParticleTypePtrList pt = the_branch.particle_types();
  ParticleData p0 = ParticleData(*pt[0]);
  p0.set_4momentum(pt[0]->mass(), 0., 0., 0.);
  ParticleData p1 = ParticleData(*pt[1]);
  Integrator integrate;
  const float Gam_in = integrate(pt[1]->minimum_mass(), m_R - pt[0]->mass(),
                                  [&](float m) {
                                    p1.set_4momentum(m, 0., 0., 0.);
                                    return t.get_partial_in_width(m_R, p0, p1)
                                            * pt[1]->spectral_function(m);
                                  });

//   COMPARE(Gam_out, Gam_in);
  printf("width comparison at m=%7.4f: %7.4f %7.4f \n", m_R, Gam_out, Gam_in);
}


TEST(Roper_in_vs_out_width) {
  // the Roper resonance, N*(1440)
  const ParticleType &t = ParticleType::find(0x202212);

  // find a decay mode with unstable daughter meson:
  // for the Roper, this should be the σN decay
  DecayBranch *the_branch = nullptr;
  for (const auto &mode : t.decay_modes().decay_mode_list()) {
    const auto pt = mode->particle_types();
    if (!pt[1]->is_stable() && !pt[1]->is_baryon()) {
      the_branch = mode.get();
      printf("found mode: %s %s \n", pt[0]->pdgcode().string().c_str(),
                                     pt[1]->pdgcode().string().c_str());
    }
  }

  // test width for N* <-> σN at different resonance masses
  compare_in_vs_out_width(t, *the_branch, 1.3);
  compare_in_vs_out_width(t, *the_branch, 1.4);
  compare_in_vs_out_width(t, *the_branch, 1.5);
  compare_in_vs_out_width(t, *the_branch, 1.6);
}
