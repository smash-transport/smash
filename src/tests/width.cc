/*
 *
 *    Copyright (c) 2016-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <vir/test.h>  // This include has to be first

#include "setup.h"

#include "../include/smash/integrate.h"

using namespace smash;

TEST(init_particle_types) {
  // enable debugging output
  create_all_loggers(Configuration(""));

  Test::create_actual_particletypes();
  Test::create_actual_decaymodes();
}

TEST(width_Delta) {
  const ParticleType &t = ParticleType::find(0x2214);
  for (int i = 0; i < 100; i++) {
    double m = 1. + i * 0.01;
    double w = t.total_width(m);
    std::printf("%7.3f %7.3f \n", m, w);
  }
}

TEST(width_Roper) {
  const ParticleType &t = ParticleType::find(0x12212);  // N(1440)
  std::printf("%7.4f \n", t.min_mass_kinematic());
  const int nModes = t.decay_modes().decay_mode_list().size();
  for (int i = 0; i < 100; i++) {
    double m = 1. + i * 0.01;
    double wtot = t.total_width(m);
    std::printf("%7.4f %7.4f ", m, wtot);
    /* Print all partial widths. */
    ProcessBranchList<DecayBranch> partial = t.get_partial_widths(
        FourVector(m, 0., 0., 0.), ThreeVector(), WhichDecaymodes::All);
    for (const auto &branch : partial) {
      std::printf("%7.4f ", branch->weight());
    }
    // print zeros for missing modes
    for (int j = partial.size(); j < nModes; j++) {
      std::printf("%7.4f ", 0.);
    }
    std::printf("\n");
  }
}

/* Compare the out-width vs the integrated in-width,
 * according to equ. (2.60) in Effenberger's thesis,
 * for a given resonance type, decay branch and resonance mass. */
static void compare_in_vs_out_width(const ParticleType &t,
                                    const DecayBranch &the_branch,
                                    const double m_R) {
  // calculate out-width
  const double Gam_out = t.partial_width(m_R, &the_branch);

  // integrate in-width
  const ParticleTypePtrList pt = the_branch.particle_types();
  ParticleData p0 = ParticleData(*pt[0]);
  p0.set_4momentum(pt[0]->mass(), 0., 0., 0.);
  ParticleData p1 = ParticleData(*pt[1]);
  Integrator integrate;
  const double Gam_in = integrate(pt[1]->min_mass_kinematic(),
                                  m_R - pt[0]->mass(), [&](double m) {
                                    p1.set_4momentum(m, 0., 0., 0.);
                                    return t.get_partial_in_width(m_R, p0, p1) *
                                           pt[1]->spectral_function(m);
                                  });

  std::printf("width comparison at m=%5.3f: %8.6f %8.6f %8.6f %8.6f \n", m_R,
              Gam_out, Gam_in, std::abs(Gam_out - Gam_in), Gam_out / Gam_in);
  COMPARE_ABSOLUTE_ERROR(Gam_out, Gam_in, 7E-4);
}

TEST(Roper_in_vs_out_width) {
  // the Roper resonance, N*(1440)
  const ParticleType &t = ParticleType::find(0x12212);

  // find a decay mode with unstable daughter meson:
  // for the Roper, this should be the σN decay
  DecayBranch *the_branch = nullptr;
  for (const auto &mode : t.decay_modes().decay_mode_list()) {
    const auto pt = mode->particle_types();
    if (!pt[1]->is_stable() && !pt[1]->is_baryon()) {
      the_branch = mode.get();
      std::printf("found mode: %s %s \n", pt[0]->pdgcode().string().c_str(),
                  pt[1]->pdgcode().string().c_str());
    }
  }

  // test width for N* <-> σN at different resonance masses
  compare_in_vs_out_width(t, *the_branch, 1.25);
  compare_in_vs_out_width(t, *the_branch, 1.3);
  compare_in_vs_out_width(t, *the_branch, 1.4);
  compare_in_vs_out_width(t, *the_branch, 1.5);
  compare_in_vs_out_width(t, *the_branch, 1.6);
  compare_in_vs_out_width(t, *the_branch, 1.7);
  compare_in_vs_out_width(t, *the_branch, 1.8);
  compare_in_vs_out_width(t, *the_branch, 1.9);
  compare_in_vs_out_width(t, *the_branch, 2.0);
}

TEST(photon_widths) {
  // Test the photon widths that are used as a baseline for the dilepton Dalitz
  // widths.
  const ParticleType &photon = ParticleType::find(0x22);
  const ParticleType &pi0 = ParticleType::find(0x111);
  const ParticleType &eta = ParticleType::find(0x221);
  const ParticleType &etap = ParticleType::find(0x331);
  const ParticleType &omega = ParticleType::find(0x223);
  const ParticleType &phi = ParticleType::find(0x333);

  const double err = 1E-10;

  // The partial widths should be given by the width of the resonance times the
  // branching ratio. This test breaks whenever the branching ratios or width of
  // the involved particles are changed.
  COMPARE_ABSOLUTE_ERROR(pi0.get_partial_width(pi0.mass(), photon, photon),
                         7.699999749e-09, err);
  COMPARE_ABSOLUTE_ERROR(eta.get_partial_width(eta.mass(), photon, photon),
                         5.203024193548386e-7, err);
  COMPARE_ABSOLUTE_ERROR(etap.get_partial_width(etap.mass(), photon, photon),
                         4.3889449263667536e-6, err);

  COMPARE_ABSOLUTE_ERROR(omega.get_partial_width(omega.mass(), pi0, photon),
                         0.0007194189448, err);
  COMPARE_ABSOLUTE_ERROR(phi.get_partial_width(phi.mass(), pi0, photon),
                         5.4068538571729e-6, err);
}
