#include "unittest.h"  // This include has to be first

#include "setup.h"

#include <cmath>
#include <fstream>
#include <typeinfo>

#include "../include/smash/cxx14compat.h"
#include "../include/smash/decayaction.h"
#include "../include/smash/decaymodes.h"

using namespace smash;

double costheta_12(FourVector f1, FourVector f2) { return 2.0; }

double kinematic_g(const double x, const double y, const double z,
                   const double u, const double v, const double w) {
  // equation IV.5.25 in TODO
  return (-2 * u * v * w + 2 * v * v * w + 2 * v * w * w + 2 * u * v * x -
          2 * v * w * x + 2 * u * w * y - 2 * v * w * y - 2 * u * x * y -
          2 * v * x * y - 2 * w * x * y + 2 * x * x * y + 2 * x * y * y +
          2 * u * u * z - 2 * u * v * z - 2 * u * w * z - 2 * v * w * z -
          2 * u * x * z + 2 * w * x * z - 2 * u * y * z + 2 * v * y * z -
          2 * x * y * z + 2 * u * z * z) /
         2.;
}

double kinematic_l_sqrt(const double x, const double y, const double z) {
  return std::sqrt(std::pow(x - y - z, 2) - 4 * y * z);
}

ParticleList sample_3body_phase_space_jonas(const ParticleData &in,
                                            const ParticleList &out) {
  ParticleData incoming = in;
  // create particledata for now with bogus id. set it later in production code

  ParticleList outgoing = out;

  const double sqrts = incoming.effective_mass();
  std::cout << sqrts << std::endl;
  const double s = sqrts * sqrts;
  const double m1 = out[0].type().mass();
  const double m2 = out[1].type().mass();
  const double m3 = out[2].type().mass();

  const double s1_min = std::pow(m1 + m2, 2);
  const double s1_max = std::pow(sqrts - m3, 2);
  const double s2_min = std::pow(m2 + m3, 2);
  const double s2_max = std::pow(sqrts - m3, 2);
  std::cout << "Limits for s1 ";
  std::cout << s1_min << " " << s1_max << std::endl;
  bool done = false;
  double s1, s2;
  while (!done) {
    s1 = random::uniform(s1_min, s1_max);
    s2 = random::uniform(s2_min, s2_max);

    if (kinematic_g(s1, s2, sqrts * sqrts, m2 * m2, m1 * m1, m3 * m3) < 0)
      done = true;
  }

  const double s3 = s + m1 * m1 + m2 * m2 + m3 * m3 - s1 - s2;

  // work in R23. Transform incoming particle to frame where 2, 3 are at rest.
  // Its energy in this frame is known, therefore we can boost in z direction

  // first transform incoming particle to frame where p = (px, 0, 0).
  // this is not necessary, but makes the math a bit easier.

  double t1, t2;
  auto vec = incoming.momentum().threevec();
  std::cout << "p0 before rotation: " << incoming.momentum() << std::endl;
  if (vec[0] != 0) {
    t1 = std::atan(vec[2] / vec[0]);
  } else
    t1 = 0;  // ???
  std::cout << t1 << std::endl;
  vec.rotate_around_y(t1);
  if (vec[0] != 0) {
    t2 = std::atan(-vec[1] / vec[0]);
  } else
    t2 = 0;  // ??
  vec.rotate_around_z(t2);
  incoming.set_3momentum(vec);
  std::cout << "p0 after rotation: " << incoming.momentum() << std::endl;

  // boost to rest frame of particle 2 and 3. See eq. V.1.7
  const double E23 = (s + s2 - m1 * m1) / (2 * std::sqrt(s2));
  const double E = incoming.momentum()[0];
  const double beta =
      (E * incoming.momentum()[1] +
       std::sqrt(-(std::pow(E, 2) * std::pow(E23, 2)) + std::pow(E23, 4) +
                 std::pow(E23, 2) * std::pow(incoming.momentum()[1], 2))) /
      (std::pow(E23, 2) + std::pow(incoming.momentum()[1], 2));
  incoming.boost_momentum(ThreeVector(beta, 0, 0));
  std::cout << "p0 after boost: " << incoming.momentum() << std::endl;

  // for testing
  const double E23_expected = (s + s2 - m1 * m1) / std::sqrt(4 * s2);
  std::cout << "Expected energy: " << E23_expected << std::endl;
  const double P23_expected =
      kinematic_l_sqrt(s, s2, m1 * m1) / std::sqrt(4 * s2);
  const double P23_actual = incoming.momentum().abs3();
  std::cout << "P0 difference " << P23_expected - P23_actual << std::endl;

  // particle energies
  const double E23_1 = (s - s2 - m1 * m1) / std::sqrt(4 * s2);
  const double E23_2 = (s2 + m2 * m2 - m3 * m3) / std::sqrt(4 * s2);
  const double E23_3 = (s2 + m3 * m3 - m2 * m2) / std::sqrt(4 * s2);
  // particle three-momenta magnitudes
  const double p23_1 = kinematic_l_sqrt(s, s2, m1 * m1) / std::sqrt(4 * s2);
  const double p23_2 =
      kinematic_l_sqrt(s2, m2 * m2, m3 * m3) / std::sqrt(4 * s2);

  // now calculate angle between p1 and p2
  const double cos_theta_12 = ((s - s2 - m1 * m1) * (s2 + m2 * m2 - m3 * m3) +
                               2 * s2 * (m1 * m1 + m2 * m2 - s1)) /
                              (kinematic_l_sqrt(s, s2, m1 * m1) *
                               kinematic_l_sqrt(s2, m2 * m2, m3 * m3));

  Angles phitheta;
  phitheta.distribute_isotropically();
  phitheta.set_costheta(cos_theta_12);

  outgoing[0].set_4momentum(E23_1, incoming.momentum().threevec());
  outgoing[1].set_4momentum(E23_2, phitheta.threevec() * p23_2);
  outgoing[2].set_4momentum(E23_3, -phitheta.threevec() * p23_2);

  std::cout << "Momenta of outgoing particles: " << std::endl;
  for (auto &out: outgoing) { 
    std::cout << out.momentum() << std::endl;
  }

  // still in r23. check momenta
  //
  // TODO: Energy is not conserved. 
  std::cout << "p0 + p1 " << outgoing[0].momentum() - incoming.momentum()
            << std::endl;
  std::cout << "p2 + p3 " << outgoing[1].momentum() + outgoing[2].momentum()
            << std::endl;
  std::cout << "p0 + p1 + p2 + p3 "
            << incoming.momentum() - outgoing[0].momentum() -
                   outgoing[1].momentum() - outgoing[2].momentum()
            << std::endl;

  // rotate back (should not matter, since theta is sampled uniformly) and boost
  // back to lab system
  for (auto &out : outgoing) {
    out.momentum().threevec().rotate_around_z(-t2);
    out.momentum().threevec().rotate_around_y(-t1);
    out.boost_momentum(ThreeVector(-beta, 0, 0));
  }

  return outgoing;
}

void sample_3body_phase_space_dima(double srts, ParticleData &a,
                                   ParticleData &b, ParticleData &c) {
  const double m_a = a.type().mass(), m_b = b.type().mass(),
               m_c = c.type().mass();
  // sample mab from pCM(sqrt, mab, mc) pCM (mab, ma, mb) <= sqrts^2/4
  double mab, r, probability, pcm_ab, pcm;
  do {
    mab = random::uniform(m_a + m_b, srts - m_c);
    r = random::canonical();
    pcm = pCM(srts, mab, m_c);
    pcm_ab = pCM(mab, m_a, m_b);
    probability = pcm * pcm_ab * 4 / (srts * srts);
  } while (r > probability);
  Angles phitheta;
  phitheta.distribute_isotropically();
  c.set_4momentum(m_c, pcm * phitheta.threevec());
  const ThreeVector beta_cm =
      pcm * phitheta.threevec() / std::sqrt(pcm * pcm + mab * mab);

  phitheta.distribute_isotropically();
  a.set_4momentum(m_a, pcm_ab * phitheta.threevec());
  b.set_4momentum(m_b, -pcm_ab * phitheta.threevec());
  std::cout << "srts dima: " << srts << std::endl;
  // a.boost_momentum(beta_cm);
  // b.boost_momentum(beta_cm);
  // std::cout << a.momentum() + b.momentum() + c.momentum() << std::endl;
}

ParticleList one_to_three(ParticleList incoming_particles_,
                          ParticleList outgoing_particles_) {
  // this is an exact copy from the function found in decayaction.cc
  // We want to write out some values, this is the easiest way.
  const auto &log = logger<LogArea::DecayModes>();
  ParticleData &outgoing_a = outgoing_particles_[0];
  ParticleData &outgoing_b = outgoing_particles_[1];
  ParticleData &outgoing_c = outgoing_particles_[2];
  const ParticleType &outgoing_a_type = outgoing_a.type();
  const ParticleType &outgoing_b_type = outgoing_b.type();
  const ParticleType &outgoing_c_type = outgoing_c.type();

  log.debug("Note: Doing 1->3 decay!");

  const double mass_a = outgoing_a_type.mass();
  const double mass_b = outgoing_b_type.mass();
  const double mass_c = outgoing_c_type.mass();

  // std::cout << mass_a << " " << mass_b << std::endl;
  const double mass_resonance = incoming_particles_[0].effective_mass();
  //  std::cout << mass_resonance << std::endl;
  // mandelstam-s limits for pairs ab and bc
  const double s_ab_max = (mass_resonance - mass_c) * (mass_resonance - mass_c);
  const double s_ab_min = (mass_a + mass_b) * (mass_a + mass_b);
  const double s_bc_max = (mass_resonance - mass_a) * (mass_resonance - mass_a);
  const double s_bc_min = (mass_b + mass_c) * (mass_b + mass_c);

  log.debug("s_ab limits: ", s_ab_min, " ", s_ab_max);
  log.debug("s_bc limits: ", s_bc_min, " ", s_bc_max);

  /* randomly pick values for s_ab and s_bc
   * until the pair is within the Dalitz plot */
  double dalitz_bc_max = 0.0, dalitz_bc_min = 1.0;
  double s_ab = 0.0, s_bc = 0.5;
  while (s_bc > dalitz_bc_max || s_bc < dalitz_bc_min) {
    s_ab = random::uniform(s_ab_min, s_ab_max);
    s_bc = random::uniform(s_bc_min, s_bc_max);
    const double e_b_rest =
        (s_ab - mass_a * mass_a + mass_b * mass_b) / (2 * std::sqrt(s_ab));
    const double e_c_rest =
        (mass_resonance * mass_resonance - s_ab - mass_c * mass_c) /
        (2 * std::sqrt(s_ab));
    dalitz_bc_max = (e_b_rest + e_c_rest) * (e_b_rest + e_c_rest) -
                    (std::sqrt(e_b_rest * e_b_rest - mass_b * mass_b) -
                     std::sqrt(e_c_rest * e_c_rest - mass_c * mass_c)) *
                        (std::sqrt(e_b_rest * e_b_rest - mass_b * mass_b) -
                         std::sqrt(e_c_rest * e_c_rest - mass_c * mass_c));
    dalitz_bc_min = (e_b_rest + e_c_rest) * (e_b_rest + e_c_rest) -
                    (std::sqrt(e_b_rest * e_b_rest - mass_b * mass_b) +
                     std::sqrt(e_c_rest * e_c_rest - mass_c * mass_c)) *
                        (std::sqrt(e_b_rest * e_b_rest - mass_b * mass_b) +
                         std::sqrt(e_c_rest * e_c_rest - mass_c * mass_c));
  }

  log.debug("s_ab: ", s_ab, " s_bc: ", s_bc, " min: ", dalitz_bc_min,
            " max: ", dalitz_bc_max);

  // std::cout << "s_ab: " << s_ab << " s_bc " << s_bc << std::endl;

  // Compute energy and momentum magnitude
  const double energy_a =
      (mass_resonance * mass_resonance + mass_a * mass_a - s_bc) /
      (2 * mass_resonance);
  const double energy_c =
      (mass_resonance * mass_resonance + mass_c * mass_c - s_ab) /
      (2 * mass_resonance);
  const double energy_b =
      (s_ab + s_bc - mass_a * mass_a - mass_c * mass_c) / (2 * mass_resonance);
  const double momentum_a = std::sqrt(energy_a * energy_a - mass_a * mass_a);
  const double momentum_c = std::sqrt(energy_c * energy_c - mass_c * mass_c);
  const double momentum_b = std::sqrt(energy_b * energy_b - mass_b * mass_b);

  // CM system. "Incoming particle" is at rest.
  // const double total_energy = incoming_particles_[0].effective_mass();
  const double total_energy = incoming_particles_[0].momentum().abs();
  std::cout << "sqrts \t" << total_energy << std::endl;
  if (std::abs(energy_a + energy_b + energy_c - total_energy) > really_small) {
    log.warn("1->3: Ea + Eb + Ec: ", energy_a + energy_b + energy_c,
             " Total E: ", total_energy);
  }
  log.debug("Calculating the angles...");

  // momentum_a direction is random
  Angles phitheta;
  phitheta.distribute_isotropically();
  // This is the angle of the plane of the three decay particles
  outgoing_a.set_4momentum(mass_a, phitheta.threevec() * momentum_a);

  // Angle between a and b
  double theta_ab = std::acos(
      (energy_a * energy_b - 0.5 * (s_ab - mass_a * mass_a - mass_b * mass_b)) /
      (momentum_a * momentum_b));
  log.debug("theta_ab: ", theta_ab, " Ea: ", energy_a, " Eb: ", energy_b,
            " sab: ", s_ab, " pa: ", momentum_a, " pb: ", momentum_b);
  bool phi_has_changed = phitheta.add_to_theta(theta_ab);
  outgoing_b.set_4momentum(mass_b, phitheta.threevec() * momentum_b);

  // Angle between b and c
  double theta_bc = std::acos(
      (energy_b * energy_c - 0.5 * (s_bc - mass_b * mass_b - mass_c * mass_c)) /
      (momentum_b * momentum_c));
  log.debug("theta_bc: ", theta_bc, " Eb: ", energy_b, " Ec: ", energy_c,
            " sbc: ", s_bc, " pb: ", momentum_b, " pc: ", momentum_c);
  /* pass information on whether phi has changed during the last adding
   * on to add_to_theta: */
  phitheta.add_to_theta(theta_bc, phi_has_changed);
  outgoing_c.set_4momentum(mass_c, phitheta.threevec() * momentum_c);

  // Momentum check
  FourVector ptot =
      outgoing_a.momentum() + outgoing_b.momentum() + outgoing_c.momentum();

  if (std::abs(ptot.x0() - total_energy) > really_small) {
    log.warn("1->3 energy not conserved! Before: ", total_energy,
             " After: ", ptot.x0());
  }
  if (std::abs(ptot.x1()) > really_small ||
      std::abs(ptot.x2()) > really_small ||
      std::abs(ptot.x3()) > really_small) {
    log.warn("1->3 momentum check failed. Total momentum: ", ptot.threevec());
  }

  log.debug("outgoing_a: ", outgoing_a.momentum(),
            "\noutgoing_b: ", outgoing_b.momentum(),
            "\noutgoing_c: ", outgoing_c.momentum());

  // sanity check that we are really working with references.
  //  assert(outgoing_a.momentum[1] == outgoing_particles_[0].momentum[0]);
  return {outgoing_a, outgoing_b, outgoing_c};
}

TEST(init_particle_types) { Test::create_actual_particletypes(); }

/*
TEST(three_body_dima) {
  ParticleData kaon{ParticleType::find(0x311)};
  ParticleData pi1{ParticleType::find(0x111)};
  ParticleData pi2{ParticleType::find(0x211)};
  ParticleData pi3{ParticleType::find(0x211)};
  kaon.set_4momentum(kaon.type().mass() + 0.0, ThreeVector(0.0, 0.0, 0.0));
  ParticleList incoming{kaon};
  ParticleList outgoing{pi1, pi2, pi3};

  constexpr int N = 10000;

  std::fstream fs;
  fs.open("/home/rothermel/Work/three_body/dalitz_results/dima.txt",
          std::fstream::out);
  for (int i = 0; i < N; i++) {
    sample_3body_phase_space_dima(kaon.momentum().abs(), pi1, pi2, pi3);

    FourVector s12_temp = pi1.momentum() + pi2.momentum();
    FourVector s23_temp = pi2.momentum() + pi3.momentum();

    const double s12 = s12_temp.sqr();
    const double s23 = s23_temp.sqr();
    // std::cout << "s12: " << s12 << "\t s23: " << s23 << std::endl;
    // std::cout << s12 << std::endl;
    // std::cout << s23 << std::endl;

    fs << s12 << '\t' << s23 << '\t'
       << costheta_12(pi1.momentum(), pi2.momentum()) << '\t'
       << costheta_12(pi2.momentum(), pi3.momentum()) << '\t' << '\n';
  }

  fs.close();
}
 */
/*
TEST(three_body_decay) {
  // K -> Pi+ Pi- Pi0
  ParticleData kaon{ParticleType::find(0x311)};
  ParticleData pi1{ParticleType::find(0x111)};
  ParticleData pi2{ParticleType::find(0x211)};
  ParticleData pi3{ParticleType::find(0x211)};
  kaon.set_4momentum(kaon.type().mass() + 0.0, ThreeVector(0.0, 0.0, 0.0));
  ParticleList incoming{kaon};
  ParticleList outgoing{pi1, pi2, pi3};

  constexpr int N = 10000;

  std::fstream fs;
  fs.open("/home/rothermel/Work/three_body/dalitz_results/v1.txt",
          std::fstream::out);
  for (int i = 0; i < N; i++) {
    ParticleList products = one_to_three(incoming, outgoing);

    FourVector s12_temp = products[0].momentum() + products[1].momentum();
    FourVector s23_temp = products[1].momentum() + products[2].momentum();

    const double s12 = s12_temp.sqr();
    const double s23 = s23_temp.sqr();
    // std::cout << "s12: " << s12 << "\t s23: " << s23 << std::endl;
    std::cout << s12 << std::endl;
    std::cout << s23 << std::endl;

    fs << s12 << '\t' << s23 << '\t'
       << costheta_12(products[0].momentum(), products[1].momentum()) << '\t'
       << costheta_12(products[1].momentum(), products[2].momentum()) << '\t'
       << '\n';
  }

  fs.close();
}
 */

ThreeVector random_vector() {
  const double x = random::uniform(0.0, 5.0);
  const double y = random::uniform(0.0, 5.0);
  const double z = random::uniform(0.0, 5.0);
  return ThreeVector(x, y, z);
}
TEST(three_body_jonas) {
  ParticleData kaon{ParticleType::find(0x311)};
  ParticleData pi1{ParticleType::find(0x111)};
  ParticleData pi2{ParticleType::find(0x211)};
  ParticleData pi3{ParticleType::find(0x211)};
  ThreeVector vec = random_vector();
  const double energy =
      std::sqrt(std::pow(kaon.type().mass(), 2) + std::pow(vec.abs(), 2));
  kaon.set_4momentum(energy, vec);
  ParticleList incoming{kaon};
  ParticleList outgoing{pi1, pi2, pi3};

  constexpr int N = 1;

  std::fstream fs;
  fs.open("/home/rothermel/Work/three_body/dalitz_results/v_jonas.txt",
          std::fstream::out);
  for (int i = 0; i < N; i++) {
    ParticleList products = sample_3body_phase_space_jonas(kaon, outgoing);

    FourVector s12_temp = products[0].momentum() + products[1].momentum();
    FourVector s23_temp = products[1].momentum() + products[2].momentum();

    const double s12 = s12_temp.sqr();
    const double s23 = s23_temp.sqr();
    // std::cout << "s12: " << s12 << "\t s23: " << s23 << std::endl;
    // std::cout << s12 << std::endl;
    // std::cout << s23 << std::endl;

    fs << s12 << '\t' << s23 << '\t';
    for (auto &out : products) {
      fs << out.momentum().threevec().get_phi() << '\t'
         << out.momentum().threevec().get_theta() << '\t';
    }
    fs << '\n';
  }

  fs.close();
}
