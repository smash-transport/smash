#include "unittest.h"  // This include has to be first

#include "setup.h"

#include <cmath>
#include <fstream>
#include <typeinfo>

#include "../include/smash/cxx14compat.h"
#include "../include/smash/decayaction.h"
#include "../include/smash/decaymodes.h"

using namespace smash;

double costheta_12(FourVector p1, FourVector p2) {
  double tmp = p1.threevec() * p2.threevec();
  return tmp / (p1.abs3() * p2.abs3());
}

double kinematic_g_function(x, y, z, u, v, w) {
  // check this with mathematica

  return 2.0;
}

ParticleList sample_3body_phase_space_jonas(const double sqrts,
                                            const ParticleList &outgoing) {

  // general considerations: 
  // we work in the R23 frame, i.e. the rest frame of particles 2 and 3. 
  // we will use relations relating s1, s2 to the energy of the 
  // mother particle in this frame. from there we know how to boost in the 
  // lab system. 
  
  // sample s1 and s2
  // check that they are in bounds
  // sample Omega1 and phi3. Think about their significance
  // construct 4 vectors

  const double m1 = outgoing[0].type().mass();
  const double m2 = outgoing[1].type().mass();
  const double m3 = outgoing[2].type().mass();

  const double s1_min =
      pow(outgoing[0].type().mass() + outgoing[1].type().mass(), 2);
  const double s1_max = pow(sqrts - outgoing[2].type().mass(), 2);
  const double s2_min =
      pow(outgoing[1].type().mass() + outgoing[2].type().mass(), 2);
  const double s2_max = pow(sqrts - outgoing[0].type().mass(), 2);

  bool done = false;
  while (!done) {
    double s1 = random::uniform(s1_min, s1_max);
    double s2 = random::uniform(s2_min, s2_max);

    if (kinematic_g_function(s1, s2, sqrts * sqrts, m2 * m2, m1 * m1, m3 * m3) <
        0)
      done = true;
  }

  // from now on work in R23. here p2 = -p3, p1 = p

  Angles phi_theta1;
  phi_theta1.distribute_isotropically();
  2

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
  assert(outgoing_a.momentum[1] == outgoing_particles_[0].momentum[0]);
  return {outgoing_a, outgoing_b, outgoing_c};
}

TEST(init_particle_types) { Test::create_actual_particletypes(); }

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
