/*
 *
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */
#include "include/resonances.h"

#include <gsl/gsl_sf_coupling.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>

#include "include/constants.h"
#include "include/DecayModes.h"
#include "include/distributions.h"
#include "include/FourVector.h"
#include "include/macros.h"
#include "include/outputroutines.h"
#include "include/Particles.h"
#include "include/ParticleData.h"
#include "include/ParticleType.h"

/* calculate_minimum_mass
 * - calculate the minimum rest energy the resonance must have
 * to be able to decay through any of its decay channels
 * NB: This function assumes stable decay products!
 */
float calculate_minimum_mass(Particles *particles, int pdgcode) {
  /* If the particle happens to be stable, just return the mass */
  if ((particles->particle_type(pdgcode)).width() < 0.0)
    return (particles->particle_type(pdgcode)).mass();
  /* Otherwise, let's find the highest mass value needed in any decay mode */
  float minimum_mass = 0.0;
  const std::vector< std::pair<std::vector<int>, float> > decaymodes
    = (particles->decay_modes(pdgcode)).decay_mode_list();
  for (std::vector< std::pair<std::vector<int>, float> >::const_iterator
         mode = decaymodes.begin(); mode != decaymodes.end(); ++mode) {
    size_t decay_particles = (mode->first).size();
    float total_mass = 0.0;
    for (size_t i = 0; i < decay_particles; i++) {
      /* Stable decay products assumed; for resonances the mass can be lower! */
      total_mass += particles->particle_type((mode->first)[i]).mass();
    }
    if (total_mass > minimum_mass)
      minimum_mass = total_mass;
  }
  return minimum_mass;
}

/* resonance_cross_section - energy-dependent cross section
 * for producing a resonance
 */
std::map<int, double> resonance_cross_section(
  const ParticleData &particle1, const ParticleData &particle2,
  const ParticleType &type_particle1, const ParticleType &type_particle2,
  Particles *particles) {
  const int charge1 = type_particle1.charge(),
    charge2 = type_particle2.charge();
  const int pdgcode1 = type_particle1.pdgcode(),
    pdgcode2 = type_particle2.pdgcode();
  /* Isospin z-component based on Gell-Mann–Nishijima formula
   * 2 * Iz = 2 * charge - (baryon number + strangeness + charm)
   * XXX: Strangeness and charm ignored for now!
   */
  const int isospin_z1 = type_particle1.spin() % 2 == 0
                         ? charge1 * 2
                         : charge1 * 2 - pdgcode1 / abs(pdgcode1);
  const int isospin_z2 = type_particle2.spin() % 2 == 0
                         ? charge2 * 2
                         : charge2 * 2 - pdgcode2 / abs(pdgcode2);
  std::map<int, double> possible_resonances;

  /* key 0 refers to total resonance production cross section */
  possible_resonances[0] = 0.0;

  /* Resonances do not form resonances */
  if (type_particle1.width() > 0.0 || type_particle2.width() > 0.0)
    return possible_resonances;

  /* No baryon-baryon interactions for now */
  if (type_particle1.spin() % 2 != 0 && type_particle2.spin() % 2 != 0)
    return possible_resonances;

  /* Isospin symmetry factor, by default 1 */
  int symmetryfactor = 1;
  /* Do the particles have the same isospin value? */
  if (type_particle1.isospin() == type_particle2.isospin()) {
    /* Do they have the same spin? */
    if (type_particle1.spin() == type_particle2.spin()) {
      /* Are their PDG codes of same length? */
      int abs_pdg1 = abs(pdgcode1), digits1 = 0;
      while (abs_pdg1) {
        abs_pdg1 /= 10;
        digits1++;
      }
      int abs_pdg2 = abs(pdgcode2), digits2 = 0;
      while (abs_pdg2) {
        abs_pdg2 /= 10;
        digits2++;
      }
      if (digits1 == digits2) {
        /* If baryons, do they have the same baryon number? */
        if (type_particle1.spin() % 2 == 0 ||
            std::signbit(pdgcode1) == std::signbit(pdgcode2)) {
          /* Ok, particles are in the same isospin multiplet,
             apply symmetry factor */
          symmetryfactor = 2;
        }
      }
    }
  }

  /* Mandelstam s = (p_a + p_b)^2 = square of CMS energy */
  const double mandelstam_s =
       (particle1.momentum() + particle2.momentum()).Dot(
         particle1.momentum() + particle2.momentum());

  /* CM momentum */
  const double cm_momentum_squared
    = (particle1.momentum().Dot(particle2.momentum())
       * particle1.momentum().Dot(particle2.momentum())
       - type_particle1.mass() * type_particle1.mass()
       * type_particle2.mass() * type_particle2.mass()) / mandelstam_s;

  /* Find all the possible resonances */
  for (std::map<int, ParticleType>::const_iterator
       i = particles->types_cbegin(); i != particles->types_cend(); ++i) {
       ParticleType type_resonance = i->second;
    /* Not a resonance, go to next type of particle */
    if (type_resonance.width() < 0.0)
      continue;

    /* Check for charge conservation */
    if (type_resonance.charge() != charge1 + charge2)
      continue;

    /* Check for baryon number conservation */
    if (type_particle1.spin() % 2 != 0 || type_particle2.spin() % 2 != 0) {
      /* Step 1: We must have fermion */
      if (type_resonance.spin() % 2 == 0) {
        continue;
      }
      /* Step 2: We must have antiparticle for antibaryon
       * (and non-antiparticle for baryon)
       */
      if (type_particle1.spin() % 2 != 0
          && (std::signbit(pdgcode1)
              != std::signbit(type_resonance.pdgcode()))) {
        continue;
      } else if (type_particle2.spin() % 2 != 0
          && (std::signbit(pdgcode2)
          != std::signbit(type_resonance.pdgcode()))) {
        continue;
      }
    }

    int isospin_z_resonance = (type_resonance.spin()) % 2 == 0
     ? type_resonance.charge() * 2
     : type_resonance.charge() * 2 - type_resonance.pdgcode()
                                    / abs(type_resonance.pdgcode());

    /* Calculate isospin Clebsch-Gordan coefficient
     * (-1)^(j1 - j2 + m3) * sqrt(2 * j3 + 1) * [Wigner 3J symbol]
     * Note that the calculation assumes that isospin values
     * have been multiplied by two
     */
    double wigner_3j =  gsl_sf_coupling_3j(type_particle1.isospin(),
       type_particle2.isospin(), type_resonance.isospin(),
       isospin_z1, isospin_z2, -isospin_z_resonance);
    double clebsch_gordan_isospin = 0.0;
    if (fabs(wigner_3j) > really_small)
      clebsch_gordan_isospin = pow(-1, type_particle1.isospin() / 2.0
      - type_particle2.isospin() / 2.0 + isospin_z_resonance / 2.0)
      * sqrt(type_resonance.isospin() + 1) * wigner_3j;

    printd("CG: %g I1: %i I2: %i IR: %i iz1: %i iz2: %i izR: %i \n",
         clebsch_gordan_isospin,
         type_particle1.isospin(), type_particle2.isospin(),
         type_resonance.isospin(),
         isospin_z1, isospin_z2, isospin_z_resonance);

    /* If Clebsch-Gordan coefficient is zero, don't bother with the rest */
    if (fabs(clebsch_gordan_isospin) < really_small)
      continue;

    /* Check the decay modes of this resonance */
    const std::vector< std::pair<std::vector<int>, float> > decaymodes
      = (particles->decay_modes(type_resonance.pdgcode())).decay_mode_list();
    bool not_enough_energy = false;
    for (std::vector< std::pair<std::vector<int>, float> >::const_iterator mode
         = decaymodes.begin(); mode != decaymodes.end(); ++mode) {
      size_t decay_particles = (mode->first).size();
      if ( decay_particles > 3 ) {
        printf("Warning: Not a 1->2 or 1->3 process!\n");
        printf("Number of decay particles: %zu \n", decay_particles);
      } else {
        /* There must be enough energy to produce all decay products */
        float mass_a, mass_b, mass_c = 0.0;
        mass_a = calculate_minimum_mass(particles, (mode->first)[0]);
        mass_b = calculate_minimum_mass(particles, (mode->first)[1]);
        if (decay_particles == 3) {
          mass_c = calculate_minimum_mass(particles, (mode->first)[2]);
        }
        if (sqrt(mandelstam_s) < mass_a + mass_b + mass_c)
          not_enough_energy = true;
      }
    }
    if (not_enough_energy)
      continue;

    /* Calculate spin factor */
    const double spinfactor = (type_resonance.spin() + 1)
      / ((type_particle1.spin() + 1) * (type_particle2.spin() + 1));

    float resonance_width = type_resonance.width();
    float resonance_mass = type_resonance.mass();
    /* Calculate resonance production cross section
     * using the Breit-Wigner distribution as probability amplitude
     */
    double resonance_xsection =  clebsch_gordan_isospin * clebsch_gordan_isospin
         * spinfactor * symmetryfactor
         * 4.0 * M_PI / cm_momentum_squared
         * breit_wigner(mandelstam_s, resonance_mass, resonance_width)
         * hbarc * hbarc / fm2_mb;

    /* If cross section is non-negligible, add resonance to the list */
    if (resonance_xsection > really_small) {
      possible_resonances[type_resonance.pdgcode()] = resonance_xsection;
      possible_resonances[0] += resonance_xsection;
      printd("Found resonance %i (%s) with mass %f and width %f.\n",
             type_resonance.pdgcode(), type_resonance.name().c_str(),
             resonance_mass, resonance_width);
      printd("Original particles: %s %s Charges: %i %i \n",
             type_particle1.name().c_str(), type_particle2.name().c_str(),
             type_particle1.charge(), type_particle2.charge());
    }
  }
  return possible_resonances;
}

/* 1->2 process kinematics */
int one_to_two(Particles *particles, int resonance_id, int type_a, int type_b) {
  /* Add two new particles */
  ParticleData new_particle_a, new_particle_b;
  new_particle_a.set_pdgcode(type_a);
  new_particle_b.set_pdgcode(type_b);

  double mass_a = particles->particle_type(type_a).mass(),
    mass_b = particles->particle_type(type_b).mass();
  const double total_energy = particles->data(resonance_id).momentum().x0();
  double energy_a = (total_energy * total_energy
                       + mass_a * mass_a - mass_b * mass_b)
    / (2.0 * total_energy);

  double momentum_radial = sqrt(energy_a * energy_a - mass_a * mass_a);
  if (momentum_radial < 0.0)
    printf("Warning: radial momenta %g \n", momentum_radial);
  /* phi in the range from [0, 2 * pi) */
  double phi = 2.0 * M_PI * drand48();
  /* cos(theta) in the range from [-1.0, 1.0) */
  double cos_theta = -1.0 + 2.0 * drand48();
  double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
  if (energy_a  < mass_a || abs(cos_theta) > 1) {
    printf("Particle %d radial momenta %g phi %g cos_theta %g\n", type_a,
           momentum_radial, phi, cos_theta);
    printf("Etot: %g m_a: %g m_b %g E_a: %g", total_energy, mass_a, mass_b,
           energy_a);
  }
  new_particle_a.set_momentum(mass_a,
                              momentum_radial * cos(phi) * sin_theta,
                              momentum_radial * sin(phi) * sin_theta,
                              momentum_radial * cos_theta);
  new_particle_b.set_momentum(mass_b,
                              - new_particle_a.momentum().x1(),
                              - new_particle_a.momentum().x2(),
                              - new_particle_a.momentum().x3());

  /* Both decay products begin from the same point */
  FourVector decay_point = particles->data(resonance_id).position();
  new_particle_a.set_position(decay_point);
  new_particle_b.set_position(decay_point);

  /* No collision yet */
  new_particle_a.set_collision(-1, 0, -1);
  new_particle_b.set_collision(-1, 0, -1);

  /* Assign IDs to new particles */
  int new_id_a = particles->id_max() + 1;
  int new_id_b = new_id_a + 1;
  new_particle_a.set_id(new_id_a);
  new_particle_b.set_id(new_id_b);

  particles->add_data(new_particle_a);
  particles->add_data(new_particle_b);

  printd("Created %s and %s with IDs %d and %d \n",
         particles->type(new_id_a).name().c_str(),
         particles->type(new_id_b).name().c_str(), new_id_a, new_id_b);

  return new_id_a;
}

/* 1->3 process kinematics */
int one_to_three(Particles *particles, int resonance_id,
                 int type_a, int type_b, int type_c) {
  /* Add three new particles */
  ParticleData new_particle_a, new_particle_b, new_particle_c;
  new_particle_a.set_pdgcode(type_a);
  new_particle_b.set_pdgcode(type_b);
  new_particle_c.set_pdgcode(type_c);

  FourVector momentum_resonance = particles->data(resonance_id).momentum();
  const double mass_a = particles->particle_type(type_a).mass(),
    mass_b = particles->particle_type(type_b).mass(),
    mass_c = particles->particle_type(type_c).mass(),
    mass_resonance = sqrt(momentum_resonance.Dot(momentum_resonance));

  /* mandelstam-s limits for pairs ab and bc */
  double s_ab_max = (mass_resonance - mass_c) * (mass_resonance - mass_c);
  double s_ab_min = (mass_a + mass_b) * (mass_a + mass_b);
  double s_bc_max = (mass_resonance - mass_a) * (mass_resonance - mass_a);
  double s_bc_min = (mass_b + mass_c) * (mass_b + mass_c);

  printd("s_ab limits: %g %g \n", s_ab_min, s_ab_max);
  printd("s_bc limits: %g %g \n", s_bc_min, s_bc_max);

  /* randomly pick values for s_ab and s_bc
   * until the pair is within the Dalitz plot */
  double dalitz_bc_max = 0.0, dalitz_bc_min = 1.0;
  double s_ab = 0.0, s_bc = 0.5;
  while (s_bc > dalitz_bc_max || s_bc < dalitz_bc_min) {
    s_ab = (s_ab_max - s_ab_min) * drand48() + s_ab_min;
    s_bc = (s_bc_max - s_bc_min) * drand48() + s_bc_min;
    double e_b_rest = (s_ab - mass_a * mass_a + mass_b * mass_b)
                           / (2 * sqrt(s_ab));
    double e_c_rest = (mass_resonance * mass_resonance - s_ab
                            - mass_c * mass_c) / (2 * sqrt(s_ab));
    dalitz_bc_max = (e_b_rest + e_c_rest) * (e_b_rest + e_c_rest)
      - (sqrt(e_b_rest * e_b_rest - mass_b * mass_b)
         - sqrt(e_c_rest * e_c_rest - mass_c * mass_c))
      * (sqrt(e_b_rest * e_b_rest - mass_b * mass_b)
         - sqrt(e_c_rest * e_c_rest - mass_c * mass_c));
    dalitz_bc_min = (e_b_rest + e_c_rest) * (e_b_rest + e_c_rest)
      - (sqrt(e_b_rest * e_b_rest - mass_b * mass_b)
         + sqrt(e_c_rest * e_c_rest - mass_c * mass_c))
      * (sqrt(e_b_rest * e_b_rest - mass_b * mass_b)
         + sqrt(e_c_rest * e_c_rest - mass_c * mass_c));
  }

  printd("s_ab: %g s_bc: %g min: %g max: %g\n",
         s_ab, s_bc, dalitz_bc_min, dalitz_bc_max);

  /* Compute energy and momentum magnitude */
  const double energy_a = (mass_resonance * mass_resonance + mass_a * mass_a
                           - s_bc) / (2 * mass_resonance);
  const double energy_c = (mass_resonance * mass_resonance + mass_c * mass_c
                           - s_ab) / (2 * mass_resonance);
  const double energy_b = (s_ab + s_bc - mass_a * mass_a - mass_c * mass_c)
                           / (2 * mass_resonance);
  const double momentum_a = sqrt(energy_a * energy_a - mass_a * mass_a);
  const double momentum_c = sqrt(energy_c * energy_c - mass_c * mass_c);
  const double momentum_b = sqrt(energy_b * energy_b - mass_b * mass_b);

  const double total_energy = momentum_resonance.x0();
  if (fabs(energy_a + energy_b + energy_c - total_energy) > really_small)
    printf("1->3 warning: Ea + Eb + Ec: %g Total E: %g",
           energy_a + energy_b + energy_c, total_energy);
  printd("Calculating the angles...\n");

  /* momentum_a direction is random */
  /* phi in the range from [0, 2 * pi) */
  /* This is the angle of the plane of the three decay particles */
  double phi = 2.0 * M_PI * drand48();
  /* cos(theta) in the range from [-1.0, 1.0) */
  double cos_theta_a = -1.0 + 2.0 * drand48();
  double sin_theta_a = sqrt(1.0 - cos_theta_a * cos_theta_a);
  new_particle_a.set_momentum(mass_a,
                              momentum_a * cos(phi) * sin_theta_a,
                              momentum_a * sin(phi) * sin_theta_a,
                              momentum_a * cos_theta_a);

  /* Angle between a and b */
  double theta_ab = acos((energy_a * energy_b - 0.5 * (s_ab - mass_a * mass_a
                          - mass_b * mass_b)) / (momentum_a * momentum_b));
  printd("theta_ab: %g Ea: %g Eb: %g sab: %g pa: %g pb: %g\n",
         theta_ab, energy_a, energy_b, s_ab, momentum_a, momentum_b);
  /* b angle is sum of a angle and ab angle */
  double theta_b = theta_ab + acos(cos_theta_a);
  new_particle_b.set_momentum(mass_b,
                              momentum_b * cos(phi) * sin(theta_b),
                              momentum_b * sin(phi) * sin(theta_b),
                              momentum_b * cos(theta_b));

  /* Angle between b and c */
  double theta_bc = acos((energy_b * energy_c - 0.5 *(s_bc - mass_b * mass_b
                         - mass_c * mass_c)) / (momentum_b * momentum_c));
  printd("theta_bc: %g Eb: %g Ec: %g sbc: %g pb: %g pc: %g\n",
         theta_bc, energy_b, energy_c, s_bc, momentum_b, momentum_c);
  /* c angle is sum of b angle and bc angle */
  double theta_c = theta_bc + theta_b;
  new_particle_c.set_momentum(mass_c,
                              momentum_c * cos(phi) * sin(theta_c),
                              momentum_c * sin(phi) * sin(theta_c),
                              momentum_c * cos(theta_c));

  /* Momentum check */
  double energy = new_particle_a.momentum().x0()
    + new_particle_b.momentum().x0() + new_particle_c.momentum().x0();
  double px = new_particle_a.momentum().x1() + new_particle_b.momentum().x1()
    + new_particle_c.momentum().x1();
  double py = new_particle_a.momentum().x2() + new_particle_b.momentum().x2()
    + new_particle_c.momentum().x2();
  double pz = new_particle_a.momentum().x3() + new_particle_b.momentum().x3()
    + new_particle_c.momentum().x3();

  if (fabs(energy - total_energy) > really_small)
    printf("1->3 energy not conserved! Before: %g After: %g",
           total_energy, energy);

  if (px > really_small || py > really_small || pz > really_small)
    printf("1->3 momentum check failed. Total momentum: %g %g %g\n",
           px, py, pz);

  /* All decay products begin from the same point */
  FourVector decay_point = particles->data(resonance_id).position();
  new_particle_a.set_position(decay_point);
  new_particle_b.set_position(decay_point);
  new_particle_c.set_position(decay_point);

  /* No collision yet */
  new_particle_a.set_collision(-1, 0, -1);
  new_particle_b.set_collision(-1, 0, -1);
  new_particle_c.set_collision(-1, 0, -1);

  /* Assign IDs to new particles */
  int new_id_a = particles->id_max() + 1;
  int new_id_b = new_id_a + 1;
  int new_id_c = new_id_b + 1;
  new_particle_a.set_id(new_id_a);
  new_particle_b.set_id(new_id_b);
  new_particle_c.set_id(new_id_c);

  particles->add_data(new_particle_a);
  particles->add_data(new_particle_b);
  particles->add_data(new_particle_c);

  printd("Created %s %s %s with IDs %d %d %d \n",
         particles->type(new_id_a).name().c_str(),
         particles->type(new_id_b).name().c_str(),
         particles->type(new_id_c).name().c_str(),
         new_id_a, new_id_b, new_id_c);

  printd("p0: %g %g %g \n", new_particle_a.momentum().x0(),
         new_particle_b.momentum().x0(), new_particle_c.momentum().x0());
  printd("p1: %g %g %g \n", new_particle_a.momentum().x1(),
         new_particle_b.momentum().x1(), new_particle_c.momentum().x1());
  printd("p2: %g %g %g \n", new_particle_a.momentum().x2(),
         new_particle_b.momentum().x2(), new_particle_c.momentum().x2());
  printd("p3: %g %g %g \n", new_particle_a.momentum().x3(),
         new_particle_b.momentum().x3(), new_particle_c.momentum().x3());

  return new_id_a;
}

/* Resonance decay process */
int resonance_decay(Particles *particles, int particle_id) {
  const int pdgcode = particles->type(particle_id).pdgcode();
  const std::vector< std::pair<std::vector<int>, float> > decaymodes
    = (particles->decay_modes(pdgcode)).decay_mode_list();
  int type_a = 0, type_b = 0, type_c = 0, new_id_a = -1;

  /* Ratios of decay channels should add to 1; pick a random number
   * between 0 and 1 to select the decay mode to be used
   */
  double random_mode = drand48();
  double cumulated_probability = 0.0;
  size_t decay_particles = 0;
  for (std::vector< std::pair<std::vector<int>, float> >::const_iterator mode
         = decaymodes.begin(); mode != decaymodes.end(); ++mode) {
    cumulated_probability += mode->second;
    if (random_mode < cumulated_probability) {
      decay_particles = (mode->first).size();
      if ( decay_particles > 3 ) {
        printf("Warning: Not a 1->2 or 1->3 process!\n");
        printf("Number of decay particles: %zu \n", decay_particles);
        printf("Decay particles: ");
        for (size_t i = 0; i < decay_particles; i++) {
          printf("%i ", (mode->first)[i]);
        }
        printf("\n");
      } else if (decay_particles == 2) {
        type_a = (mode->first)[0];
        type_b = (mode->first)[1];
        if (abs(type_a) < 100 || abs(type_b) < 100)
          printf("Warning: decay products A: %i B: %i\n", type_a, type_b);
      } else if (decay_particles == 3) {
        type_a = (mode->first)[0];
        type_b = (mode->first)[1];
        type_c = (mode->first)[2];
        if (abs(type_a) < 100 || abs(type_b) < 100 || abs(type_c) < 100)
          printf("Warning: decay products A: %i B: %i C: %i\n",
                 type_a, type_b, type_c);
      }
    }
  }
  if (decay_particles == 2) {
    new_id_a = one_to_two(particles, particle_id, type_a, type_b);
  } else if (decay_particles == 3) {
    printd("Note: Doing 1->3 decay!\n");
    new_id_a = one_to_three(particles, particle_id, type_a, type_b, type_c);
  }
  return new_id_a;
}

/* 2->1 resonance formation process */
int resonance_formation(Particles *particles, int particle_id, int other_id,
  int pdg_resonance) {
  /* Add a new particle */
  ParticleData resonance;
  resonance.set_pdgcode(pdg_resonance);

  /* Center-of-momentum frame of initial particles
   * is the rest frame of the resonance
   */
  const double energy = particles->data(particle_id).momentum().x0()
    + particles->data(other_id).momentum().x0();
  /* We use fourvector to set 4-momentum, as setting it
   * with doubles requires that particle is on
   * mass shell, which is not generally true for resonances
   */
  FourVector resonance_momentum(energy, 0.0, 0.0, 0.0);
  resonance.set_momentum(resonance_momentum);

  printd("Momentum of the new particle: %g %g %g %g \n",
    resonance.momentum().x0(),
    resonance.momentum().x1(),
    resonance.momentum().x2(),
    resonance.momentum().x3());

  /* The real position should be between parents in the computational frame! */
  resonance.set_position(1.0, 0.0, 0.0, 0.0);

  /* No collision yet */
  resonance.set_collision(-1, 0, -1);
  int new_id = particles->id_max() + 1;
  resonance.set_id(new_id);
  particles->add_data(resonance);
  printd("Created %s with ID %i \n", particles->type(new_id).name().c_str(),
         new_id);

  return new_id;
}
