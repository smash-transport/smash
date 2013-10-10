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

/* resonance_cross_section - energy-dependent cross section
 * for producing a resonance
 */
std::map<int, double> resonance_cross_section(
  const ParticleData &particle1, const ParticleData &particle2,
  const ParticleType &type_particle1, const ParticleType &type_particle2,
  const Particles &particles) {
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
       i = particles.types_cbegin(); i != particles.types_cend(); ++i) {
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

/* 1->2 resonance decay process */
int resonance_decay(Particles *particles, int particle_id) {
  const int pdgcode = particles->type(particle_id).pdgcode();
  const std::vector< std::pair<std::vector<int>, float> > decaymodes
    = (particles->decay_modes(pdgcode)).decay_mode_list();
  const double total_energy = particles->data(particle_id).momentum().x0();
  int type_a = 0, type_b = 0;

  /* Ratios of decay channels should add to 1; pick a random number
   * between 0 and 1 to select the decay mode to be used
   */
  double random_mode = drand48();
  double cumulated_probability = 0.0;
  for (std::vector< std::pair<std::vector<int>, float> >::const_iterator mode
         = decaymodes.begin(); mode != decaymodes.end(); ++mode) {
    cumulated_probability += mode->second;
    if (random_mode < cumulated_probability) {
      if ( (mode->first).size() != 2 ) {
        printf("Warning: Not a 1->2 process! Number of decay particles: %zu \n",
               (mode->first).size());
        printf("Decay particles: ");
        for (size_t i = 0; i < (mode->first).size(); i++) {
          printf("%i ", (mode->first)[i]);
        }
        printf("\n");
      } else {
        type_a = (mode->first)[0];
        type_b = (mode->first)[1];
        if (abs(type_a) < 100 || abs(type_b) < 100)
          printf("Warning: decay products A: %i B: %i\n", type_a, type_b);
     }
    }
  }

  if (abs(type_a) < 100 || abs(type_b) < 100)
    printf("Warning: decay products A: %i B: %i\n", type_a, type_b);
  /* Add two new particles */
  ParticleData new_particle_a, new_particle_b;
  new_particle_a.set_pdgcode(type_a);
  new_particle_b.set_pdgcode(type_b);

  double mass_a = particles->particle_type(type_a).mass(),
    mass_b = particles->particle_type(type_b).mass();
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
  FourVector decay_point = particles->data(particle_id).position();
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
