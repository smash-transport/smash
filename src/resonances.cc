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
