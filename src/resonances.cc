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
#include <cstring>
#include <map>
#include <vector>

#include "include/constants.h"
#include "include/distributions.h"
#include "include/FourVector.h"
#include "include/macros.h"
#include "include/outputroutines.h"
#include "include/ParticleData.h"
#include "include/ParticleType.h"

/* resonance_cross_section - energy-dependent cross section
 * for producing a resonance
 */
std::map<int, double> resonance_cross_section(
  ParticleData *particle1, ParticleData *particle2,
  ParticleType *type_particle1, ParticleType *type_particle2,
  std::vector<ParticleType> *type_list) {
  const int charge1 = (*type_particle1).charge(),
    charge2 = (*type_particle2).charge();

  /* Isospin z-component based on Gell-Mann–Nishijima formula
   * 2 * Iz = 2 * charge - (baryon number + strangeness + charm)
   * XXX: Strangeness and charm ignored for now!
   */
  const int isospin_z1 = type_particle1->spin() % 2 == 0
                         ? charge1 * 2
                         : charge1 * 2 - type_particle1->pdgcode()
                                       / abs(type_particle1->pdgcode());
  const int isospin_z2 = type_particle2->spin() % 2 == 0
                         ? charge2 * 2
                         : charge2 * 2 - type_particle2->pdgcode()
                                       / abs(type_particle2->pdgcode());
  std::map<int, double> possible_resonances;

  /* key 0 refers to total resonance production cross section */
  possible_resonances[0] = 0.0;

  /* Resonances do not form resonances */
  if (type_particle1->width() > 0.0 || type_particle2->width() > 0.0)
    return possible_resonances;

  /* No baryon-baryon interactions for now */
  if (type_particle1->spin() % 2 != 0 && type_particle2->spin() % 2 != 0)
    return possible_resonances;

  /* Symmetry factor If initial state particles are identical,
   *  multiply by two. */
  int symmetryfactor = 1;
  if (unlikely(type_particle1->pdgcode() == type_particle2->pdgcode()))
    symmetryfactor = 2;

  /* Mandelstam s = (p_a + p_b)^2 = square of CMS energy */
  const double mandelstam_s =
       ( (*particle1).momentum() + (*particle2).momentum()).Dot(
         (*particle1).momentum() + (*particle2).momentum() );

  /* CM momentum */
  const double cm_momentum_squared
    = (particle1->momentum().Dot(particle2->momentum())
       * particle1->momentum().Dot(particle2->momentum())
       - type_particle1->mass() * type_particle1->mass()
       * type_particle2->mass() * type_particle2->mass()) / mandelstam_s;

  /* Find all the possible resonances */
  for (std::vector<ParticleType>::iterator type_resonance = type_list->begin();
       type_resonance != type_list->end(); ++type_resonance) {
    /* Not a resonance, go to next type of particle */
    if (type_resonance->width() < 0.0)
      continue;

    /* Check for charge conservation */
    if (type_resonance->charge() != charge1 + charge2)
      continue;

    /* Check for baryon number conservation */
    if (type_particle1->spin() % 2 != 0 || type_particle2->spin() % 2 != 0) {
      /* Step 1: We must have fermion */
      if (type_resonance->spin() % 2 == 0) {
        continue;
      }
      /* Step 2: We must have antiparticle for antibaryon
       * (and non-antiparticle for baryon)
       */
      if (type_particle1->spin() % 2 != 0
          && !(std::signbit(type_particle1->pdgcode())
          && std::signbit(type_resonance->pdgcode()))) {
        continue;
      } else if (type_particle2->spin() % 2 != 0
          && !(std::signbit(type_particle2->pdgcode())
          && std::signbit(type_resonance->pdgcode()))) {
        continue;
      }
    }

    int isospin_z_resonance = (type_resonance->spin()) % 2 == 0
     ? type_resonance->charge() * 2
     : type_resonance->charge() * 2 - type_resonance->pdgcode()
                                    / abs(type_resonance->pdgcode());

    /* Calculate isospin Clebsch-Gordan coefficient
     * (-1)^(j1 - j2 + m3) * sqrt(2 * j3 + 1) * [Wigner 3J symbol]
     * Note that the calculation assumes that isospin values
     * have been multiplied by two
     */
    double wigner_3j =  gsl_sf_coupling_3j(type_particle1->isospin(),
       type_particle2->isospin(), type_resonance->isospin(),
       isospin_z1, isospin_z2, -isospin_z_resonance);
    double clebsch_gordan_isospin = 0.0;
    if (fabs(wigner_3j) > really_small)
      clebsch_gordan_isospin = pow(-1, type_particle1->isospin() / 2.0
      - type_particle2->isospin() / 2.0 + isospin_z_resonance / 2.0)
      * sqrt(type_resonance->isospin() + 1) * wigner_3j;

    printd("CG: %g I1: %i I2: %i IR: %i iz1: %i iz2: %i izR: %i \n",
         clebsch_gordan_isospin,
         type_particle1->isospin(), type_particle2->isospin(),
         type_resonance->isospin(),
         isospin_z1, isospin_z2, isospin_z_resonance);

    /* If Clebsch-Gordan coefficient is zero, don't bother with the rest */
    if (fabs(clebsch_gordan_isospin) < really_small)
      continue;

    /* Calculate spin factor */
    const double spinfactor = (type_resonance->spin() + 1)
      / ((type_particle1->spin() + 1) * (type_particle2->spin() + 1));

    float resonance_width = type_resonance->width();
    float resonance_mass = type_resonance->mass();
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
      possible_resonances[type_resonance->pdgcode()] = resonance_xsection;
      possible_resonances[0] += resonance_xsection;
      printd("Found resonance %i (%s) with mass %f and width %f.\n",
             type_resonance->pdgcode(), type_resonance->name().c_str(),
             resonance_mass, resonance_width);
      printd("Original particles: %s %s Charges: %i %i \n",
             type_particle1->name().c_str(), type_particle2->name().c_str(),
             type_particle1->charge(), type_particle2->charge());
    }
  }
  return possible_resonances;
}

/* 1->2 resonance decay process */
size_t resonance_decay(std::map<int, ParticleData> *particles,
  std::vector<ParticleType> *types, std::map<int, int> *map_type,
  int *particle_id, int *id_max) {
  /* Add two new particles */
  size_t new_id_a = *id_max + 1;
  (*id_max)++;
  size_t new_id_b = new_id_a + 1;
  (*id_max)++;
  {
  ParticleData new_particle_a, new_particle_b;
  (*particles)[new_id_a] = new_particle_a;
  (*particles)[new_id_a].set_id(new_id_a);
  (*particles)[new_id_b] = new_particle_b;
  (*particles)[new_id_b].set_id(new_id_b);
  }

  const int charge = (*types)[(*map_type)[*particle_id]].charge();
  int type_a = 0, type_b = 0;
  /* XXX: Can the hardcoding of decay channels be avoided? */
  if ((*types)[(*map_type)[*particle_id]].spin() % 2 == 0) {
    /* meson resonance decays into pions */
    if (charge == 0) {
      type_a = 211;
      type_b = -211;
    } else if (charge == 1) {
      type_a = 211;
      type_b = 111;
    } else if (charge == -1) {
      type_a = -211;
      type_b = 111;
    }
  } else if ((*types)[(*map_type)[*particle_id]].pdgcode() > 0) {
    /* Baryon resonance decays into pion and baryon */
    if (charge == 0) {
      type_a = 2212;
      type_b = -211;
    } else if (charge == 1) {
      type_a = 2112;
      type_b = 211;
    } else if (charge == -1) {
      type_a = 2112;
      type_b = -211;
    } else if (charge == 2) {
      type_a = 2212;
      type_b = 211;
    }
  } else {
    /* Antibaryon resonance decays into pion and antibaryon */
    if (charge == 0) {
      type_a = -2212;
      type_b = 211;
    } else if (charge == 1) {
      type_a = -2112;
      type_b = 211;
    } else if (charge == -1) {
      type_a = -2112;
      type_b = -211;
    } else if (charge == -2) {
      type_a = -2212;
      type_b = -211;
    }
  }

  /* Find the desired particle types */
  bool not_found_a = true, not_found_b = true;
  size_t type_index = 0;
  while ( (not_found_a || not_found_b) && type_index < (*types).size() ) {
    if ((*types)[type_index].pdgcode() == type_a) {
      printd("Found particle %i.\n", type_a);
      (*map_type)[new_id_a] = type_index;
      not_found_a = false;
    }
    if ((*types)[type_index].pdgcode() == type_b) {
      printd("Found particle %i.\n", type_b);
      (*map_type)[new_id_b] = type_index;
      not_found_b = false;
    }
    type_index++;
  }

  const double total_energy = ((*particles)[*particle_id]).momentum().x0();
  double mass_a = (*types)[(*map_type)[new_id_a]].mass(),
    mass_b = (*types)[(*map_type)[new_id_b]].mass();
  double energy_a = (total_energy * total_energy
                     + mass_a * mass_a - mass_b * mass_b)
                    / (2.0 * total_energy);

  double momentum_radial = sqrt(energy_a * energy_a - mass_a * mass_a);
  /* phi in the range from [0, 2 * pi) */
  double phi = 2.0 * M_PI * drand48();
  /* cos(theta) in the range from [-1.0, 1.0) */
  double cos_theta = -1.0 + 2.0 * drand48();
  double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
  if (energy_a  < mass_a || abs(cos_theta) > 1) {
    printf("Particle %lu radial momenta %g phi %g cos_theta %g\n", new_id_a,
         momentum_radial, phi, cos_theta);
    printf("Etot: %g m_a: %g m_b %g E_a: %g", total_energy, mass_a, mass_b,
           energy_a);
  }
  (*particles)[new_id_a].set_momentum(mass_a,
      momentum_radial * cos(phi) * sin_theta,
      momentum_radial * sin(phi) * sin_theta,
      momentum_radial * cos_theta);
  (*particles)[new_id_b].set_momentum(mass_b,
    - (*particles)[new_id_a].momentum().x1(),
    - (*particles)[new_id_a].momentum().x2(),
    - (*particles)[new_id_a].momentum().x3());

  /* Both decay products begin from the same point */
  FourVector decay_point = (*particles)[*particle_id].position();
  (*particles)[new_id_a].set_position(decay_point);
  (*particles)[new_id_b].set_position(decay_point);

  /* No collision yet */
  (*particles)[new_id_a].set_collision(-1, 0, -1);
  (*particles)[new_id_b].set_collision(-1, 0, -1);

  printd("Created %s and %s with IDs %lu and %lu \n",
  (*types)[(*map_type)[new_id_a]].name().c_str(),
  (*types)[(*map_type)[new_id_b]].name().c_str(), new_id_a, new_id_b);

  return new_id_a;
}

/* 2->1 resonance formation process */
size_t resonance_formation(std::map<int, ParticleData> *particles,
  std::vector<ParticleType> *types, std::map<int, int> *map_type,
  int *particle_id, int *other_id, int type_resonance, int *id_max) {
  /* Add a new particle */
  size_t new_id = *id_max + 1;
  (*id_max)++;
  {
  ParticleData new_particle;
  (*particles)[new_id] = new_particle;
  (*particles)[new_id].set_id(new_id);
  }
  /* Find the desired resonance */
  bool not_found = true;
  size_t type_index = 0;
  while (not_found && type_index < types->size()) {
    if ((*types)[type_index].pdgcode() == type_resonance) {
      printd("Found resonance %i (%s).\n", type_resonance,
             (*types)[type_index].name().c_str());
      printd("Parent particles: %s %s \n",
             (*types)[(*map_type)[*particle_id]].name().c_str(),
             (*types)[(*map_type)[*other_id]].name().c_str());
      (*map_type)[new_id] = type_index;
      not_found = false;
    }
    type_index++;
  }

  /* Center-of-momentum frame of initial particles
   * is the rest frame of the resonance
   */
  const double energy = (*particles)[*particle_id].momentum().x0()
    + (*particles)[*other_id].momentum().x0();
  /* We use fourvector to set 4-momentum, as setting it
   * with doubles requires that particle is on
   * mass shell, which is not generally true for resonances
   */
  FourVector resonance_momentum(energy, 0.0, 0.0, 0.0);
  (*particles)[new_id].set_momentum(resonance_momentum);

  printd("Momentum of the new particle: %g %g %g %g \n",
    (*particles)[new_id].momentum().x0(), (*particles)[new_id].momentum().x1(),
    (*particles)[new_id].momentum().x2(), (*particles)[new_id].momentum().x3());

  /* The real position should be between parents in the computational frame! */
  (*particles)[new_id].set_position(1.0, 0.0, 0.0, 0.0);

  /* No collision yet */
  (*particles)[new_id].set_collision(-1, 0, -1);

  printd("Created %s with ID %i \n",
         (*types)[(*map_type)[new_id]].name().c_str(),
         (*particles)[new_id].id());

  return new_id;
}
