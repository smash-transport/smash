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

double resonance_cross_section(ParticleData *particle1, ParticleData *particle2,
  ParticleType *type_particle1, ParticleType *type_particle2,
  std::vector<ParticleType> *type_list) {
  const int charge1 = (*type_particle1).charge(),
    charge2 = (*type_particle2).charge();

  /* Resonances do not form resonances */
  if (type_particle1->width() > 0.0 || type_particle2->width() > 0.0)
    return 0.0;

  /* Total charge defines the type of resonance */
  /* We have no resonances with charge > 1 */
  if (abs(charge1 + charge2) > 1)
    return 0.0;

  int type_resonance;
  if (charge1 + charge2 == 1)
    type_resonance = 213;
  else if (charge1 + charge2 == -1)
    type_resonance = -213;
  else
    type_resonance = 113;

  /* Find the width and mass of the desired resonance */
  float resonance_width = -1.0, resonance_mass = 0.0;
  size_t type_index = 0;
  while (resonance_width < 0 && type_index < (*type_list).size()) {
    if ((*type_list)[type_index].pdgcode() == type_resonance) {
      resonance_width = (*type_list)[type_index].width();
      resonance_mass = (*type_list)[type_index].mass();
      printd("Found resonance %i with mass %f and width %f.\n",
             type_resonance, resonance_mass, resonance_width);
      printd("Original particles: %s %s Charges: %i %i \n",
             (*type_particle1).name().c_str(), (*type_particle2).name().c_str(),
             (*type_particle1).charge(), (*type_particle2).charge());
    }
    type_index++;
  }

  /* If there was no such resonance in the list, return 0 */
  if (resonance_width < 0.0)
    return 0.0;

  /* XXX: Calculate isospin Clebsch-Gordan coefficient */
  const double clebsch_gordan_isospin = 10;

  /* If Clebsch-Gordan coefficient is zero, don't bother with the rest */
  if (clebsch_gordan_isospin < really_small)
    return 0.0;

  /* Calculate spin factor */
  const double spinfactor = (2 * (*type_list)[type_index].spin() + 1)
    / ((2 * type_particle1->spin() + 1) * (2 * type_particle2->spin() + 1));

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


  /* Calculate resonance production cross section
   * using the Breit-Wigner distribution as probability amplitude
   */
  return clebsch_gordan_isospin * spinfactor * symmetryfactor
         * 4.0 * M_PI / cm_momentum_squared
         * breit_wigner(mandelstam_s, resonance_mass, resonance_width);
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
  if (charge == 0) {
    type_a = 211;
    type_b = -211;
  } else if ( charge == 1 ) {
    type_a = 211;
    type_b = 111;
  } else if ( charge == -1 ) {
    type_a = -211;
    type_b = 111;
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
  printd("Particle %lu radial momenta %g phi %g cos_theta %g\n", new_id_a,
         momentum_radial, phi, cos_theta);
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
  int *particle_id, int *other_id, int *id_max) {
  /* Add a new particle */
  size_t new_id = *id_max + 1;
  (*id_max)++;
  {
  ParticleData new_particle;
  (*particles)[new_id] = new_particle;
  (*particles)[new_id].set_id(new_id);
  }
  /* Which resonance is formed */
  const int charge1 = (*types)[(*map_type)[*particle_id]].charge(),
    charge2 = (*types)[(*map_type)[*other_id]].charge();

  int type_resonance;
  if (charge1 + charge2 == 1)
    type_resonance = 213;
  else if (charge1 + charge2 == -1)
    type_resonance = -213;
  else
    type_resonance = 113;

  /* Find the desired resonance */
  bool not_found = true;
  size_t type_index = 0;
  while (not_found && type_index < (*types).size()) {
    if ((*types)[type_index].pdgcode() == type_resonance) {
      printd("Found resonance %i.\n", type_resonance);
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
