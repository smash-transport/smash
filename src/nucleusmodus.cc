/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#include <cmath>
#include <cstdlib>
#include <list>

#include "include/nucleusmodus.h"
#include "include/angles.h"
#include "include/experimentparameters.h"
#include "include/outputroutines.h"
#include "include/parameters.h"

void NucleusModus::assign_params(std::list<Parameters> *configuration) {
  bool match = false;
  std::list<Parameters>::iterator i = configuration->begin();
  while (i != configuration->end()) {
    char *key = i->key();
    char *value = i->value();
    printd("%s %s\n", key, value);
    /* integer values */
    if (strcmp(key, "PROJECTILE") == 0) {
      // projectile_.add_particle(atoi(value));
      match = true;
    }
    if (strcmp(key, "TARGET") == 0) {
      // target_.push_back(atoi(value));
      match = true;
    }
    /* float values */
    if (strcmp(key, "SQRTS") == 0) {
      sqrts_ = (fabs(atof(value)));
      match = true;
    }
    /* remove processed entry */
    if (match) {
      i = configuration->erase(i);
      match = false;
    } else {
      ++i;
    }
  }
}

/* print_startup - console output on startup of box specific parameters */
void NucleusModus::print_startup() {
//  for (std::vector<int>::iterator p = projectile_.begin();
//                       p != projectile_.end(); p++) {
//    printf("Particle in projectile: %d\n", *p);
//  }
//  for (std::vector<int>::iterator t = target_.begin();
//                       t != target_.end(); t++) {
//    printf("Particle in target: %d\n", *t);
//  }
//  printf("Center-of-mass energy %10.3f GeV\n", sqrts_);
}

/* initial_conditions - sets particle data for @particles */
void NucleusModus::initial_conditions(Particles *particles,
                                       const ExperimentParameters &) {
  // PHYSICS:
  //    First, calculate the masses of both nuclei MA and MB.
  //    Second: gamma*beta = sqrt( s-(MA+MB)**2 )/sqrt(4*MA*MB)
  //    Third: Get momenta from p_i = +-gamma*beta*mass_i
  //
  //    That is if a sqrt(s)_total is given and we want calculations to
  //    happen in frame of equal velocity.
  //
  //    if sqrt(s)_(ij) is given (i and j being the indices of two
  //    particles), then
  //    s_tot = (s_ij - mi**2 - mj**2) * MA*MB/(mi*mj) + MA**2 + MB**2.
  //
  //    One could also think of defining sqrt(s)_NN as the collision of
  //    two particles whose mass is mA = MA/NA, where NA is the number
  //    of particles in nucleus A, and correspondingly for mB.
  //
  //    Momenta/positions
  //    First, initialize nuclei at (0,0,0) at rest
  //    Second, boost nuclei
  //    Third, shift them so that they barely touch each other
  //    Fourth, set the time when they /will/ touch to 0.0.
  //
  // NUMERICS:
  //
  //    Maybe there should be a new class "Nucleus" that carries the
  //    parameters of that nucleus like particle list (not the particles
  //    themselves, but the information "5 protons, 7 Lambdas, 2
  //    neutrons"), mass and initial displacement (as well as initial
  //    velocity). This would only be needed here, though, so it is not
  //    clear if there is a benefit from a new class.
  float sqrt_s_NN = 24.f;
  float mass_projec = projectile_.get_mass();
  float mass_target = target_.get_mass();
  float mass_proton = particles->particle_type(2212).mass();
  float total_mandelstam_s = (sqrt_s_NN - 2.0*mass_proton*mass_proton)
                             * mass_projec*mass_target
                             / (mass_proton*mass_proton)
                           + mass_projec*mass_projec + mass_target*mass_target;
  float velocity_squared = (total_mandelstam_s - (mass_projec+mass_target))
                         / (total_mandelstam_s - (mass_projec-mass_target));
  // populate the nuclei with appropriately distributed nucleons
  projectile_.arrange_nucleons();
  target_.arrange_nucleons();
  // boost the nuclei to the appropriate velocity (target is in opposite
  // direction!
  projectile_.boost(velocity_squared);
  target_.boost(-velocity_squared);
}
