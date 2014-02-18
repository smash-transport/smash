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
#include "include/configuration.h"
#include "include/experimentparameters.h"
#include "include/outputroutines.h"

NucleusModus::NucleusModus(Configuration &config)
    : sqrt_s_NN_(config.take({"Nucleus", "SQRTSNN"})) {
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
  float mass_projec = projectile_.mass();
  float mass_target = target_.mass();
  float mass_1 = particles->particle_type(pdg_sNN_1_).mass();
  float mass_2 = particles->particle_type(pdg_sNN_2_).mass();
  float total_mandelstam_s = (sqrt_s_NN_ - mass_1*mass_1 - mass_2*mass_2)
                             * mass_projec*mass_target
                             / (mass_1*mass_2)
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
