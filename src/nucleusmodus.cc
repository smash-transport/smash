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
      projectile_.push_back(atoi(value));
      match = true;
    }
    if (strcmp(key, "TARGET") == 0) {
      target_.push_back(atoi(value));
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
  for (std::vector<int>::iterator p = projectile_.begin();
                       p != projectile_.end(); p++) {
    printf("Particle in projectile: %d\n", *p);
  }
  for (std::vector<int>::iterator t = target_.begin();
                       t != target_.end(); t++) {
    printf("Particle in target: %d\n", *t);
  }
  printf("Center-of-mass energy %10.3f GeV\n", sqrts_);
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

                                       
  // make a temporary list to do identical things to projectile and
  // target.
  std::vector<int>* nuclei[2] = { &projectile_, &target_ };
  for (int pt = 0; pt < 1; pt++) {
    int direction = (pt == 0) ? 1 : -1;
    std::vector<int>* nucleus = nuclei[pt];
    // get nuclear radius which should depend on the number of particles
    // inside.
    float radius = nuclear_radius(nucleus->size());
    for (std::vector<int>::iterator i = nucleus->begin();
                         i != nucleus->end(); i++) {
      float position_of_nucleon = woods_saxon(radius);
      Angles phitheta;
      phitheta.distribute_isotropically();
      ParticleData new_particle();
      // position of particle is within a woods-saxon distribution
      // sphere centered around (0,0,-D*I*R), where I is a parameter and R
      // the radius. D is +- 1, depending on which nucleus we are at.
      new_particle.set_position(0.0, position_of_nucleon*phitheta.x()
                                   , position_of_nucleon*phitheta.y()
                                   , position_of_nucleon*phitheta.z()
                                   - direction*initial_displacement_*radius );
      // TODO(baeuchle) particles of different mass have a different momentum!
      float p_lab = sqrts_;
      new_particle.set_momentum(0.0, 0.0, 0.0, direction*p_lab);
    }
  }

  for (auto i = particles->begin(); i != particles->end(); i++) {
    float mass = particles->type(i->first).mass();
    printf("id %d pdgcode %d mass %f\n", i->first, i->second.pdgcode(), mass);
    /* velocity of particles */
    double cms_gamma = sqrts_ / mass;
    double cms_beta = sqrt(sqrts_ * sqrts_ - mass * mass / sqrts_);
    // Sample impact parameter
    double impact_parameter = drand48() * 5.0;
    if (i->first == 0) {
      i->second.set_position(1.0, impact_parameter, 0.0, -1.0);
      i->second.set_momentum(mass, 0.0, 0.0, cms_gamma * cms_beta * mass);
    } else if (i->first == 1) {
      i->second.set_position(1.0, 0.0, 0.0, 1.0);
      i->second.set_momentum(mass, 0.0, 0.0, -cms_gamma * cms_beta * mass);
    }
  }
}

float NucleusModus::nuclear_radius(const int &number_of_particles) {
  return 1.2*pow(number_of_particles,1./3.);
}

/* returns a length r with a relative probability of $$\frac{dN}{dr} =
 * \frac{r^2}{\exp\left(\frac{r-R}{d}\right) + 1}$$ where $d$ is a class
 * parameter and $R$ is the function parameter ''radius''. */
float NucleusModus::woods_saxon(const float &radius) {
  float radius_scaled = radius/steepness_;
  float prob_range1 = 1.0;
  float prob_range2 = 3. / radius_scaled;
  float prob_range3 = 2. * prob_range2 / radius_scaled;
  float prob_range4 = 1. * prob_range3 / radius_scaled;
  float all_ranges = prob_range1 + prob_range2 + prob_range3 + prob_range4;
  float t;
  do {
    float which_range = drand48() * all_ranges - prob_range1;
    if (which_range < 0.0) {
      t = radius_scaled * (pow(drand48(),1./3.) - 1.) ;
    }
    else {
      t = -log(drand48());
      if (which_range >= prob_range2) {
        t += -log(drand48());
        if (which_range >= prob_range2 + prob_range3) {
          t += -log(drand48());
        }
      }
    }
  } while (drand48() > 1./(1. + exp(-fabs(t)) ) );
  float position_scaled = t + radius_scaled;
  float position = position_scaled * steepness_;
  return position;
}
