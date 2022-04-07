#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <streambuf>

#include "smash/particles.h"
#include "smash/random.h"
#include "smash/decaymodes.h"
#include "smash/angles.h"
#include "smash/setup_particles_decaymodes.h"

using namespace smash;

static Integrator integrate;

int main(int argc, char *argv[]) {
  const int example_number = (argc > 1) ? std::stoi(argv[1]) : 1000;

  if (example_number > 0) {
    std::cout << "\nExample 1\n---------\n" << std::endl;
    std::cout << "Using SMASH wrapper of random number generator" << std::endl;
    random::set_seed(random::generate_63bit_seed());

    Angles phitheta;
    phitheta.distribute_isotropically();
    double r_length = 5.0;
    ThreeVector r(phitheta.threevec() * r_length);
    std::cout << "A random vector of length " << r_length
              << ": " << r << std::endl;
    const double mean_exp = 2.0;
    std::cout << "Drawing random number from exponential distribution of mean "
              << mean_exp << ": "
              << random::exponential(1.0 / mean_exp) << std::endl;
  }
  if (example_number > 1) {
    std::cout << "\nExample 2\n---------\n" << std::endl;
    std::cout << "Loading SMASH particle types and decay modes" << std::endl;
    smash::initialize_default_particles_and_decaymodes();
    std::cout << "Print all strange mesons lighter than 1 GeV" << std::endl;
    for (const ParticleType &ptype : ParticleType::list_all()) {
      if (ptype.is_meson() && ptype.strangeness() != 0 && ptype.mass() < 1.0) {
        std::cout << ptype << std::endl;
      }
    }
  }
  if (example_number > 2) {
    std::cout << "\nExample 3\n---------\n" << std::endl;
    std::cout << "List particles with pole width < 0.1 GeV and mass < 1.5 GeV,"
              << " which decay into Lambda" << std::endl;
    for (const ParticleType &ptype : ParticleType::list_all()) {
      if (ptype.width_at_pole() > 0.1 || ptype.mass() > 1.5) {
        continue;
      }
      const auto &modes = ptype.decay_modes().decay_mode_list();
      for (const auto &decay_branch : modes) {
        bool decay_mode_has_rho = false;
        for (const ParticleTypePtr decay_into : decay_branch->particle_types()) {
          if (decay_into->pdgcode() == pdg::Lambda) {
            decay_mode_has_rho = true;
          }
        }
        if (decay_mode_has_rho) {
          std::cout << ptype.name() << "->";
          for (const ParticleTypePtr decay_product :
                     decay_branch->particle_types()) {
            std::cout << decay_product->name();
          }
          std::cout << " " << decay_branch->weight() << std::endl;

        }
      }
    }
  }

}
