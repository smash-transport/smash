#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <streambuf>

#include "example.h"

#include "smash/particles.h"
#include "smash/random.h"
#include "smash/decaymodes.h"
#include "smash/angles.h"

using namespace smash;

static Integrator integrate;

void initialize_random_number_generator() {
  // Seed with a truly random 63-bit value, if possible
  std::random_device rd;
  static_assert(std::is_same<decltype(rd()), uint32_t>::value,
               "random_device is assumed to generate uint32_t");
  uint64_t unsigned_seed = (static_cast<uint64_t>(rd()) << 32) |
                            static_cast<uint64_t>(rd());
  // Discard the highest bit to make sure it fits into a positive int64_t
  int64_t seed = static_cast<int64_t>(unsigned_seed >> 1);
  random::set_seed(seed);
}

int main(int argc, char *argv[]) {
  const int example_number = (argc > 1) ? std::stoi(argv[1]) : 1000;

  if (example_number > 0) {
    std::cout << "\nExample 1\n---------\n" << std::endl;
    std::cout << "Using SMASH wrapper of random number generator" << std::endl;
    initialize_random_number_generator();

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
    std::ifstream particles_input_file("../particles.txt");
    std::stringstream buffer;
    if (particles_input_file) {
      buffer << particles_input_file.rdbuf();
      ParticleType::create_type_list(buffer.str());
    } else {
      std::cout << "File with SMASH particle list not found." << std::endl;
      return 1;
    }
    std::ifstream decaymodes_input_file("../decaymodes.txt");
    if (decaymodes_input_file) {
      buffer.clear();
      buffer.str(std::string());
      buffer << decaymodes_input_file.rdbuf();
      DecayModes::load_decaymodes(buffer.str());
      ParticleType::check_consistency();
    } else {
      std::cout << "File with SMASH decaymodes not found." << std::endl;
      return 1;
    }
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
  
}
