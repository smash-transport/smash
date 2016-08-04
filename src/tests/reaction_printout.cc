/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "setup.h"

#include "../include/isoparticletype.h"
#include "../include/particletype.h"
#include "../include/particledata.h"
#include "../include/scatteraction.h"
#include "../include/scatteractionsfinder.h"

using namespace Smash;

static ScatterActionPtr construct_scatter_action(const ParticleData &a,
                                                 const ParticleData &b) {
  constexpr float elastic_parameter = -1.0f;
  constexpr int ntest = 1;
  ScatterActionsFinder finder(elastic_parameter, ntest);
  ScatterActionPtr act = finder.construct_scatter(a, b);
  return act;
}

static void remove_substr(std::string& s, const std::string& p) { 
  std::string::size_type n = p.length();
  for (std::string::size_type i = s.find(p); i != std::string::npos; i = s.find(p)) {
    s.erase(i, n);
  }
}

static std::string isoclean(std::string s) {
  remove_substr(s, "⁺");
  remove_substr(s, "⁻");
  remove_substr(s, "⁰");
  return s;
}

TEST(init_particle_types) {
  Test::create_actual_particletypes();
}

TEST(init_decaymodes) {
  Test::create_actual_decaymodes();
}

TEST(printout_possible_channels) {
  const size_t N_isotypes = IsoParticleType::list_all().size();
  const size_t N_pairs = N_isotypes * (N_isotypes - 1) / 2;
  // We do not consider decays here, only 2->n processes, where n > 1, 
  constexpr bool two_to_one = false;
  // We do not set all elastic cross-sections to fixed value
  constexpr float elastic_parameter = -1.0f;
  constexpr bool two_to_two = true, strings_switch = true;
  std::cout << N_isotypes << " iso-particle types." << std::endl;
  std::cout << "They can make " << N_pairs << " pairs." << std::endl;
  
  for (const IsoParticleType &A_isotype : IsoParticleType::list_all()) {
    for (const IsoParticleType &B_isotype : IsoParticleType::list_all()) {
      if (&A_isotype > &B_isotype) {
        continue;
      }
      bool any_nonzero_cs = false;
      std::vector<std::string> r_list;
      for (const ParticleTypePtr A_type : A_isotype.get_states()) {
        for (const ParticleTypePtr B_type : B_isotype.get_states()) {
          if (A_type > B_type) {
            continue;
          }
          ParticleData A(*A_type), B(*B_type);
          // Momentum should be enough to get non-zero string cross-section
          A.set_4momentum(A.pole_mass(), 3.0, 0.0, 0.0);
          B.set_4momentum(B.pole_mass(), -3.0, 0.0, 0.0);
          ScatterActionPtr act = construct_scatter_action(A, B);
          act->add_all_processes(elastic_parameter, two_to_one,
                                 two_to_two, strings_switch);
          const float total_cs = act->cross_section();
          if (total_cs <= 0.0) {
            continue;
          }
          any_nonzero_cs = true;
          for (const auto& channel : act->collision_channels()) {
            std::string r;
            if (channel->get_type() == ProcessType::String) {
              r =  A_type->name() + B_type->name()
                   + std::string("->strings");
            } else {
              std::string r_type =
                (channel->get_type() == ProcessType::Elastic) ?
                std::string(" (el)") :
                      (channel->get_type() == ProcessType::TwoToTwo) ?
                      std::string(" (inel)") :
                           std::string(" (?)");
              r = A_type->name() + B_type->name()
                    + std::string("->")
                    + channel->particle_types()[0]->name()
                    + channel->particle_types()[1]->name()
                    + r_type;
            }
            r_list.push_back(isoclean(r));
          }  
        }
      }
      std::sort(r_list.begin(), r_list.end());
      r_list.erase(std::unique(r_list.begin(), r_list.end()), r_list.end() );
      if (any_nonzero_cs) {
        for (auto r : r_list) {
          std::cout << r << ", ";
        }
        std::cout << std::endl;
      }
    }
  }
}
