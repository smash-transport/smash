/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"

#include "../include/particledata.h"
#include "../include/decaymodes.h"
#include "../include/action.h"

namespace particles_txt {
#include <particles.txt.h>
}

namespace decaymodes_txt {
#include <decaymodes.txt.h>
}

using namespace Smash;

TEST(init_particle_types) {
  ParticleType::create_type_list(particles_txt::data);
  DecayModes::load_decaymodes(decaymodes_txt::data);
}


static ScatterAction* set_up_action(const ParticleData &proj,
                                    const ParticleData &targ,
                                    CollisionBranchList &proc_list) {
  ScatterAction *act = new ScatterActionBaryonBaryon(proj, targ, 0.);
  proc_list = act->two_to_two_cross_sections();
//   act->add_processes(proc_list);

  std::printf("%s+ %s, sqrt(s) = %f GeV, sigma = %f mb, %lu Channels \n",
              proj.type().name().c_str(), targ.type().name().c_str(),
              act->sqrt_s(), act->cross_section(), proc_list.size());

  for (auto &proc: proc_list) {
    std::printf("-> %s %s (%f mb) \n",
                proc->particle_list()[0].pdgcode().string().c_str(),
                proc->particle_list()[1].pdgcode().string().c_str(),
                proc->weight());
  };

  return act;
}


TEST(NN_NDelta) {
  // test isospin symmetry for N N -> N Delta
  ParticleTypePtr proton  = &ParticleType::find(0x2212);
  ParticleTypePtr neutron = &ParticleType::find(0x2112);

  ParticleData p1 = ParticleData{*proton, 1};
  ParticleData p2 = ParticleData{*proton, 2};
  ParticleData n  = ParticleData{*neutron,3};

  const double ptot=0.7;
  p1.set_4momentum(proton->mass(), 0., 0., ptot);
  p2.set_4momentum(proton->mass(), 0., 0., -ptot);

  ScatterAction *act_pp, *act_pn, *act_np;
  CollisionBranchList proc_list_pp, proc_list_pn, proc_list_np;

  // p p -> N Delta
  act_pp = set_up_action (p1, p2, proc_list_pp);
  // p n -> N Delta
  n.set_4momentum(proton->mass(), 0., 0., -ptot);
  act_pn = set_up_action (p1, n, proc_list_pn);
  // n p -> N Delta
  n.set_4momentum(proton->mass(), 0., 0., ptot);
  act_np = set_up_action (n, p2, proc_list_np);

  VERIFY(proc_list_pp.size()==2);
  VERIFY(proc_list_pn.size()==2);
  VERIFY(proc_list_np.size()==2);

  // check isospin ratios
  FUZZY_COMPARE(3*proc_list_pp[0]->weight(),proc_list_pp[1]->weight());  // ratio 1:3
  FUZZY_COMPARE(proc_list_pn[0]->weight(),proc_list_pn[1]->weight());    // ratio 1:1
  FUZZY_COMPARE(proc_list_pn[0]->weight(),proc_list_np[1]->weight());    // ratio 1:1

  FUZZY_COMPARE(proc_list_pp[0]->weight(),proc_list_pn[0]->weight());    // ratio 1:1
  FUZZY_COMPARE(proc_list_pp[0]->weight(),proc_list_np[0]->weight());    // ratio 1:1

//   FUZZY_COMPARE(act_pp->weight(),2*act_pn->weight());                  // ratio 2:1

  delete act_pp;
  delete act_pn;
  delete act_np;
}


TEST(NDelta_NN) {
  // test isospin symmetry for N Delta -> N N
  ParticleTypePtr proton   = &ParticleType::find(0x2212);
  ParticleTypePtr neutron  = &ParticleType::find(0x2112);
  ParticleTypePtr Delta_pp = &ParticleType::find(0x2224);
  ParticleTypePtr Delta_p  = &ParticleType::find(0x2214);
//   ParticleTypePtr Delta_z  = &ParticleType::find(0x2114);
//   ParticleTypePtr Delta_m  = &ParticleType::find(0x1114);

  ParticleData Dp  = ParticleData{*Delta_p, 1};
  ParticleData Dpp = ParticleData{*Delta_pp,2};
  ParticleData p   = ParticleData{*proton,  3};
  ParticleData n   = ParticleData{*neutron, 4};

  const double ptot=0.7;
  Dp.set_4momentum (Delta_p->mass(),  0., 0., ptot);
  Dpp.set_4momentum(Delta_pp->mass(), 0., 0., ptot);
  p.set_4momentum  (proton->mass(),   0., 0., -ptot);
  n.set_4momentum  (proton->mass(),   0., 0., -ptot);

  ScatterAction *act_Dp, *act_Dn, *act_DDn;
  CollisionBranchList proc_list_Dp, proc_list_Dn, proc_list_DDn;

  // Delta+ p -> N N
  act_Dp = set_up_action (Dp, p, proc_list_Dp);
  // Delta+ n -> N N
  act_Dn = set_up_action (Dp, n, proc_list_Dn);
  // Delta++ n -> N N
  act_DDn = set_up_action (Dpp, n, proc_list_DDn);

  VERIFY(proc_list_Dp.size()==1);
  VERIFY(proc_list_Dn.size()==2);
  VERIFY(proc_list_DDn.size()==1);

  // check isospin ratios
  FUZZY_COMPARE(proc_list_Dn[0]->weight(),proc_list_Dn[1]->weight());     // ratio 1:1
  FUZZY_COMPARE(proc_list_Dp[0]->weight(),proc_list_Dn[0]->weight());     // ratio 1:1
  FUZZY_COMPARE(3*proc_list_Dp[0]->weight(),proc_list_DDn[0]->weight());  // ratio 1:3

//   FUZZY_COMPARE(2*act_Dp->weight(),act_Dn->weight());                  // ratio 1:2

  delete act_Dp;
  delete act_Dn;
  delete act_DDn;
}
