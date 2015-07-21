/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */


#include "include/decayactiondilepton.h"

#include "include/angles.h"
#include "include/kinematics.h"

namespace Smash {

DecayActionDilepton::DecayActionDilepton(const ParticleData &p,
                                         float time,
                                         float shining_weight,
                                         float dilepton_mass)
  : DecayAction({p}, time), shining_weight_(shining_weight), dilepton_mass_(dilepton_mass) {}  // #CleanUp


void DecayActionDilepton::one_to_three() {

// find the non lepton particle position
  int non_lepton_position = -1;
  for (int i=0; i<3; ++i) {
    if (!(outgoing_particles_[i].type().is_lepton())) {
      non_lepton_position = i;
      break;
    }
  }

  if (non_lepton_position == -1) {
    throw std::runtime_error("Error in DecayActionDilepton::one_to_three.");
  }

  ParticleData &nl = outgoing_particles_[non_lepton_position];
  ParticleData &l1 = outgoing_particles_[(non_lepton_position+1)%3];
  ParticleData &l2 = outgoing_particles_[(non_lepton_position+2)%3];

  const ParticleType &l1_type = l1.type();
  const ParticleType &l2_type = l2.type();
  const ParticleType &nl_type = nl.type();

  const double mass_l1 = l1_type.mass();
  const double mass_l2 = l2_type.mass();
  const double mass_nl = nl_type.mass();

  double cms_energy = sqrt_s();


// perform non_dilepton/virtual photon decay (checks from onetotwo missing)

  double momentum_radial = pCM(cms_energy, double(dilepton_mass_), mass_nl);

  const double p0_v_pho = std::sqrt(dilepton_mass_*dilepton_mass_ + momentum_radial * momentum_radial);

  /* Here we assume an isotropic angular distribution. */
  Angles phitheta;
  phitheta.distribute_isotropically();

  FourVector virtual_pho_4mom(p0_v_pho,  phitheta.threevec() * momentum_radial);
  nl.set_4momentum(mass_nl, -phitheta.threevec() * momentum_radial);

// perform virtual photon -> dilepton decay

  cms_energy = virtual_pho_4mom.abs();

  momentum_radial = pCM(cms_energy, mass_l1, mass_l2);

  /* Here we assume an isotropic angular distribution. */
  phitheta.distribute_isotropically();

  l1.set_4momentum(mass_l1,  phitheta.threevec() * momentum_radial);
  l2.set_4momentum(mass_l2, -phitheta.threevec() * momentum_radial);

// Boost Dileptons back in parent particle rest frame

  ThreeVector velocity_CM = virtual_pho_4mom.velocity();
  l1.boost_momentum(-velocity_CM);
  l2.boost_momentum(-velocity_CM);

}

}  // namespace Smash
