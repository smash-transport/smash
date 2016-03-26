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
                                         float shining_weight)
  : DecayAction({p}, time), shining_weight_(shining_weight) {}


void DecayActionDilepton::one_to_three() {
  // find the non-lepton particle position
  int non_lepton_position = -1;
  for (int i = 0; i < 3; ++i) {
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

  const double cms_energy = sqrt_s();

  // randomly select a dilepton mass
  const double dil_mass = Random::uniform(mass_l1 + mass_l2,
                                          cms_energy - mass_nl);
  const float delta_m = cms_energy - mass_nl - mass_l1 - mass_l2;

  const float diff_width = ThreeBodyDecayDilepton::diff_width(
                                      cms_energy, dil_mass, mass_nl,
                                      incoming_particles_[0].type().pdgcode());

  // correct shining weight for differential width at particular dilepton mass
  shining_weight_ *= delta_m * diff_width / decay_channels_[0]->weight();


  // perform decay into non-lepton and virtual photon (dilepton)
  const double dil_mom = pCM(cms_energy, dil_mass, mass_nl);

  /* Here we assume an isotropic angular distribution. */
  Angles phitheta;
  phitheta.distribute_isotropically();

  FourVector dil_4mom(std::sqrt(dil_mass*dil_mass + dil_mom*dil_mom),
                      phitheta.threevec() * dil_mom);
  nl.set_4momentum(mass_nl, -phitheta.threevec() * dil_mom);


  // perform decay of virtual photon into two leptons
  const double mom_lep = pCM(dil_mass, mass_l1, mass_l2);

  /* Here we assume an isotropic angular distribution. */
  phitheta.distribute_isotropically();

  l1.set_4momentum(mass_l1,  phitheta.threevec() * mom_lep);
  l2.set_4momentum(mass_l2, -phitheta.threevec() * mom_lep);

  // Boost Dileptons back in parent particle rest frame
  ThreeVector velocity_CM = dil_4mom.velocity();
  l1.boost_momentum(-velocity_CM);
  l2.boost_momentum(-velocity_CM);
}


}  // namespace Smash
