/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/decayactiondilepton.h"

#include "smash/angles.h"
#include "smash/kinematics.h"

namespace smash {

DecayActionDilepton::DecayActionDilepton(const ParticleData &p, double time,
                                         double shining_weight)
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
  ParticleData &l1 = outgoing_particles_[(non_lepton_position + 1) % 3];
  ParticleData &l2 = outgoing_particles_[(non_lepton_position + 2) % 3];

  const double mass_l1 = l1.type().mass();  // (pole) mass of first lepton
  const double mass_l2 = l2.type().mass();  // (pole) mass of second lepton
  const double mass_nl = nl.type().mass();  // (pole) mass of non-lepton

  const double cms_energy = kinetic_energy_cms();

  // randomly select a dilepton mass
  const double dil_mass =
      random::uniform(mass_l1 + mass_l2, cms_energy - mass_nl);
  const double delta_m = cms_energy - mass_nl - mass_l1 - mass_l2;

  const double diff_width = ThreeBodyDecayDilepton::diff_width(
      cms_energy, mass_l1, dil_mass, mass_nl, &nl.type(),
      &incoming_particles_[0].type());

  /* Branching factor, which corrects the shining weight for the differential
   * width at a particular dilepton mass. We do an implicit Monte-Carlo
   * integration over the dilepton mass here, and delta_m is simply the
   * integration volume. */
  branching_ = delta_m * diff_width / decay_channels_[0]->weight();

  // perform decay into non-lepton and virtual photon (dilepton)
  const double dil_mom = pCM(cms_energy, dil_mass, mass_nl);

  // Here we assume an isotropic angular distribution.
  Angles phitheta;
  phitheta.distribute_isotropically();

  FourVector dil_4mom(std::sqrt(dil_mass * dil_mass + dil_mom * dil_mom),
                      phitheta.threevec() * dil_mom);
  nl.set_4momentum(mass_nl, -phitheta.threevec() * dil_mom);

  // perform decay of virtual photon into two leptons
  const double mom_lep = pCM(dil_mass, mass_l1, mass_l2);

  // Here we assume an isotropic angular distribution.
  phitheta.distribute_isotropically();

  l1.set_4momentum(mass_l1, phitheta.threevec() * mom_lep);
  l2.set_4momentum(mass_l2, -phitheta.threevec() * mom_lep);

  // Boost Dileptons back in parent particle rest frame
  ThreeVector velocity_CM = dil_4mom.velocity();
  l1.boost_momentum(-velocity_CM);
  l2.boost_momentum(-velocity_CM);
}

}  // namespace smash
