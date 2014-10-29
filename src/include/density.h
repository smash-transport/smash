/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#ifndef SRC_INCLUDE_DENSITY_H_
#define SRC_INCLUDE_DENSITY_H_

#include <vector>

#include "fourvector.h"
#include "threevector.h"
#include "particledata.h"
#include "pdgcode.h"
#include "forwarddeclarations.h"

namespace Smash {

  enum Density_type {baryon, proton, neutron};

  bool particle_in_denstype(const PdgCode pdg, Density_type dens_type);

  /** Calculates 4-current in the computational frame.
   *  \f[j^{\mu} = (\sqrt{2\pi} \sigma )^{-3} \sum_{i=1}^N C_i u^{\mu}_i
   *  exp \left(- \frac{(\vec r -\vec r_i + \frac{\gamma_i^2}{1 + \gamma_i}
   *  \vec \beta_i (\vec \beta_i, \vec r - \vec r_i))^2}{2\sigma^2} \right)\f]
   *  \f[ \rho^{Eckart} = \sqrt{j^{\mu} j_{\mu}} \f]
   *  Here \f$ C_i \f$ is a corresponding value of "charge". If baryon
   *  current option is selected then \f$ C_i \f$ is 1 for baryons,
   *  -1 for antibaryons and 0 otherwize. For proton/neutron
   *  current \f$ C_i = 1\f$ for proton/neutron and 0 otherwize.
   *
   * \param[in] r Arbitrary space point where 4-current is calculated
   * \param[in] plist List of all particles to be used in \f$j^{\mu}\f$
   *            calculation. If the distance between particle and calculation
   *            point r, \f$ |r-r_i| > 6 \sigma \f$ then particle input
   *            to density will be ignored.
   * \param[in] gs_sigma Width of the gaussian (\f$ \sigma \f$),
   *  which represents particle Wigner density.
   * \param[in] type of four-currect to be calculated: baryon, proton or
   *            neutron options are currently available
   */
  FourVector four_current(ThreeVector r, const ParticleList &plist,
                          double gs_sigma, Density_type dens_type);

  std::pair<double, ThreeVector> rho_eckart_gradient(ThreeVector r,
                                              const ParticleList &plist,
                                  double gs_sigma, Density_type dens_type);

  

}  // namespace Smash

#endif  // SRC_INCLUDE_DENSITY_H_
